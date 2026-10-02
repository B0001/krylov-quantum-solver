#!/usr/bin/env python3
"""
Does centering change how shot noise in <H^n> amplifies into PDS(K)?  (specs/SPEC_centered_pds_shots.md)

Circuit-level noise model: ONE classical-shadow dataset of the Hartree-Fock state (random-Pauli
snapshots, Huang-Kueng-Preskill; same estimator as ``classical_shadows.py``) feeds EVERY moment, raw
and centered. Each snapshot is the product operator rho_s = (x)_q (I + 3 s_q P_q)/2, so every
moment <O> is estimated by Tr(O rho_avg) from the same shots -- the noise in mu_1..mu_15 is
correlated exactly as it would be on hardware, not independent Gaussian per moment.

Symmetry projection (a choice of estimator, not shot post-selection, which X/Y bases do not allow): H
commutes with the (N_alpha, N_beta) projector P, so the same snapshots estimate <P H^n> / <P>, which
needs only the sector block of rho_avg. That block is 400x400 for N2 CAS(6,6); snapshots are grouped
by their (basis, sign) pattern on each half of the qubits, so 10^6 snapshots cost a few seconds.
Moments then come from the sector eigenbasis: mu_n = sum_i (lam_i - mu)^n d_i, d = diag(V^dag rho V).

Larger budgets (``scale_to``) rescale a sampled deviation d - d_exact by sqrt(S0/S). Covariance is
exact under that rescaling (cov of a mean is 1/S); higher cumulants are not.

Raw-frame PDS at K >= 6 depends on the last bits of H (cond(M) > 1e17), which depend on BLAS threads
and on PYTHONHASHSEED (upstream term order; SPEC_centered_pds_roots). ``pinned_n2_study`` fixes
both; ``paths`` re-sums H in other term orders to measure how much the raw frame moves.
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np
from scipy import sparse

from moment_expansion import pds_energy

_PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=complex)  # X,Y,Z
N2 = dict(atom="N 0 0 0; N 0 0 1.10", basis="sto3g", active_electrons=6, active_orbitals=6)
PIN = {"PYTHONHASHSEED": "0", **dict.fromkeys(
    ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"), "1")}
PATHS = (None, *range(1, 9))     # float64 paths: H's own term order, then 8 seeded re-orderings


def _popcount_states(n_bits: int, k: int) -> np.ndarray:
    return np.array([x for x in range(2 ** n_bits) if bin(x).count("1") == k])


def sector_eigensystem(mh, order_seed=None):
    """(lam, V, hf_bits, lo_states, hi_states) of H restricted to HF's (N_alpha, N_beta) sector.

    JW ordering (qiskit-nature): alpha spin-orbitals are qubits 0..n-1 (lo), beta n..2n-1 (hi).
    Sector basis order is (hi, lo) row-major, matching ``shadow_sector_rho``. ``order_seed`` sums
    H's Pauli terms in a seeded random order: the same operator on another float64 path.
    """
    n = mh.num_spatial_orbitals
    lo = _popcount_states(n, mh.num_particles[0])
    hi = _popcount_states(n, mh.num_particles[1])
    idx = (hi[:, None] << n | lo[None, :]).ravel()
    op = mh.qubit_hamiltonian
    if order_seed is not None:
        op = op[np.random.default_rng(order_seed).permutation(op.size)]
    H = op.to_matrix(sparse=True).tocsr()[idx][:, idx].toarray()
    lam, V = np.linalg.eigh(H)
    hf_index = int(np.argmax(np.abs(np.asarray(mh.hf_state().data))))
    hf_bits = np.array([(hf_index >> q) & 1 for q in range(mh.num_qubits)])
    return lam, V, hf_bits, lo, hi


def sample_hf_shadow(hf_bits: np.ndarray, n_shots: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """Random-Pauli snapshots of a computational-basis state (same output as
    ``classical_shadows.collect_classical_shadow`` on it, without the statevector loop): Z-measured
    qubits return their bit deterministically, X/Y-measured qubits a fair +/-1."""
    n = len(hf_bits)
    bases = rng.integers(0, 3, size=(n_shots, n))
    signs = np.where(bases == 2, 1 - 2 * hf_bits[None, :], rng.choice([-1, 1], size=(n_shots, n)))
    return bases, signs


def _half_block(bases, signs, qubits, states):
    """(S, m, m) block of the snapshot product operator on ``qubits`` between ``states`` (m of them)."""
    out = np.ones((bases.shape[0], len(states), len(states)), dtype=complex)
    for j, q in enumerate(qubits):
        A = 0.5 * (np.eye(2) + 3 * signs[:, q, None, None] * _PAULI[bases[:, q]])   # (S,2,2)
        b = (states >> j) & 1
        out *= A[:, b[:, None], b[None, :]]
    return out


def _patterns(bases, signs, qubits):
    """Distinct (basis, sign) patterns of the snapshots on ``qubits``: (bases_u, signs_u, inverse)."""
    k = len(qubits)
    code = (2 * bases[:, qubits] + (signs[:, qubits] < 0)) @ 6 ** np.arange(k)
    u, inv = np.unique(code, return_inverse=True)
    digits = u[:, None] // 6 ** np.arange(k) % 6
    return digits // 2, 1 - 2 * (digits % 2), inv


def shadow_sector_rho(bases, signs, lo, hi) -> np.ndarray:
    """Sector block of the mean snapshot operator, indexed [(h,l),(h',l')]. Snapshot = (hi block) x
    (lo block), so the sum over shots is Bh^T C Bl with C the (hi pattern, lo pattern) count matrix."""
    n = bases.shape[1] // 2
    bh, sh, ih = _patterns(bases, signs, range(n, 2 * n))
    bl, sl, il = _patterns(bases, signs, range(n))
    counts = sparse.csr_matrix((np.ones(len(ih)), (ih, il)), shape=(len(bh), len(bl)))
    Bh = _half_block(bh, sh, range(n), hi).reshape(len(bh), -1)
    Bl = _half_block(bl, sl, range(n), lo).reshape(len(bl), -1)
    acc = (Bh.T @ (counts @ Bl)).reshape(len(hi), len(hi), len(lo), len(lo))
    return acc.transpose(0, 2, 1, 3).reshape(len(hi) * len(lo), -1) / bases.shape[0]


def eigen_weights(V, rho) -> np.ndarray:
    """d_i = <v_i|rho|v_i>."""
    return np.real(np.sum(V.conj() * (rho @ V), axis=0))


def shadow_sector_weights(bases, signs, V, lo, hi) -> np.ndarray:
    """Eigenbasis weights from one shadow dataset, UNnormalized: sum(d) is the shadow estimate of
    <P_sector>. PDS is invariant under a common scale of all moments."""
    return eigen_weights(V, shadow_sector_rho(bases, signs, lo, hi))


def exact_weights(V, lo, hi, hf_bits) -> np.ndarray:
    n = len(hf_bits) // 2
    hf_index = sum(int(b) << q for q, b in enumerate(hf_bits))
    pos = list((hi[:, None] << n | lo[None, :]).ravel()).index(hf_index)
    return np.abs(V[pos]) ** 2


def moments(lam, d, mu: float, max_order: int) -> np.ndarray:
    """mu_n = sum_i (lam_i - mu)^n d_i, n = 0..max_order (mu=0: raw frame)."""
    return np.array([np.sum((lam - mu) ** k * d) for k in range(max_order + 1)])


def pds_or_nan(mom, K, offset):
    try:
        return pds_energy(mom, K, offset)
    except (ValueError, np.linalg.LinAlgError):
        return np.nan


def hankel_cond(mom, K):
    M = np.array([[mom[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])
    return np.linalg.cond(M)


def scale_to(d, d_exact, s0: int, s: float) -> np.ndarray:
    """Weights at budget s from normalized weights sampled at s0 (exact covariance, not higher cumulants)."""
    return d_exact + (d - d_exact) * np.sqrt(s0 / s)


SAMPLED = (10**3, 10**4, 10**5, 10**6)          # directly sampled snapshot budgets (nested prefixes)
SCALED = (10**8, 10**10, 10**12)                # covariance-rescaled from the 10^6 sample (lead only)


def noise_study(mh, trials: int = 16, max_k: int = 8, seed: int = 2000, paths=(None,)):
    """PDS(K) - E0 for raw and centered frames, per budget x trial x K, from shared-shot shadows.

    One result dict per float64 path in ``paths`` (``sector_eigensystem`` order seeds), all from the
    same snapshots. Also an INDEPENDENT-per-moment contrast with identical marginals: moment n of
    trial t is taken from trial (t + n) % trials, in each frame separately (the "measure each power
    on its own shots" model a prior attempt used). Arrays are shaped (len(budgets), trials, max_k);
    "exact" is the noiseless PDS(K) - E0.
    """
    if trials < 2 * max_k:
        raise ValueError(f"the independent contrast needs trials >= {2 * max_k} (one per moment)")
    systems = [sector_eigensystem(mh, p) for p in paths]
    hf, lo, hi = systems[0][2:]
    sampled = np.empty((len(paths), len(SAMPLED), trials, len(lo) * len(hi)))
    for t in range(trials):
        bases, signs = sample_hf_shadow(hf, SAMPLED[-1], np.random.default_rng(seed + t))
        for b, S in enumerate(SAMPLED):
            rho = shadow_sector_rho(bases[:S], signs[:S], lo, hi)
            for p, (_, V, *_) in enumerate(systems):
                d = eigen_weights(V, rho)
                sampled[p, b, t] = d / d.sum()
    return [_frames(lam, exact_weights(V, lo, hi, hf), mh.energy_offset, w, max_k)
            for (lam, V, *_), w in zip(systems, sampled)]


def _frames(lam, d0, off, sampled, max_k):
    """Result dict for one eigensystem from normalized sampled weights (len(SAMPLED), trials, dim)."""
    reach = lam[d0 > 1e-8]
    mu = 0.5 * (reach.max() + reach.min())          # device_odmd.centered_frame's mu, sector-local
    e0 = reach.min() + off
    budgets, order, trials = SAMPLED + SCALED, 2 * max_k - 1, sampled.shape[1]
    weights = np.concatenate([sampled, [scale_to(sampled[-1], d0, SAMPLED[-1], S) for S in SCALED]])
    mom = {f: np.array([[moments(lam, d, shift, order) for d in w] for w in weights])
           for f, shift in (("raw", 0.0), ("cen", mu))}
    shuffle = (np.arange(trials)[:, None] + np.arange(order + 1)[None, :]) % trials
    for f in ("raw", "cen"):
        mom[f + "_ind"] = np.take_along_axis(mom[f], shuffle[None], axis=1)
    out = {"budgets": budgets, "e0": e0, "mu": mu}
    for f, shift in (("exact", mu), ("exact_raw", 0.0)):          # noiseless PDS(K) - E0
        out[f] = np.array([pds_or_nan(moments(lam, d0, shift, order), K, off + shift) - e0
                           for K in range(1, max_k + 1)])
    for f, m in mom.items():
        shift = mu if f.startswith("cen") else 0.0
        out[f] = np.array([[[pds_or_nan(m[b, t], K, off + shift) - e0 for K in range(1, max_k + 1)]
                            for t in range(trials)] for b in range(len(budgets))])
    out["cond_raw"] = np.array([[[hankel_cond(mom["raw"][b, t], K) for K in range(1, max_k + 1)]
                                 for t in range(trials)] for b in range(len(budgets))])
    return out


def n2_study(trials: int = 16, seed: int = 2000):
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    return noise_study(build_molecular_hamiltonian(**N2), trials, seed=seed, paths=PATHS)


def pinned_n2_study(trials: int = 16, seed: int = 2000):
    """``n2_study`` in a fresh interpreter with ``PIN`` set (hash seed and BLAS threads are read at
    start-up), so every number is reproducible on a given platform. The caller's env is restored."""
    saved = os.environ.copy()
    os.environ.update(PIN)                          # inherited by the spawned child
    try:
        with ProcessPoolExecutor(1, mp_context=get_context("spawn")) as ex:
            return ex.submit(n2_study, trials, seed).result()
    finally:
        os.environ.clear()
        os.environ.update(saved)


def hit_rate(err: np.ndarray, tol: float = 1.6e-3) -> np.ndarray:
    """Fraction of trials with |PDS - E0| < tol (NaN = miss). err: (..., trials, K) -> (..., K)."""
    return np.mean(np.nan_to_num(np.abs(err), nan=np.inf) < tol, axis=-2)


def survives(res, b: int, tol: float = 1.6e-3, margin: float = 0.25) -> bool:
    """Pre-registered rule (spec G2): centering's float64 gain survives at budget index b iff some
    K in 5..8 has hit_cen(K) >= hit_cen(4) + margin AND hit_cen(K) >= hit_raw(K) + margin."""
    hc, hr = hit_rate(res["cen"][b], tol), hit_rate(res["raw"][b], tol)
    return any(hc[K - 1] >= hc[3] + margin and hc[K - 1] >= hr[K - 1] + margin for K in range(5, 9))


if __name__ == "__main__":
    import sys
    import time

    t0 = time.time()
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    base, *alt = pinned_n2_study(trials)
    print(f"N2 CAS(6,6)  E0={base['e0']:.8f}  mu={base['mu']:.4f}  trials={trials} (seeds 2000..)  "
          f"pinned {PIN}  wall {time.time() - t0:.0f} s")
    print("noiseless centered PDS(K)-E0 mHa, K=1..8:", " ".join(f"{e * 1e3:.4g}" for e in base["exact"]))
    for p, r in zip(PATHS, (base, *alt)):
        print(f"  noiseless raw K=6..8 mHa, path {p}:", " ".join(f"{e * 1e3:.4g}" for e in r["exact_raw"][5:]))
    print(f"hits = trials (of {trials}) with |PDS(K) - E0| < 1.6 mHa; 'all paths' = min-max over {PATHS}")
    for b, S in enumerate(base["budgets"]):
        tag = "sampled" if S in SAMPLED else "scaled"
        print(f"S={S:.0e} ({tag})  survives on pinned path: {survives(base, b)}  "
              f"by path: {''.join('TF'[not survives(r, b)] for r in (base, *alt))}")
        hit = {f: np.array([hit_rate(r[f][b]) for r in (base, *alt)]) * trials for f in ("raw", "cen")}
        hci, hri = (hit_rate(base[f][b]) * trials for f in ("cen_ind", "raw_ind"))
        for K in range(1, 9):
            med = {f: np.median(np.nan_to_num(np.abs(base[f][b, :, K - 1]), nan=np.inf)) * 1e3
                   for f in ("raw", "cen", "raw_ind", "cen_ind")}
            rng_ = {f: f"{hit[f][:, K-1].min():2.0f}-{hit[f][:, K-1].max():2.0f}" for f in hit}
            print(f"  K={K} median|err| mHa raw {med['raw']:9.4g} cen {med['cen']:9.4g} | "
                  f"indep raw {med['raw_ind']:9.4g} cen {med['cen_ind']:9.4g} | hits raw "
                  f"{hit['raw'][0, K-1]:2.0f} cen {hit['cen'][0, K-1]:2.0f} indep raw {hri[K-1]:2.0f} "
                  f"cen {hci[K-1]:2.0f} | all paths raw {rng_['raw']} cen {rng_['cen']}")
    n = len(SAMPLED)
    ok = base["cond_raw"][:n] < 1e10
    raw, cen = base["raw"][:n], base["cen"][:n]
    dev = np.nan_to_num(np.abs(cen - base["exact"]), nan=np.inf)       # shot-noise deviation
    shot = min(np.median(dev[b][ok[b][:, k], k]) for b in range(n) for k in range(3))
    print(f"G1: where raw cond < 1e10, max |raw - cen| {np.nanmax(np.abs(raw - cen)[ok]):.2e} Ha; "
          f"no real root in either frame at (budget, trial, K-1) "
          f"{np.argwhere(np.isnan(raw) & np.isnan(cen) & ok).tolist()}; "
          f"{int((np.isnan(raw) ^ np.isnan(cen))[ok].sum())} pairs disagree on a root; "
          f"K<=3 pairs with raw cond >= 1e10: {np.argwhere(~ok[..., :3]).tolist()}; "
          f"min median shot-noise deviation at K<=3 {shot:.2e} Ha; median |<H> - <H>_exact| mHa "
          f"per sampled budget {[float(f'{np.median(dev[b][:, 0]) * 1e3:.4g}') for b in range(n)]}")
    i6 = SAMPLED.index(10**6)
    g3 = [np.median(np.nan_to_num(np.abs(base[f][i6, :, 7]), nan=np.inf)) for f in ("cen_ind", "cen")]
    print(f"G3: 1e6 K=8 centered median |err| independent {g3[0] * 1e3:.4g} mHa / shared "
          f"{g3[1] * 1e3:.4g} mHa = {g3[0] / g3[1]:.0f}x")
