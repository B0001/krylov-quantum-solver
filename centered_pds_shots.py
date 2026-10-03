#!/usr/bin/env python3
"""
Does centering change how shot noise in <H^n> amplifies into PDS(K)?  (specs/SPEC_centered_pds_shots.md)

Circuit-level noise model: ONE classical-shadow dataset of the Hartree-Fock state (random-Pauli
snapshots, Huang-Kueng-Preskill; same estimator as ``classical_shadows.py``) feeds EVERY moment, raw
and centered. Each snapshot is the product operator rho_s = (x)_q (I + 3 s_q P_q)/2, so every
moment <O> is estimated by Tr(O rho_avg) from the same shots -- the noise in mu_1..mu_15 is
correlated exactly as it would be on hardware, not independent Gaussian per moment.

Symmetry projection (standard post-selection, classical post-processing of the same shots): H
commutes with the (N_alpha, N_beta) projector P, so <P H^n>/<P> is estimated from the same snapshots
using only the sector block of rho_avg. That block is a 400x400 matrix for N2 CAS(6,6), built by a
hi/lo-qubit GEMM, which is what makes 10^6 snapshots cheap. Moments then come from the sector
eigenbasis: mu_n = sum_i (lam_i - mu)^n d_i with d = diag(V^dag rho_sec V) / tr.

Larger budgets (``scale_to``) rescale a sampled deviation d - d_exact by sqrt(S0/S). Covariance is
exact under that rescaling (cov of a mean is 1/S); higher cumulants are not.
"""
from __future__ import annotations

import numpy as np

from moment_expansion import pds_energy

_PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=complex)  # X,Y,Z


def _popcount_states(n_bits: int, k: int) -> np.ndarray:
    return np.array([x for x in range(2 ** n_bits) if bin(x).count("1") == k])


def sector_eigensystem(mh):
    """(lam, V, hf_bits, lo_states, hi_states) of H restricted to HF's (N_alpha, N_beta) sector.

    JW ordering (qiskit-nature): alpha spin-orbitals are qubits 0..n-1 (lo), beta n..2n-1 (hi).
    Sector basis order is (hi, lo) row-major, matching ``shadow_sector_weights``.
    """
    n = mh.num_spatial_orbitals
    lo = _popcount_states(n, mh.num_particles[0])
    hi = _popcount_states(n, mh.num_particles[1])
    idx = (hi[:, None] << n | lo[None, :]).ravel()
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsr()[idx][:, idx].toarray()
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


def shadow_sector_weights(bases, signs, V, lo, hi, chunk: int = 50_000) -> np.ndarray:
    """Eigenbasis weights d_i = <v_i|rho_sec|v_i> from one shadow dataset, UNnormalized: sum(d) is the
    shadow estimate of <P_sector>. PDS is invariant under a common scale of all moments."""
    n = bases.shape[1] // 2
    m_hi, m_lo = len(hi), len(lo)
    acc = np.zeros((m_hi, m_hi, m_lo, m_lo), dtype=complex)
    for s0 in range(0, bases.shape[0], chunk):
        b, s = bases[s0:s0 + chunk], signs[s0:s0 + chunk]
        Bh = _half_block(b, s, range(n, 2 * n), hi).reshape(len(b), -1)
        Bl = _half_block(b, s, range(n), lo).reshape(len(b), -1)
        acc += (Bh.T @ Bl).reshape(m_hi, m_hi, m_lo, m_lo)
    rho = acc.transpose(0, 2, 1, 3).reshape(m_hi * m_lo, m_hi * m_lo)   # [(h,l),(h',l')]
    return np.real(np.einsum("ji,jk,ki->i", V.conj(), rho, V)) / bases.shape[0]


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


def noise_study(mh, trials: int = 16, max_k: int = 8, seed: int = 2000):
    """PDS(K) - E0 for raw and centered frames, per budget x trial x K, from shared-shot shadows.

    Also an INDEPENDENT-per-moment contrast with identical marginals: moment n of trial t is taken
    from trial (t + n) % trials, in each frame separately (the "measure each power on its own shots"
    model a prior attempt used). Returns a dict of arrays shaped (len(budgets), trials, max_k).
    """
    lam, V, hf, lo, hi = sector_eigensystem(mh)
    d0 = exact_weights(V, lo, hi, hf)
    reach = lam[d0 > 1e-8]
    mu = 0.5 * (reach.max() + reach.min())          # device_odmd.centered_frame's mu, sector-local
    off, e0 = mh.energy_offset, reach.min() + mh.energy_offset
    budgets = SAMPLED + SCALED
    order = 2 * max_k - 1
    mom = {f: np.empty((len(budgets), trials, order + 1)) for f in ("raw", "cen")}
    for t in range(trials):
        bases, signs = sample_hf_shadow(hf, SAMPLED[-1], np.random.default_rng(seed + t))
        for b, S in enumerate(budgets):
            if S <= SAMPLED[-1]:
                d = shadow_sector_weights(bases[:S], signs[:S], V, lo, hi)
                d = d_top = d / d.sum()
            else:
                d = scale_to(d_top, d0, SAMPLED[-1], S)
            mom["raw"][b, t] = moments(lam, d, 0.0, order)
            mom["cen"][b, t] = moments(lam, d, mu, order)
    shuffle = (np.arange(trials)[:, None] + np.arange(order + 1)[None, :]) % trials
    for f in ("raw", "cen"):
        mom[f + "_ind"] = np.take_along_axis(mom[f], shuffle[None], axis=1)
    out = {"budgets": budgets, "e0": e0, "mu": mu}
    for f, m in mom.items():
        shift = mu if f.startswith("cen") else 0.0
        out[f] = np.array([[[pds_or_nan(m[b, t], K, off + shift) - e0 for K in range(1, max_k + 1)]
                            for t in range(trials)] for b in range(len(budgets))])
    out["cond_raw"] = np.array([[[hankel_cond(mom["raw"][b, t], K) for K in range(1, max_k + 1)]
                                 for t in range(trials)] for b in range(len(budgets))])
    return out


def hit_rate(err: np.ndarray, tol: float = 1.6e-3) -> np.ndarray:
    """Fraction of trials with |PDS - E0| < tol (NaN = miss). err: (..., trials, K) -> (..., K)."""
    return np.mean(np.nan_to_num(np.abs(err), nan=np.inf) < tol, axis=-2)


def survives(res, b: int, tol: float = 1.6e-3, margin: float = 0.25) -> bool:
    """Pre-registered rule (spec G2): centering's float64 gain survives at budget index b iff some
    K in 5..8 has hit_cen(K) >= hit_cen(4) + margin AND hit_cen(K) >= hit_raw(K) + margin."""
    hc, hr = hit_rate(res["cen"][b], tol), hit_rate(res["raw"][b], tol)
    return any(hc[K - 1] >= hc[3] + margin and hc[K - 1] >= hr[K - 1] + margin for K in range(5, 9))


if __name__ == "__main__":
    import os
    import sys

    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    mh = build_molecular_hamiltonian(atom="N 0 0 0; N 0 0 1.10", basis="sto3g",
                                     active_electrons=6, active_orbitals=6)
    res = noise_study(mh, trials=int(sys.argv[1]) if len(sys.argv) > 1 else 16)
    print(f"N2 CAS(6,6)  E0={res['e0']:.8f}  mu={res['mu']:.4f}  threads={os.environ.get('OMP_NUM_THREADS')}")
    for b, S in enumerate(res["budgets"]):
        tag = "sampled" if S in SAMPLED else "scaled"
        print(f"S={S:.0e} ({tag})  survives={survives(res, b)}")
        hc, hr = hit_rate(res["cen"][b]), hit_rate(res["raw"][b])
        hci, hri = hit_rate(res["cen_ind"][b]), hit_rate(res["raw_ind"][b])
        for K in range(1, 9):
            med = {f: np.median(np.nan_to_num(np.abs(res[f][b, :, K - 1]), nan=np.inf)) * 1e3
                   for f in ("raw", "cen", "raw_ind", "cen_ind")}
            print(f"  K={K} median|err| mHa raw {med['raw']:9.4g} cen {med['cen']:9.4g} | "
                  f"indep raw {med['raw_ind']:9.4g} cen {med['cen_ind']:9.4g} | "
                  f"hit raw {hr[K-1]:.2f} cen {hc[K-1]:.2f} indep raw {hri[K-1]:.2f} cen {hci[K-1]:.2f}")
