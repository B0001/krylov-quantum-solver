#!/usr/bin/env python3
"""
Observable dynamic mode decomposition (ODMD) -- ground-state energy from the survival amplitude
ALONE (arXiv:2306.01858, Shen/Klymko/Sud/Williams-Young/de Jong/Van Beeumen).

Every other subspace rung in this repo (QKSD, SKQD, MSD) measures a projected Hamiltonian matrix,
whose per-element sampling variance carries the Pauli 1-norm lambda (or MSD's stencil 1-norm).
ODMD needs neither: the complex overlap time series

    s_k = <phi_0| e^{-i k tau H} |phi_0>,   k = 0..K-1

-- the *first row* of the S matrix the QKSD pipeline already measures, K numbers total -- is a sum
of unimodular modes  s_k = sum_n p_n e^{-i E_n k tau}  with p_n = |<phi_0|E_n>|^2 >= 0, so a linear
one-step propagator fitted to it (dynamic mode decomposition on Hankel data matrices) has
eigenvalues e^{-i E_n tau}. The ground-state energy is the minimum eigenphase. The SVD-truncated
pseudoinverse in the DMD least-squares is the noise-robustness mechanism: directions of the Hankel
matrix below the shot-noise floor are discarded instead of amplified.

HONEST SCOPE (see specs/SPEC_odmd.md): reproduction, exact-statevector overlaps, idealized i.i.d.
per-element shot noise (msd.py Hadamard-test conventions). We use the complex signal (both
Hadamard quadratures); the paper's Re-only variant is out of scope. ODMD is NOT variational --
eigenphases can dip below FCI at small K (gated in tests/test_odmd_spec.py G2). The matched-budget
advantage over KQD is measured against this repo's LCU/Hadamard noise model for H elements.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eig, hankel, svd
from scipy.sparse import identity
from scipy.sparse.linalg import expm_multiply

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian


@dataclass
class ODMDProblem:
    """Noiseless survival amplitudes in the centered frame (H - mu), msd.py conventions.

    Energies from :func:`odmd_energy` are in this frame; add ``offset`` to recover the total
    energy. ``ref`` is the exact (FCI) ground energy in the frame, so ``ref + offset`` is the
    total FCI energy.
    """
    n: int
    tau: float             # pi / W, W = HF-reachable spectral width
    mu: float              # spectral center (the energy-level shift)
    s: np.ndarray          # complex s_k, k = 0..n-1 (s[0] = 1)
    dim: int               # Hilbert-space dimension d (sets the Hadamard-test variance constant)
    offset: float          # mh.energy_offset + mu  -> total energy
    ref: float             # exact ground energy in the centered frame


def build_odmd_problem(mh: MolecularHamiltonian, n: int = 20) -> ODMDProblem:
    """Exact survival amplitudes s_0..s_{n-1} with the msd.py energy-level shift and tau = pi/W."""
    H_full = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    psi0 = np.asarray(mh.hf_state().data, dtype=complex)
    w_eig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    reach = w_eig[(np.abs(V.conj().T @ psi0) ** 2) > 1e-8].real
    width = float(reach.max() - reach.min())
    mu = float(0.5 * (reach.max() + reach.min()))
    H_s = (H_full - mu * identity(H_full.shape[0], format="csc")).tocsc()
    tau = float(np.pi / width)
    s = np.array([psi0.conj() @ expm_multiply(-1j * (k * tau) * H_s, psi0) for k in range(n)])
    return ODMDProblem(n=n, tau=tau, mu=mu, s=s, dim=int(H_full.shape[0]),
                       offset=mh.energy_offset + mu, ref=float(w_eig[0] - mu))


def _dmd_modes(s, tau, mod_window, cutoff_abs=0.0, cutoff_rel=0.0):
    """Near-unimodular DMD eigenvalues of the Hankel data matrices; returns (lambdas, rank).

    Hankel DMD: X[i,j] = s_{i+j}, X'[i,j] = s_{i+j+1} with d = K//2 rows; the reduced operator
    A = U_r^H X' V_r Sigma_r^{-1} keeps only singular directions with sigma above the larger of
    an absolute cutoff (``cutoff_abs`` -- the noise edge, see :func:`noise_edge`) and a relative
    floor (``cutoff_rel * sigma_max``). Eigenvalues ~ e^{-i E_n tau}; modes with ||lambda|-1| >=
    ``mod_window`` are noise artifacts and dropped (the centered frame keeps |E_n| tau <= pi/2 --
    no phase wrapping).
    """
    s = np.asarray(s, dtype=complex)
    K = len(s)
    if K < 4:
        raise ValueError("need at least 4 signal points")
    d = K // 2
    X = hankel(s[:d], s[d - 1:K - 1])                 # (d, K-d)
    Xp = hankel(s[1:d + 1], s[d:K])
    U, sig, Vh = svd(X, full_matrices=False)
    cutoff = max(cutoff_abs, cutoff_rel * sig[0])
    rank = max(int(np.sum(sig > cutoff)), 1)
    A = U[:, :rank].conj().T @ Xp @ Vh[:rank, :].conj().T / sig[:rank]
    lam = eig(A, right=False)
    keep = np.abs(np.abs(lam) - 1.0) < mod_window
    if not np.any(keep):                              # degenerate guard
        keep = np.ones_like(lam, dtype=bool)
    return lam[keep], rank


def odmd_energy(s, tau: float, svd_threshold: float = 1e-10, mod_window: float = 0.2):
    """Ground-state energy (centered frame) from the signal ``s``; returns (energy, rank).

    The minimum retained eigenphase of :func:`_dmd_modes` with the relative noise floor
    ``svd_threshold * sigma_max`` (the SPEC_odmd.md ground-state semantics, gates pinned).
    """
    lam, rank = _dmd_modes(s, tau, mod_window, cutoff_rel=svd_threshold)
    return float(np.min(-np.angle(lam) / tau)), rank


def noise_edge(sigma: float, d: int, m: int, c: float = 1.2) -> float:
    """Absolute Hankel singular-value cutoff c * sigma * (sqrt(d) + sqrt(m)).

    The largest singular value of a d x m matrix with i.i.d. noise of scale ``sigma`` sits at the
    Marchenko-Pastur edge ~ sigma*(sqrt d + sqrt m); singular directions above it carry signal.
    Hankel noise entries are antidiagonally *correlated*, so ``c = 1.2`` is a calibrated
    heuristic, not a theorem (see specs/SPEC_odmd_excited.md R1). Unlike SPEC_odmd.md's relative
    ``5 sigma * sigma_max`` floor -- which the dominant ground mode (p0 ~ 0.95) inflates until it
    swallows the excited singular value ~ p1*sqrt(dm) -- this cutoff keeps ANY mode whose
    amplitude clears the noise: mode n is visible iff p_n*sqrt(dm) > c*sigma*(sqrt d + sqrt m),
    so signal depth buys excited-state visibility as ~ sqrt(K).
    """
    return float(c * sigma * (np.sqrt(d) + np.sqrt(m)))


def odmd_spectrum(s, tau: float, cutoff: float = 0.0, mod_window: float = 0.2,
                  amp_floor: float = 0.0):
    """Low-lying spectrum (centered frame) from the SAME signal; returns (energies, amps, rank).

    All retained DMD eigenphases sorted ascending -- ``energies[0]`` is the ground state,
    ``energies[1]`` the first HF-visible excited state (no extra measurements beyond the
    :func:`odmd_energy` signal). ``cutoff`` is the *absolute* singular-value truncation (use
    :func:`noise_edge` when a noise scale is known; a 1e-10 relative floor is always applied).
    ``amps`` are |a_n| from a Vandermonde least-squares refit s_k ~ sum_n a_n lambda_n^k
    (a_n ~ p_n = |<phi_0|E_n>|^2); modes with amplitude <= ``amp_floor`` are dropped.
    NOT variational -- see SPEC_odmd.md G2 / SPEC_odmd_excited.md.
    """
    s = np.asarray(s, dtype=complex)
    lam, rank = _dmd_modes(s, tau, mod_window, cutoff_abs=cutoff, cutoff_rel=1e-10)
    V = np.vander(lam, N=len(s), increasing=True).T   # (K, n_modes): V[k, i] = lam_i^k
    a, *_ = np.linalg.lstsq(V, s, rcond=None)
    E = -np.angle(lam) / tau
    order = np.argsort(E)
    E, a = E[order], np.abs(a[order])
    m = a > amp_floor
    return E[m], a[m], rank


def sample_odmd_spectrum(prob: ODMDProblem, shots: int, seed: int, n: int | None = None,
                         c: float = 1.2, amp_floor: float = 0.0) -> np.ndarray:
    """One shot-noisy spectrum draw (centered frame) with noise-edge truncation.

    Same per-element Hadamard-test noise as :func:`sample_odmd_energy` (s_0 = 1 exactly); the
    truncation is the absolute :func:`noise_edge` for the Hankel shape (d, n-d). Returns the
    ascending retained energies -- fewer than the exact count when a mode's amplitude is below
    the visibility law (e.g. length 1 when the excited state is invisible at this budget/depth).
    """
    rng = np.random.default_rng(seed)
    n = prob.n if n is None else n
    sigma = np.sqrt(2.0 * (2.0 - 1.0 / prob.dim) / shots)
    g = rng.normal(0, sigma / np.sqrt(2), n) + 1j * rng.normal(0, sigma / np.sqrt(2), n)
    g[0] = 0.0
    d = n // 2
    energies, _, _ = odmd_spectrum(prob.s[:n] + g, prob.tau,
                                   cutoff=noise_edge(sigma, d, n - d, c), amp_floor=amp_floor)
    return energies


def sample_odmd_energy(prob: ODMDProblem, shots: int, seed: int,
                       n: int | None = None, svd_threshold: float | None = None) -> float:
    """One shot-noisy ODMD estimate (centered frame) from the first ``n`` signal elements.

    Per-element noise follows msd.py's Hadamard-test analysis: each complex overlap gets i.i.d.
    Gaussian noise of total variance 2(2 - 1/d)/shots, s_0 = 1 exactly. The default SVD threshold
    is noise-aware, 5 * sigma_s (relative to sigma_max) -- the ODMD analogue of the Krylov
    solver's noise-aware overlap floor.
    """
    rng = np.random.default_rng(seed)
    n = prob.n if n is None else n
    sigma = np.sqrt(2.0 * (2.0 - 1.0 / prob.dim) / shots)
    if svd_threshold is None:
        svd_threshold = 5.0 * sigma
    g = rng.normal(0, sigma / np.sqrt(2), n) + 1j * rng.normal(0, sigma / np.sqrt(2), n)
    g[0] = 0.0
    return odmd_energy(prob.s[:n] + g, prob.tau, svd_threshold)[0]


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from msd import build_msd_problem, sample_ground_energy

    cases = {
        "H2": build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74"),
        "N2 CAS(6,6)": build_molecular_hamiltonian(
            atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
    }
    for name, mh in cases.items():
        prob = build_odmd_problem(mh, n=20)
        fci = mh.ground_state_energy()
        print("=" * 72)
        print(f"{name}: tau={prob.tau:.3f}  FCI={fci:.6f} Ha")
        for K in (8, 12, 16, 20):
            e, r = odmd_energy(prob.s[:K], prob.tau)
            print(f"  noiseless K={K:2d} rank={r:2d}  err={(e - prob.ref) * 1e3:+8.4f} mHa")
        kqd = build_msd_problem(mh, n=8, order=8)
        for shots in (1e4, 1e5):
            m_o = np.median([abs(sample_odmd_energy(prob, int(shots), sd, n=16) - prob.ref)
                             for sd in range(100)])
            m_k = np.median([abs(sample_ground_energy(kqd, int(shots), "kqd", sd) - prob.ref)
                             for sd in range(100)])
            print(f"  shots={int(shots):>7} (matched 16 elements): median ODMD={m_o * 1e3:8.3f} "
                  f"mHa   KQD={m_k * 1e3:8.3f} mHa   advantage={m_k / m_o:.1f}x")
    print("=" * 72)
