#!/usr/bin/env python3
"""
quantum_krylov_solver.py -- a real real-time quantum Krylov subspace solver.

This REPLACES the broken subspace path (the near-identity qDRIFT "basis" in
``quantum_sampler.py`` and the noise-injecting ``StabilizedSubspaceShifter`` /
``execute_subspace_sweep`` in ``orchestrate_hybrid_pipeline.py``). The original code
(a) started from the empty vacuum |00..0>, (b) built each "Krylov vector" by appending
ONE single-Pauli rotation with an infinitesimal angle, so all states were ~identical
(overlap matrix rank ~1), and (c) solved a generalized eigenproblem with asymmetric
Gaussian noise added to H -- which is unbounded below and routinely returned energies
hundreds of Hartree *below* the true ground state.

Real-time quantum Krylov subspace diagonalisation (QKSD) instead builds the genuine
Krylov space from the Hartree-Fock reference and real-time propagation:

    |phi_k> = U^k |phi_HF>,   U = exp(-i H dt),   k = 0, 1, ..., M-1

    H_ij = <phi_i| H |phi_j>,   S_ij = <phi_i|phi_j>        (both Hermitian)

    solve   H c = E S c   (generalised eigenproblem)

The generalised problem is solved by a *thresholded canonical orthogonalisation*:
diagonalise S, drop directions whose eigenvalue is below ``threshold * lambda_max``
(these are the linearly dependent / numerically singular Krylov directions), and
diagonalise H in the well-conditioned remainder. Because every retained vector is a
genuine state in Hilbert space, the returned energy is a Rayleigh quotient of H and is
therefore **variationally bounded below by the true ground-state energy** -- it can
never dip below FCI (contrast the original code).

References: Parrish & McMahon 2019; Stair, Huang & Evangelista, JCTC 2020 (multireference
quantum Krylov); Klymko et al., PRX Quantum 3, 020323 (2022, real-time evolution Krylov);
Cortes & Gray, PRA 105, 022417 (2022).

Scope: time evolution here is performed *exactly* on the statevector (``expm_multiply``),
which isolates and validates the subspace algorithm itself. Trotter/qDRIFT circuit
synthesis and a hardware noise model are the next step (Phase 3 in REFACTOR_PLAN.md);
this solver is the correct, validated target they must reproduce.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, expm_multiply

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian

import math


def expm_multiply_taylor(H, v, t, lambda_bound, theta: float = 0.9, order: int = 18):
    """Backend-agnostic action of ``exp(-i t H) @ v`` via a scaled Taylor series.

    Uses ONLY sparse mat-vec (``H.dot``) and vector axpy, so the identical code runs on
    NumPy + scipy.sparse OR CuPy + cupyx.scipy.sparse arrays -- this is what lets the GPU
    backend evolve statevectors without a CuPy ``expm_multiply`` (which cupyx does not provide).

    ``lambda_bound`` must be an upper bound on ``||H||_2`` (the sum of |Pauli coefficients|
    works: each Pauli has spectral norm 1). It sets the number of sub-steps ``s`` so each
    Taylor block has spectral radius <= ``theta`` and the truncation at ``order`` terms is at
    machine precision.  Validated on CPU against ``scipy.sparse.linalg.expm_multiply`` in
    ``tests/test_gpu_backend.py``.
    """
    s = max(1, int(math.ceil(abs(t) * float(lambda_bound) / theta)))
    c = -1j * t / s
    w = v
    for _ in range(s):
        term = w
        out = w
        for k in range(1, order + 1):
            term = H.dot(term) * (c / k)
            out = out + term
        w = out
    return w


def ritz_spectrum(H, S, threshold: float = 1e-10, noise_floor: float = 0.0):
    """Thresholded canonical orthogonalisation for ``H c = E S c``; returns (ritz, rank).

    ``ritz`` is the full ascending array of Ritz values (the eigenvalues of H projected onto the
    well-conditioned subspace) -- ``ritz[0]`` is the ground-state estimate, ``ritz[1:]`` the
    excited-state estimates (see ``solve_excited``).

    The overlap cutoff is the larger of a relative floor (``threshold * lambda_max``) and a
    noise-aware absolute floor (``noise_floor``): overlap directions buried below the
    sampling-noise level carry no signal and must be dropped, otherwise dividing by their tiny
    eigenvalues amplifies noise and lets the estimate drift below the true minimum
    (cf. Epperly, Lin & Nakatsukasa, SIAM J. Matrix Anal. Appl. 43, 1263, 2022).
    """
    s_vals, s_vecs = eigh(S)
    cutoff = max(threshold * s_vals.max(), noise_floor)
    keep = s_vals > cutoff
    if not np.any(keep):                              # degenerate guard
        keep = s_vals >= s_vals.max()
    s_keep, V = s_vals[keep], s_vecs[:, keep]
    X = V / np.sqrt(s_keep)                           # canonical S^{-1/2} on kept subspace
    H_proj = X.conj().T @ H @ X
    H_proj = 0.5 * (H_proj + H_proj.conj().T)
    ritz = np.linalg.eigvalsh(H_proj).real            # ascending
    return ritz, int(keep.sum())


def solve_generalized_eig(H, S, threshold: float = 1e-10, noise_floor: float = 0.0):
    """Ground-state Ritz value of ``H c = E S c``; returns (energy, rank).

    Thin wrapper over :func:`ritz_spectrum` (returns its lowest Ritz value). Shared by the
    exact-evolution solver and the Trotter-circuit solver; the variational floor argument applies
    to every Ritz value, not just the lowest.
    """
    ritz, rank = ritz_spectrum(H, S, threshold, noise_floor)
    return float(ritz[0]), rank


@dataclass
class KrylovStep:
    """Result at one Krylov dimension."""
    dim: int          # requested Krylov dimension M
    rank: int         # effective subspace rank kept after thresholding
    energy: float     # total energy in Ha (electronic eigenvalue + offset)


@dataclass
class ExcitedKrylovStep:
    """Low-lying spectrum at one Krylov dimension (ground + excited Ritz values)."""
    dim: int                # requested Krylov dimension M
    rank: int               # effective subspace rank kept after thresholding
    energies: List[float]   # ascending total energies in Ha (offset included); [0] is the ground


class QuantumKrylovSolver:
    """Real-time quantum Krylov subspace diagonalisation on top of a MolecularHamiltonian."""

    def __init__(
        self,
        molecular_hamiltonian: MolecularHamiltonian,
        dt: Optional[float] = None,
        threshold: float = 1e-10,
        noise_sigma: float = 0.0,
        seed: Optional[int] = None,
        device: str = "cpu",
    ):
        """
        Args:
            molecular_hamiltonian: built by ``build_molecular_hamiltonian`` -- supplies the
                qubit Hamiltonian, the Hartree-Fock reference, and the energy offset.
            dt: real-time evolution step. If ``None``, a robust default ``pi / spectral_width``
                is used, with the width estimated by Lanczos (no full diagonalisation).
            threshold: relative cutoff for canonical orthogonalisation; Krylov directions with
                overlap-eigenvalue below ``threshold * lambda_max(S)`` are dropped.
            noise_sigma: standard deviation of statistical shot noise on the measured H and S
                matrix elements (0 = exact statevector). Models the finite-sampling error of a
                Hadamard-test estimate; perturbations are added *Hermitian-symmetrically*
                (unlike the original code's symmetry-breaking additive noise on H). A useful
                scale is ``1/sqrt(shots)`` (see ``noise.shot_noise_sigma``).
            seed: RNG seed for the shot noise (reproducible studies).
            device: ``"cpu"`` (default, scipy.sparse + scipy.expm_multiply -- the validated path
                exercised by the test suite) or ``"gpu"`` (CuPy + cupyx.scipy.sparse, with the
                CPU-validated ``expm_multiply_taylor`` for the time step). The GPU path lets the
                statevector simulation reach larger qubit counts on an NVIDIA GPU (an A100 80GB
                holds ~32 qubits); it requires ``cupy`` and is validated only on its CPU-fallback
                math here -- run it on a GPU node to confirm end to end.
        """
        self.mh = molecular_hamiltonian
        self.offset = molecular_hamiltonian.energy_offset
        self.threshold = threshold
        self.noise_sigma = float(noise_sigma)
        self._rng = np.random.default_rng(seed)
        self.device = device

        H_sparse = molecular_hamiltonian.qubit_hamiltonian.to_matrix(sparse=True)
        if device == "gpu":
            try:
                import cupy as cp
                import cupyx.scipy.sparse as cxsp
            except Exception as exc:  # pragma: no cover - exercised only with a CuPy/GPU node
                raise ImportError(
                    "device='gpu' requires cupy (pip install cupy-cuda12x) and an NVIDIA GPU."
                ) from exc
            self._cp = cp
            self._H = cxsp.csr_matrix(H_sparse.astype(complex))
            self._psi0 = cp.asarray(molecular_hamiltonian.hf_state().data, dtype=complex)
            # Upper bound on ||H||_2 for the Taylor sub-step count (sum of |Pauli coeffs|).
            self._lambda = float(np.sum(np.abs(molecular_hamiltonian.qubit_hamiltonian.coeffs)))
        else:
            self._cp = None
            self._H = H_sparse.tocsc()
            self._psi0 = np.asarray(molecular_hamiltonian.hf_state().data, dtype=complex)

        self.dt = float(dt) if dt is not None else self._default_dt()
        self._basis: List[np.ndarray] = [self._psi0.copy()]  # cached Krylov vectors

    def _hermitian_shot_noise(self, A: np.ndarray, unit_diagonal: bool) -> np.ndarray:
        """Add Hermiticity-preserving Gaussian shot noise of scale ``self.noise_sigma``."""
        if self.noise_sigma <= 0.0:
            return A
        g = (self._rng.normal(scale=self.noise_sigma, size=A.shape)
             + 1j * self._rng.normal(scale=self.noise_sigma, size=A.shape))
        perturbed = A + 0.5 * (g + g.conj().T)          # Hermitian perturbation
        if unit_diagonal:
            np.fill_diagonal(perturbed, 1.0)            # <phi_i|phi_i> = 1 exactly
        return 0.5 * (perturbed + perturbed.conj().T)

    # -- time step -----------------------------------------------------------
    def _spectral_radius(self) -> float:
        """Largest |eigenvalue| via power iteration (device-agnostic; only sparse mat-vec)."""
        xp = self._cp if self.device == "gpu" else np
        v = xp.asarray(self._psi0, dtype=complex)
        v = v / xp.linalg.norm(v)
        r = 0.0
        for _ in range(50):
            w = self._H.dot(v)
            nrm = float(xp.linalg.norm(w))
            if nrm == 0.0:
                break
            v, r = w / nrm, nrm
        return r

    def _default_dt(self) -> float:
        """dt = pi / spectral_width."""
        if self.device == "gpu":
            # cupyx lacks a robust complex-Hermitian eigsh; estimate the width from the spectral
            # radius (power iteration). Pass dt explicitly for large GPU runs to avoid this.
            width = 2.0 * self._spectral_radius()
            return np.pi / width if width > 0 else np.pi
        try:
            e_max = float(eigsh(self._H, k=1, which="LA", return_eigenvectors=False)[0])
            e_min = float(eigsh(self._H, k=1, which="SA", return_eigenvectors=False)[0])
            width = e_max - e_min
        except Exception:
            width = float(np.linalg.eigvalsh(self._H.toarray())[[0, -1]] @ [-1, 1])
        if not np.isfinite(width) or width <= 0:
            width = 2.0 * float(np.sum(np.abs(self.mh.qubit_hamiltonian.coeffs)))
        return np.pi / width

    # -- Krylov basis --------------------------------------------------------
    def _ensure_basis(self, dim: int) -> None:
        """Extend the cached Krylov basis to at least ``dim`` vectors."""
        while len(self._basis) < dim:
            if self.device == "gpu":
                nxt = expm_multiply_taylor(self._H, self._basis[-1], self.dt, self._lambda)
            else:
                nxt = expm_multiply(-1j * self.dt * self._H, self._basis[-1])
            self._basis.append(nxt)

    def _subspace_matrices(self, dim: int):
        self._ensure_basis(dim)
        xp = self._cp if self.device == "gpu" else np
        B = xp.array(self._basis[:dim])              # (M, N)
        S = B.conj() @ B.T                           # <phi_i|phi_j>
        H = B.conj() @ self._H.dot(B.T)              # <phi_i|H|phi_j>
        if self.device == "gpu":                     # M x M is tiny -> finish on the CPU
            H, S = self._cp.asnumpy(H), self._cp.asnumpy(S)
        # Hermitise away numerical asymmetry (a *symmetric* correction, unlike the old
        # code which ADDED asymmetric noise to H). Optional shot noise is also added
        # Hermitian-symmetrically, so the eigenproblem stays well-posed.
        H = self._hermitian_shot_noise(0.5 * (H + H.conj().T), unit_diagonal=False)
        S = self._hermitian_shot_noise(0.5 * (S + S.conj().T), unit_diagonal=True)
        return H, S

    # -- generalised eigenproblem -------------------------------------------
    def _solve_generalized(self, H, S):
        """Thresholded canonical orthogonalisation with a noise-aware overlap cutoff."""
        return solve_generalized_eig(H, S, self.threshold, 5.0 * self.noise_sigma)

    # -- public API ----------------------------------------------------------
    def solve(self, krylov_dim: int) -> KrylovStep:
        """Estimate the ground-state energy from an ``krylov_dim``-vector Krylov space."""
        if krylov_dim < 1:
            raise ValueError("krylov_dim must be >= 1")
        H, S = self._subspace_matrices(krylov_dim)
        energy, rank = self._solve_generalized(H, S)
        return KrylovStep(dim=krylov_dim, rank=rank, energy=energy + self.offset)

    def convergence(self, max_dim: int) -> List[KrylovStep]:
        """Energies for M = 1 .. max_dim (basis built once and reused)."""
        self._ensure_basis(max_dim)
        return [self.solve(m) for m in range(1, max_dim + 1)]

    def solve_excited(self, krylov_dim: int, n_states: Optional[int] = None) -> ExcitedKrylovStep:
        """Estimate the low-lying spectrum (ground + excited) from a Krylov space.

        Returns the lowest ``n_states`` Ritz values (all of them if ``n_states`` is None) of the
        same subspace ``solve`` uses -- ``energies[0]`` is identical to ``solve(krylov_dim).energy``.
        Each Ritz value is variationally above the corresponding exact eigenvalue (Cauchy
        interlacing); excited states are reachable only insofar as |HF> overlaps them. See
        specs/SPEC_qksd_excited.md.
        """
        if krylov_dim < 1:
            raise ValueError("krylov_dim must be >= 1")
        H, S = self._subspace_matrices(krylov_dim)
        ritz, rank = ritz_spectrum(H, S, self.threshold, 5.0 * self.noise_sigma)
        if n_states is not None:
            ritz = ritz[:n_states]
        energies = [float(e) + self.offset for e in ritz]
        return ExcitedKrylovStep(dim=krylov_dim, rank=rank, energies=energies)


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    for name, spec in {
        "H2": dict(atom="H 0 0 0; H 0 0 0.74"),
        "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    }.items():
        mh = build_molecular_hamiltonian(**spec)
        fci = mh.ground_state_energy()
        solver = QuantumKrylovSolver(mh)
        print(f"\n{name}: FCI={fci:.6f} Ha, dt={solver.dt:.4f}")
        for step in solver.convergence(10):
            print(f"  M={step.dim:2d} rank={step.rank:2d}  "
                  f"E={step.energy:.6f}  err={(step.energy - fci)*1e3:+.3f} mHa")
