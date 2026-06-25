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


def solve_generalized_eig(H, S, threshold: float = 1e-10, noise_floor: float = 0.0):
    """Thresholded canonical orthogonalisation for ``H c = E S c``; returns (energy, rank).

    The overlap cutoff is the larger of a relative floor (``threshold * lambda_max``) and a
    noise-aware absolute floor (``noise_floor``): overlap directions buried below the
    sampling-noise level carry no signal and must be dropped, otherwise dividing by their tiny
    eigenvalues amplifies noise and lets the estimate drift below the true minimum
    (cf. Epperly, Lin & Nakatsukasa, SIAM J. Matrix Anal. Appl. 43, 1263, 2022).

    Shared by the exact-evolution solver and the Trotter-circuit solver.
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
    energy = float(np.linalg.eigvalsh(H_proj)[0].real)
    return energy, int(keep.sum())


@dataclass
class KrylovStep:
    """Result at one Krylov dimension."""
    dim: int          # requested Krylov dimension M
    rank: int         # effective subspace rank kept after thresholding
    energy: float     # total energy in Ha (electronic eigenvalue + offset)


class QuantumKrylovSolver:
    """Real-time quantum Krylov subspace diagonalisation on top of a MolecularHamiltonian."""

    def __init__(
        self,
        molecular_hamiltonian: MolecularHamiltonian,
        dt: Optional[float] = None,
        threshold: float = 1e-10,
        noise_sigma: float = 0.0,
        seed: Optional[int] = None,
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
        """
        self.mh = molecular_hamiltonian
        self.offset = molecular_hamiltonian.energy_offset
        self.threshold = threshold
        self.noise_sigma = float(noise_sigma)
        self._rng = np.random.default_rng(seed)

        self._H = molecular_hamiltonian.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
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
    def _default_dt(self) -> float:
        """dt = pi / (E_max - E_min); spectral extent estimated via sparse Lanczos."""
        dim = self._H.shape[0]
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
            self._basis.append(expm_multiply(-1j * self.dt * self._H, self._basis[-1]))

    def _subspace_matrices(self, dim: int):
        self._ensure_basis(dim)
        B = np.array(self._basis[:dim])              # (M, N)
        S = B.conj() @ B.T                           # <phi_i|phi_j>
        H = B.conj() @ self._H.dot(B.T)              # <phi_i|H|phi_j>
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
