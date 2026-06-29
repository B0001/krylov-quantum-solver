#!/usr/bin/env python3
"""
skqd.py -- sample-based Krylov quantum diagonalisation (SKQD).

Bridges the project's two near-term pillars: real-time quantum Krylov
(``quantum_krylov_solver.py``) and sample-based diagonalisation (the SQD path). Instead of
forming the dense Krylov matrix elements <phi_i|H|phi_j> (whose off-diagonal terms need
controlled time evolution / Hadamard tests -- see ``hardware_krylov.py``), SKQD *samples*
computational-basis determinants from the real-time-evolved states

    |Psi_k> = e^{-i k dt H} |HF>,    k = 0, 1, ..., depth-1,

unions the sampled determinants into a subspace D (the Hartree-Fock determinant is always
included), and diagonalises H projected onto D. This is selected-CI in the qubit
computational basis, with the selection driven by quantum real-time evolution.

Because the projected Hamiltonian is the exact H restricted to a subspace spanned by genuine
basis states, its lowest eigenvalue is a Rayleigh-Ritz estimate and is therefore
**variationally bounded below by the true ground-state energy** (it can never dip below FCI).
As the sample budget and Krylov depth grow, the determinant subspace fills out the
ground-state support and the energy converges to FCI from above (the SKQD convergence
guarantee of Yu et al. 2025, arXiv:2501.09702; the qDRIFT-compiled variant SqDRIFT is
Piccinelli et al., arXiv:2508.02578).

Scope (see specs/SPEC_skqd.md): time evolution is exact on the statevector and sampling is
from the exact |Psi_k> amplitudes -- this isolates and validates the *subspace* algorithm.
qDRIFT/Trotter circuit compilation and hardware-noise sampling are out of scope here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver


@dataclass
class SKQDStep:
    """Result of one sample-based Krylov diagonalisation."""
    depth: int        # number of Krylov states sampled (k = 0 .. depth-1)
    n_dets: int       # determinants in the diagonalisation subspace (incl. HF)
    energy: float     # total energy in Ha (subspace eigenvalue + offset)


class SampleKrylovSolver:
    """Sample-based Krylov quantum diagonalisation on top of a ``MolecularHamiltonian``."""

    def __init__(
        self,
        molecular_hamiltonian: MolecularHamiltonian,
        dt: Optional[float] = None,
        n_shots: int = 10_000,
        depth: int = 10,
        seed: Optional[int] = None,
    ):
        """
        Args:
            molecular_hamiltonian: built by ``build_molecular_hamiltonian`` / from integrals.
            dt: real-time step for the Krylov states. ``None`` uses the validated default of
                ``QuantumKrylovSolver`` (pi / spectral width). Pass the same ``dt`` as an
                exact ``QuantumKrylovSolver`` to compare the two on a matched Krylov space.
            n_shots: determinants sampled per Krylov state in ``solve`` (the HF determinant is
                always added, so the subspace is never empty).
            depth: number of Krylov states |Psi_k>, k = 0 .. depth-1.
            seed: RNG seed for reproducible sampling.
        """
        if depth < 1:
            raise ValueError("depth must be >= 1")
        self.mh = molecular_hamiltonian
        self.offset = molecular_hamiltonian.energy_offset
        self.n_shots = int(n_shots)
        self.depth = int(depth)
        self._rng = np.random.default_rng(seed)

        # Reuse the validated real-time propagation (and the robust default dt) verbatim.
        self._krylov = QuantumKrylovSolver(molecular_hamiltonian, dt=dt)
        self.dt = self._krylov.dt
        self._H = self._krylov._H.tocsr()

        # The Hartree-Fock determinant is a single computational basis state.
        hf = np.asarray(molecular_hamiltonian.hf_state().data)
        self._hf_index = int(np.argmax(np.abs(hf)))

    # -- sampling ------------------------------------------------------------
    def _state_probabilities(self, k: int) -> np.ndarray:
        """Born-rule distribution |<b|Psi_k>|^2 over computational basis states."""
        psi = np.asarray(self._krylov._basis[k])
        p = np.abs(psi) ** 2
        total = p.sum()
        return p / total if total > 0 else p

    def _draw(self, k: int, n_shots: int) -> np.ndarray:
        """Sample ``n_shots`` basis-state indices from |Psi_k>."""
        p = self._state_probabilities(k)
        return self._rng.choice(p.size, size=n_shots, p=p)

    # -- diagonalisation -----------------------------------------------------
    def _diagonalize(self, dets: Iterable[int]) -> SKQDStep:
        """Lowest eigenvalue of H restricted to the determinant subspace ``dets``."""
        idx = np.fromiter(sorted(dets), dtype=int)
        H_sub = self._H[idx][:, idx].toarray()
        H_sub = 0.5 * (H_sub + H_sub.conj().T)          # exact H is Hermitian; guard rounding
        energy = float(np.linalg.eigvalsh(H_sub)[0].real)
        return SKQDStep(depth=self.depth, n_dets=idx.size, energy=energy + self.offset)

    # -- public API ----------------------------------------------------------
    def solve(self) -> SKQDStep:
        """Sample ``n_shots`` per Krylov state, union determinants, diagonalise."""
        self._krylov._ensure_basis(self.depth)
        dets = {self._hf_index}
        for k in range(self.depth):
            dets.update(int(b) for b in np.unique(self._draw(k, self.n_shots)))
        return self._diagonalize(dets)

    def convergence(self, shot_schedule: Iterable[int]) -> List[SKQDStep]:
        """Energies for an increasing shot budget at fixed depth.

        The draws are *nested*: the largest budget is sampled once and smaller budgets read a
        prefix, so a larger budget's determinant set is a superset of a smaller one's. This
        makes the subspaces nested, hence ``n_dets`` non-decreasing and ``energy``
        non-increasing in the schedule (a clean monotone-convergence check).
        """
        schedule = sorted(int(s) for s in shot_schedule)
        if not schedule:
            return []
        self._krylov._ensure_basis(self.depth)
        max_shots = schedule[-1]
        draws = [self._draw(k, max_shots) for k in range(self.depth)]
        steps = []
        for s in schedule:
            dets = {self._hf_index}
            for d in draws:
                dets.update(int(b) for b in np.unique(d[:s]))
            steps.append(self._diagonalize(dets))
        return steps


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")
    e_fci = mh.ground_state_energy()
    print(f"H4 chain: FCI={e_fci:.6f} Ha")
    for step in SampleKrylovSolver(mh, depth=8, seed=0).convergence((2_000, 10_000, 50_000)):
        print(f"  n_dets={step.n_dets:3d}  E={step.energy:.6f}  "
              f"err={(step.energy - e_fci) * 1e3:+.4f} mHa")
