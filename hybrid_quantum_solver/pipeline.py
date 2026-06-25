#!/usr/bin/env python3
"""
pipeline.py -- the live hybrid quantum-classical pipeline (Krylov Quantum Solver).

Wires the corrected components together:

    classical integrals  ->  molecular_hamiltonian (vetted Jordan-Wigner)
                          ->  quantum_krylov_solver  (real-time Krylov + optional shot noise)

This SUPERSEDES ``EnterprisePipelineOrchestrator`` in ``orchestrate_hybrid_pipeline.py``,
which used the broken ``AdvancedStochasticCompactor`` mapping, the near-identity "Krylov"
sampler, the asymmetric-Gaussian-on-H "noise", and the QCIVET symmetry stamp. See
``REFACTOR_PLAN.md``.

Two entry points:
  * ``run_geometry`` -- from a molecular geometry string (uses Qiskit Nature's PySCFDriver).
  * ``run_from_integrals`` -- from precomputed CASCI active-space integrals (the CIF/materials
    path via ``chemistry_gateway``); this is how NbN-style systems flow through the solver.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import (
    MolecularHamiltonian,
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
)
from hybrid_quantum_solver.noise import shot_noise_sigma
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

# Dense exact-reference diagonalisation is only attempted below this qubit count
# (2**14 = 16384 -> ~4 GB dense complex matrix). Larger systems must supply a reference.
_EXACT_REFERENCE_MAX_QUBITS = 12


@dataclass
class PipelineResult:
    computed_energy: float
    hf_energy: float
    reference_energy: Optional[float]            # exact diag / CASCI total (same frame), if known
    error_vs_reference: Optional[float]
    krylov_dim: int
    rank: int
    n_qubits: int
    dt: float
    shots: Optional[int]
    convergence: List[Tuple[int, int, float]] = field(default_factory=list)  # (M, rank, energy)

    def summary(self) -> str:
        ref = "n/a" if self.reference_energy is None else f"{self.reference_energy:.6f}"
        err = "n/a" if self.error_vs_reference is None else f"{self.error_vs_reference:.2e}"
        return (
            f"E={self.computed_energy:.6f} Ha  (HF={self.hf_energy:.6f}, ref={ref}, "
            f"|err|={err} Ha)  M={self.krylov_dim} rank={self.rank} "
            f"qubits={self.n_qubits} dt={self.dt:.4f} shots={self.shots}"
        )


def _run(
    mh: MolecularHamiltonian,
    krylov_dim: int,
    dt: Optional[float],
    shots: Optional[int],
    seed: Optional[int],
    reference_energy: Optional[float],
    track_convergence: bool,
) -> PipelineResult:
    sigma = shot_noise_sigma(shots) if shots else 0.0
    solver = QuantumKrylovSolver(mh, dt=dt, noise_sigma=sigma, seed=seed)

    if track_convergence:
        steps = solver.convergence(krylov_dim)
    else:
        steps = [solver.solve(krylov_dim)]
    final = steps[-1]

    err = None if reference_energy is None else abs(final.energy - reference_energy)
    return PipelineResult(
        computed_energy=final.energy,
        hf_energy=mh.hf_energy,
        reference_energy=reference_energy,
        error_vs_reference=err,
        krylov_dim=final.dim,
        rank=final.rank,
        n_qubits=mh.num_qubits,
        dt=solver.dt,
        shots=shots,
        convergence=[(s.dim, s.rank, s.energy) for s in steps],
    )


def run_geometry(
    atom: str,
    basis: str = "sto3g",
    charge: int = 0,
    spin: int = 0,
    active_electrons: Optional[int] = None,
    active_orbitals: Optional[int] = None,
    krylov_dim: int = 8,
    dt: Optional[float] = None,
    shots: Optional[int] = None,
    seed: Optional[int] = None,
    reference: str = "exact",
    track_convergence: bool = True,
) -> PipelineResult:
    """Run the solver from a molecular geometry string (PySCFDriver under the hood).

    ``reference="exact"`` adds a dense-diagonalisation reference when the qubit count is small
    enough; otherwise no reference is computed.
    """
    mh = build_molecular_hamiltonian(
        atom=atom, basis=basis, charge=charge, spin=spin,
        active_electrons=active_electrons, active_orbitals=active_orbitals,
    )
    reference_energy = None
    if reference == "exact" and mh.num_qubits <= _EXACT_REFERENCE_MAX_QUBITS:
        reference_energy = mh.ground_state_energy()
    return _run(mh, krylov_dim, dt, shots, seed, reference_energy, track_convergence)


def run_from_integrals(
    h1: np.ndarray,
    eri: np.ndarray,
    num_particles: Tuple[int, int],
    e_core: float,
    krylov_dim: int = 8,
    dt: Optional[float] = None,
    shots: Optional[int] = None,
    seed: Optional[int] = None,
    reference_energy: Optional[float] = None,
    track_convergence: bool = True,
) -> PipelineResult:
    """Run the solver from precomputed CASCI active-space integrals (CIF/materials path).

    Pass ``reference_energy = cas.e_tot`` (the active-space CASCI total) as the internal
    convergence target -- in a complete active space the Krylov estimate should approach it.
    """
    mh = build_hamiltonian_from_integrals(
        h1, eri, num_particles=num_particles, energy_offset=e_core
    )
    return _run(mh, krylov_dim, dt, shots, seed, reference_energy, track_convergence)
