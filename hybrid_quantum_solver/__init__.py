"""Krylov Quantum Solver -- corrected hybrid quantum-classical chemistry pipeline.

Public API (the validated, live path):
    build_molecular_hamiltonian      -- geometry -> qubit Hamiltonian (vetted Jordan-Wigner)
    build_hamiltonian_from_integrals -- CASCI active-space integrals -> qubit Hamiltonian
    QuantumKrylovSolver              -- real-time quantum Krylov subspace diagonalisation
    run_geometry / run_from_integrals-- end-to-end pipeline

The original ``EnterprisePipelineOrchestrator`` (orchestrate_hybrid_pipeline.py) and
``QiskitKrylovSampler`` (quantum_sampler.py) are retained only as quarantined regression
fixtures -- see REFACTOR_PLAN.md. They are intentionally NOT exported here.
"""
from .molecular_hamiltonian import (
    MolecularHamiltonian,
    build_dipole_operators,
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
)
from .pipeline import PipelineResult, run_from_integrals, run_geometry
from .qksd_properties import oscillator_strengths, property_matrix, transition_dipoles
from .quantum_krylov_solver import ExcitedKrylovStep, KrylovStep, QuantumKrylovSolver
from .trotter_krylov import TrotterKrylovSolver, build_trotter_step, estimate_energy_aer
from .hardware_krylov import HardwareKrylovSolver
from .dmrg_reference import dmrg_available, dmrg_energy, fci_energy, reference_energy

__all__ = [
    "MolecularHamiltonian",
    "build_molecular_hamiltonian",
    "build_hamiltonian_from_integrals",
    "build_dipole_operators",
    "QuantumKrylovSolver",
    "KrylovStep",
    "ExcitedKrylovStep",
    "property_matrix",
    "transition_dipoles",
    "oscillator_strengths",
    "TrotterKrylovSolver",
    "build_trotter_step",
    "estimate_energy_aer",
    "HardwareKrylovSolver",
    "fci_energy",
    "dmrg_energy",
    "reference_energy",
    "dmrg_available",
    "PipelineResult",
    "run_geometry",
    "run_from_integrals",
]
