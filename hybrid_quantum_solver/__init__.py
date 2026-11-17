"""Krylov Quantum Solver -- corrected hybrid quantum-classical chemistry pipeline.

Public API (the validated, live path):
    build_molecular_hamiltonian      -- geometry -> qubit Hamiltonian (vetted Jordan-Wigner)
    build_hamiltonian_from_integrals -- CASCI active-space integrals -> qubit Hamiltonian
    QuantumKrylovSolver              -- real-time quantum Krylov subspace diagonalisation
    run_geometry / run_from_integrals-- end-to-end pipeline

Every ground-state energy these return is certificate-or-abstention: ``.energy`` is the number
only if certkit's checker VERIFIED it out of process and raises ``EnergyAbstained`` otherwise;
``.unverified_energy`` is the raw estimate (specs/SPEC_certificate_or_abstention.md).

The original ``EnterprisePipelineOrchestrator`` (orchestrate_hybrid_pipeline.py) and
``QiskitKrylovSampler`` (quantum_sampler.py) are retained only as quarantined regression
fixtures -- see REFACTOR_PLAN.md. They are intentionally NOT exported here.
"""
from .certification import EnergyAbstained, Verdict
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
from .model_hamiltonians import (
    ModelIntegrals,
    hubbard_chain_integrals,
    hubbard_dimer_energy,
    hubbard_dimer_gap,
    hubbard_integrals,
    lieb_wu_energy,
    load_from_nb3x8_database,
)

__all__ = [
    "MolecularHamiltonian",
    "build_molecular_hamiltonian",
    "build_hamiltonian_from_integrals",
    "build_dipole_operators",
    "QuantumKrylovSolver",
    "KrylovStep",
    "ExcitedKrylovStep",
    "EnergyAbstained",
    "Verdict",
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
    "ModelIntegrals",
    "hubbard_integrals",
    "hubbard_chain_integrals",
    "hubbard_dimer_energy",
    "hubbard_dimer_gap",
    "lieb_wu_energy",
    "load_from_nb3x8_database",
]
