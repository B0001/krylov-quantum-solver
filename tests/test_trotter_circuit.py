#!/usr/bin/env python3
"""
Phase 3 (hardware path) validation gate (see REFACTOR_PLAN.md).

The Krylov basis built from real Suzuki-Trotter evolution circuits must reproduce FCI (up to
Trotter error) -- in contrast to the original "qDRIFT" sampler, whose single-Pauli infinitesimal
rotations collapsed the basis to rank 1. The qiskit-aer expectation path must run exactly,
under shot noise, and under a device NoiseModel, returning finite, bounded energies.

Run:  pytest tests/test_trotter_circuit.py -v
"""
import numpy as np
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.trotter_krylov import (
    TrotterKrylovSolver,
    build_trotter_step,
    estimate_energy_aer,
)

H2 = dict(atom="H 0 0 0; H 0 0 0.74")


def test_trotter_step_is_a_genuine_entangling_circuit():
    """The evolution step is a real multi-qubit circuit, not a single near-identity rotation."""
    mh = build_molecular_hamiltonian(**H2)
    step = build_trotter_step(mh.qubit_hamiltonian, dt=1.0, order=2, reps=2).decompose(reps=3)
    ops = step.count_ops()
    assert step.depth() > 10
    assert ops.get("cx", 0) > 0          # genuine entanglement (the old qDRIFT step had ~none)


def test_trotter_krylov_converges_to_fci():
    mh = build_molecular_hamiltonian(**H2)
    fci = mh.ground_state_energy()
    steps = TrotterKrylovSolver(mh, trotter_order=2, trotter_reps=2).convergence(6)
    energies = [s.energy for s in steps]

    assert abs(energies[0] - mh.hf_energy) < 1e-9        # M=1 is the Hartree-Fock reference
    assert min(abs(e - fci) for e in energies) < 1e-3    # reaches FCI within Trotter error
    assert min(energies) > fci - 1e-6                    # respects the variational floor


def test_aer_exact_expectation_matches_rhf():
    mh = build_molecular_hamiltonian(**H2)
    e = estimate_energy_aer(mh.hf_circuit, mh.qubit_hamiltonian, mh.energy_offset)
    assert abs(e - mh.hf_energy) < 1e-6


def test_aer_device_noise_is_finite_and_bounded():
    pytest.importorskip("qiskit_aer")
    from hybrid_quantum_solver.noise import build_depolarizing_noise_model

    mh = build_molecular_hamiltonian(**H2)
    nm = build_depolarizing_noise_model(2e-3, 2e-2, 1e-2)
    e = estimate_energy_aer(
        mh.hf_circuit, mh.qubit_hamiltonian, mh.energy_offset, noise_model=nm, shots=8192
    )
    assert np.isfinite(e)
    assert abs(e - mh.hf_energy) < 0.2     # perturbed by noise, but bounded (not -800 Ha)
