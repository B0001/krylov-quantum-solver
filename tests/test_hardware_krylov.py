#!/usr/bin/env python3
"""
On-hardware quantum Krylov validation gate (see REFACTOR_PLAN.md, Phase 3).

The subspace matrices Sᵢⱼ, Hᵢⱼ are MEASURED by ancilla Hadamard tests with controlled Trotter
evolution (qiskit-aer Estimator), not taken from statevector inner products. The estimate must:

  1. reproduce FCI in the exact (noiseless) limit, respecting the variational floor;
  2. stay finite and bounded under finite-shot sampling and a device NoiseModel
     (the original code's "noise" path returned values hundreds of Ha below the true minimum).

Run:  pytest tests/test_hardware_krylov.py -v
"""
import numpy as np
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.hardware_krylov import HardwareKrylovSolver

pytest.importorskip("qiskit_aer")

H2 = dict(atom="H 0 0 0; H 0 0 0.74")


def test_hardware_krylov_matches_fci_exact():
    mh = build_molecular_hamiltonian(**H2)
    fci = mh.ground_state_energy()
    steps = HardwareKrylovSolver(mh).convergence(4)
    energies = [s.energy for s in steps]

    assert abs(energies[0] - mh.hf_energy) < 1e-6        # M=1 is the measured HF energy
    assert min(abs(e - fci) for e in energies) < 1e-3    # reaches FCI
    assert min(energies) > fci - 1e-6                    # variational floor


def test_hardware_krylov_under_shot_noise_is_bounded():
    mh = build_molecular_hamiltonian(**H2)
    fci = mh.ground_state_energy()
    energies = [HardwareKrylovSolver(mh, shots=8192, seed=s).solve(3).energy for s in range(3)]
    assert all(np.isfinite(e) for e in energies)
    assert max(abs(e - fci) for e in energies) < 0.05    # bounded (mHa-scale), not -800 Ha


def test_hardware_krylov_under_device_noise_runs():
    from hybrid_quantum_solver.noise import build_depolarizing_noise_model

    mh = build_molecular_hamiltonian(**H2)
    fci = mh.ground_state_energy()
    nm = build_depolarizing_noise_model(1e-3, 1e-2, 1e-2)
    energy = HardwareKrylovSolver(mh, shots=4096, noise_model=nm, seed=0).solve(2).energy
    assert np.isfinite(energy)
    assert abs(energy - fci) < 0.2                        # perturbed by device noise, still bounded


def test_zne_reduces_device_noise_error():
    """Zero-noise extrapolation lowers the device-noise energy error (good-hardware regime)."""
    from hybrid_quantum_solver.noise import build_depolarizing_noise_model

    mh = build_molecular_hamiltonian(**H2)
    fci = mh.ground_state_energy()
    # Gate noise only (global folding amplifies gate noise, not readout); good 2-qubit error rate.
    nm = build_depolarizing_noise_model(5e-5, 5e-4, 0.0)

    err_noisy = abs(HardwareKrylovSolver(mh, noise_model=nm).solve(2).energy - fci)
    err_zne = abs(
        HardwareKrylovSolver(mh, noise_model=nm, zne_scale_factors=[1, 3, 5]).solve(2).energy - fci
    )

    assert err_noisy > 1e-3        # device noise genuinely present (~9 mHa)
    assert err_zne < err_noisy     # ZNE reduces the error
    assert err_zne < 0.01          # mitigated to single-digit mHa


def test_resource_report_scales_sensibly():
    mh = build_molecular_hamiltonian(**H2)
    solver = HardwareKrylovSolver(mh)
    r2 = solver.resource_report(2, shots=1000)
    r4 = solver.resource_report(4, shots=1000)

    assert r2["qubits"] == mh.num_qubits + 1                      # +1 ancilla
    assert r2["distinct_pair_circuits"] == 3 and r4["distinct_pair_circuits"] == 10  # M(M+1)/2
    assert r4["deepest_circuit_cx"] > r2["deepest_circuit_cx"]    # deeper for larger M
    assert r2["total_shots"] == r2["observable_evaluations"] * 1000

    # ZNE triples the number of observable evaluations
    rz = HardwareKrylovSolver(mh, zne_scale_factors=[1, 3, 5]).resource_report(2, shots=1000)
    assert rz["observable_evaluations"] == 3 * r2["observable_evaluations"]
