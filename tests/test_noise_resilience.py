#!/usr/bin/env python3
"""
Phase 3 validation gate (see REFACTOR_PLAN.md).

Replaces the original "noise" path (asymmetric Gaussian added straight onto H, which was
unbounded below and produced energies hundreds of Ha under the true minimum) and the QCIVET
stamp. Under a *real* finite-sampling shot-noise model with Hermitian perturbations and a
noise-aware overlap cutoff, the quantum Krylov estimate must:

  1. be EXACT in the noiseless limit;
  2. DEGRADE GRACEFULLY -- more shots => smaller error;
  3. stay BOUNDED -- nothing remotely like the old -800 Ha blow-ups.

Also checks the legitimate replacements for QCIVET: Hermitisation and a real Aer noise model.

Run:  pytest tests/test_noise_resilience.py -v
"""
import numpy as np
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.noise import (
    extrapolate_zero_noise,
    fold_global_circuit,
    hermitize,
    shot_noise_sigma,
)
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

H2 = dict(atom="H 0 0 0; H 0 0 0.74")


def _energies(mh, shots, krylov_dim=6, seeds=6):
    sigma = 0.0 if shots is None else shot_noise_sigma(shots)
    return [QuantumKrylovSolver(mh, noise_sigma=sigma, seed=s).solve(krylov_dim).energy
            for s in range(seeds)]


def test_noiseless_is_exact():
    mh = build_molecular_hamiltonian(**H2)
    assert abs(QuantumKrylovSolver(mh).solve(6).energy - mh.ground_state_energy()) < 1e-6


def test_shot_noise_is_bounded_and_improves_with_shots():
    mh = build_molecular_hamiltonian(**H2)
    fci = mh.ground_state_energy()

    e_few = _energies(mh, shots=1024)
    e_many = _energies(mh, shots=16384)

    mean_err_few = np.mean([abs(e - fci) for e in e_few])
    mean_err_many = np.mean([abs(e - fci) for e in e_many])

    # 2. graceful degradation: 16x more shots -> clearly smaller mean error
    assert mean_err_many < mean_err_few

    # 3. bounded: even at a modest shot budget the error is chemical-scale, NOT -800 Ha
    assert max(abs(e - fci) for e in e_few) < 0.1
    # controlled: never dips far below the variational floor
    assert min(e_few + e_many) > fci - 0.1


def test_hermitize_projects_to_hermitian():
    a = np.array([[1.0, 2.0 + 1.0j], [0.0, 3.0]])
    h = hermitize(a)
    assert np.allclose(h, h.conj().T)


def test_shot_noise_sigma_scales():
    assert shot_noise_sigma(10000) == pytest.approx(0.01)
    with pytest.raises(ValueError):
        shot_noise_sigma(0)


def test_aer_noise_model_builds():
    pytest.importorskip("qiskit_aer")
    from qiskit_aer.noise import NoiseModel
    from hybrid_quantum_solver.noise import build_depolarizing_noise_model

    nm = build_depolarizing_noise_model(1e-3, 1e-2, 1e-2)
    assert isinstance(nm, NoiseModel)
    assert len(nm.to_dict().get("errors", [])) > 0


def test_fold_global_preserves_ideal_state():
    """ZNE folding C -> C C^dagger C must leave the ideal (noiseless) state unchanged."""
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector

    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)
    qc.rz(0.37, 1)
    assert Statevector(fold_global_circuit(qc, 3)).equiv(Statevector(qc))
    assert fold_global_circuit(qc, 5).depth() > qc.depth()      # noise is amplified
    with pytest.raises(ValueError):
        fold_global_circuit(qc, 2)                              # must be odd


def test_extrapolate_zero_noise_linear():
    # points on the line y = 2 - 0.5 x  ->  intercept 2.0 at x = 0
    assert extrapolate_zero_noise([1, 3, 5], [1.5, 0.5, -0.5]) == pytest.approx(2.0)
