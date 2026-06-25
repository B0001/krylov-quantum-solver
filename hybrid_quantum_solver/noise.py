#!/usr/bin/env python3
"""
noise.py -- legitimate noise modelling for the quantum Krylov pipeline.

Replaces two pieces of the original code that were not real noise handling:

  * ``h_matrix += np.random.normal(0, noise_variance, h_matrix.shape)`` -- asymmetric
    Gaussian numbers added straight onto the Hamiltonian. That breaks Hermiticity, makes
    the generalised eigenproblem unbounded below, and produced energies hundreds of Ha
    *below* the true ground state. (It also used ``noise_variance`` as a standard deviation.)
  * "QCIVET" -- a SHA-256 stamp plus a ``max|A - A^T|`` symmetry check that certified
    collapsed zero-energy runs as ``VERIFIED_SAFE`` and rejected correct physics over a
    1e-5 asymmetry. The hash verified nothing about the returned matrices.

What this module provides instead:

  * ``shot_noise_sigma(shots)`` -- the statistical scale (~1/sqrt(shots)) of a Hadamard-test
    estimate of a subspace matrix element, the noise the Krylov solver actually faces.
  * ``hermitize(A)`` -- enforce the Hermitian structure the physics guarantees. This is the
    legitimate "integrity" step: quantum subspace methods are noise-resilient through
    Hermitisation + singular-value thresholding, not through cryptographic hashing.
  * ``build_depolarizing_noise_model(...)`` -- a real qiskit-aer ``NoiseModel`` (depolarizing
    gate error + symmetric readout error) for device-level simulation of the evolution
    circuits, i.e. the hardware path for which the exact-statevector solver is the target.
"""
from __future__ import annotations

import numpy as np


def shot_noise_sigma(shots: int) -> float:
    """Statistical standard deviation of a Hadamard-test estimate from ``shots`` samples."""
    if shots <= 0:
        raise ValueError("shots must be a positive integer")
    return 1.0 / np.sqrt(float(shots))


def hermitize(matrix: np.ndarray) -> np.ndarray:
    """Project a matrix onto its Hermitian part: ``(A + A^dagger) / 2``."""
    a = np.asarray(matrix)
    return 0.5 * (a + a.conj().T)


def fold_global_circuit(circuit, scale_factor: int):
    """Global unitary folding for zero-noise extrapolation: C -> C (C^dagger C)^n.

    For an odd ``scale_factor = 2n+1`` the ideal operation is unchanged (C^dagger C = I) while
    the gate noise is amplified by ~scale_factor. ``scale_factor = 1`` returns the circuit as-is.
    The circuit must contain no measurements (state-prep only), as in the Hadamard-test pairs.
    """
    if scale_factor == 1:
        return circuit.copy()
    if scale_factor < 1 or scale_factor % 2 == 0:
        raise ValueError("scale_factor must be an odd integer >= 1 (1, 3, 5, ...)")
    inverse = circuit.inverse()
    folded = circuit.copy()
    for _ in range((scale_factor - 1) // 2):
        folded = folded.compose(inverse).compose(circuit)
    return folded


def extrapolate_zero_noise(scale_factors, values, order: int = 1) -> float:
    """Extrapolate measured ``values`` at ``scale_factors`` back to the zero-noise limit.

    A least-squares polynomial of the given ``order`` (1 = linear Richardson) evaluated at 0.
    """
    coeffs = np.polyfit(np.asarray(scale_factors, dtype=float),
                        np.asarray(values, dtype=float), order)
    return float(np.polyval(coeffs, 0.0))


def build_depolarizing_noise_model(
    one_qubit_error: float = 1e-3,
    two_qubit_error: float = 1e-2,
    readout_error: float = 1e-2,
):
    """Construct a simple but real qiskit-aer ``NoiseModel``.

    Depolarizing error on common 1- and 2-qubit gates plus a symmetric readout error.
    Requires ``qiskit-aer``. Use it to simulate the Trotterised evolution circuits at the
    device level (the noiseless statevector solver is the reference these should converge to).
    """
    from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error

    model = NoiseModel()
    model.add_all_qubit_quantum_error(
        depolarizing_error(one_qubit_error, 1),
        ["u1", "u2", "u3", "rz", "rx", "ry", "sx", "x", "h"],
    )
    model.add_all_qubit_quantum_error(
        depolarizing_error(two_qubit_error, 2), ["cx", "cz"]
    )
    p = readout_error
    model.add_all_qubit_readout_error(ReadoutError([[1 - p, p], [p, 1 - p]]))
    return model
