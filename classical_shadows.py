#!/usr/bin/env python3
"""
Classical shadows -- estimate energies and Pauli observables from randomized single-qubit
measurements (Huang, Kueng & Preskill, "Predicting many properties of a quantum system from very few
measurements", Nature Physics 16, 1050, 2020; arXiv:2002.08953).

A different primitive from the energy *methods* in this repo: instead of a tailored measurement per
observable, measure each qubit in a uniformly random Pauli basis (X/Y/Z) and reconstruct any Pauli
expectation classically. For a single random-Pauli snapshot, the unbiased single-qubit estimator of
a Pauli P_q is ``3 * s_q`` when the measured basis matches P_q's (s_q = +/-1 the outcome) and 0
otherwise; a Pauli string's estimate is the product over its support, and the energy estimate is the
coefficient-weighted sum over the Hamiltonian's Pauli terms. Averaging over snapshots converges to
the exact expectation.

The catch (the honest finding): the single-shot variance is bounded by the **shadow norm**
``sum_k |c_k|^2 3^{w_k}`` with ``w_k`` the Pauli weight -- the ``3^{w_k}`` factor makes high-weight
terms sample-expensive (an XXYY-type weight-4 term costs 81x its coefficient). Grouped / derandomized
shadows mitigate this and are out of scope.

Reuses the validated qubit Hamiltonian (``SparsePauliOp``) and any statevector (HF, FCI, ...). The
reference is the exact expectation ``<psi|H|psi>``. Exact statevector simulation of the measurement;
small systems.
"""
from __future__ import annotations

import numpy as np

# Single-qubit unitaries that rotate the measurement basis onto Z (so outcomes are sampled in Z).
_HAD = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
_SDG = np.array([[1, 0], [0, -1j]], dtype=complex)
_BASIS_ROT = {0: _HAD, 1: _HAD @ _SDG, 2: np.eye(2, dtype=complex)}   # 0=X, 1=Y, 2=Z


def _apply_1q(psi: np.ndarray, U: np.ndarray, q: int, n: int) -> np.ndarray:
    """Apply a 2x2 gate to qubit ``q`` of an ``n``-qubit statevector (qiskit ordering: qubit 0 LSB)."""
    t = psi.reshape([2] * n)
    t = np.tensordot(U, t, axes=([1], [n - 1 - q]))
    return np.moveaxis(t, 0, n - 1 - q).reshape(-1)


def collect_classical_shadow(statevector: np.ndarray, n_qubits: int, n_shots: int, seed: int = 0):
    """Collect ``n_shots`` random-Pauli snapshots of ``statevector``.

    Returns ``(bases, signs)``, each ``(n_shots, n_qubits)``: ``bases`` in {0,1,2} (the measured
    X/Y/Z basis per qubit) and ``signs`` in {+1,-1} (the measurement outcome).
    """
    rng = np.random.default_rng(seed)
    psi = np.asarray(statevector, dtype=complex)
    bases = rng.integers(0, 3, size=(n_shots, n_qubits))
    signs = np.empty((n_shots, n_qubits), dtype=int)
    for shot in range(n_shots):
        rot = psi
        for q in range(n_qubits):
            rot = _apply_1q(rot, _BASIS_ROT[bases[shot, q]], q, n_qubits)
        probs = np.abs(rot) ** 2
        probs /= probs.sum()
        idx = rng.choice(probs.size, p=probs)
        signs[shot] = [1 - 2 * ((idx >> q) & 1) for q in range(n_qubits)]
    return bases, signs


def _pauli_basis(pauli) -> np.ndarray:
    """Per-qubit measurement basis a Pauli term needs: 0=X, 1=Y, 2=Z, -1=identity (no constraint)."""
    n = pauli.num_qubits
    pb = np.full(n, -1, dtype=int)
    for q in range(n):
        x, z = bool(pauli.x[q]), bool(pauli.z[q])
        if x and not z:
            pb[q] = 0
        elif x and z:
            pb[q] = 1
        elif z and not x:
            pb[q] = 2
    return pb


def shadow_energy_samples(bases: np.ndarray, signs: np.ndarray, pauli_op) -> np.ndarray:
    """Per-snapshot energy estimates for ``pauli_op`` (a ``SparsePauliOp``); the mean estimates ``<H>``.

    Returned in the operator's own (electronic) frame -- add any energy offset separately.
    """
    per_shot = np.zeros(bases.shape[0])
    for coeff, pauli in zip(np.real(pauli_op.coeffs), pauli_op.paulis):
        pb = _pauli_basis(pauli)
        support = np.flatnonzero(pb >= 0)
        if support.size == 0:                                  # identity term
            per_shot += coeff
            continue
        want = pb[support]
        match = np.all(bases[:, support] == want, axis=1)
        value = np.prod(3 * signs[:, support], axis=1)
        per_shot += coeff * match * value
    return per_shot


def shadow_norm(pauli_op) -> float:
    """HKP shadow norm ``sum_k |c_k|^2 3^{w_k}`` -- an upper bound on the single-shot energy variance.

    ``w_k`` is the Pauli weight (number of non-identity qubits); the ``3^{w_k}`` growth is why
    high-weight observables are sample-expensive under random-Pauli shadows.
    """
    weights = np.array([int(np.sum(p.x | p.z)) for p in pauli_op.paulis])
    return float(np.sum(np.abs(pauli_op.coeffs) ** 2 * 3.0 ** weights))


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    op = mh.qubit_hamiltonian
    psi = np.asarray(mh.hf_state().data, dtype=complex)
    exact = float(np.real(psi.conj() @ op.to_matrix() @ psi))
    for shots in (1000, 4000, 16000):
        bases, signs = collect_classical_shadow(psi, mh.num_qubits, shots, seed=1)
        s = shadow_energy_samples(bases, signs, op)
        print(f"shots={shots:>6}: <H>={s.mean():+.5f} (exact {exact:+.5f})  "
              f"|err|={abs(s.mean()-exact):.4f}  stderr={s.std()/np.sqrt(shots):.4f}")
    print(f"shadow norm = {shadow_norm(op):.3f}  (variance bound; 3^w growth with Pauli weight)")
