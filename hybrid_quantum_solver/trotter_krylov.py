#!/usr/bin/env python3
"""
trotter_krylov.py -- quantum Krylov with Trotterised real-time evolution circuits.

The exact-statevector ``QuantumKrylovSolver`` validates the *algorithm*. This module builds the
same Krylov space from actual quantum circuits: proper Suzuki-Trotter evolution of the FULL
Hamiltonian (``PauliEvolutionGate``). That is the correct replacement for the original code's
"qDRIFT" step, which applied a single Pauli rotation with an infinitesimal angle and produced a
rank-1 (collapsed) basis.

It also provides ``estimate_energy_aer`` -- an expectation-value path through qiskit-aer that can
run exactly, with finite-shot sampling, or under a real device ``NoiseModel`` -- so the hardware
behaviour can be studied. The noiseless statevector solver is the reference these circuits should
reproduce, up to Trotter error (which shrinks with more Trotter ``reps`` / higher ``order``).
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.synthesis import SuzukiTrotter
from scipy.sparse.linalg import eigsh

from hybrid_quantum_solver.certification import certify_energies
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import KrylovStep, ritz_pairs


def canonical_term_order(hamiltonian: SparsePauliOp) -> SparsePauliOp:
    """Reorder Pauli terms canonically: largest |coefficient| first, label as tiebreak.

    Suzuki-Trotter products depend on the term order, and the order coming out of the
    Jordan-Wigner/simplify pipeline follows Python's randomized hashing -- so without this,
    the synthesized unitary (and its Trotter error) differed BETWEEN PROCESSES (observed:
    Nb3F8 sector-2 reps=1 deviation 5.96e-3 or 1.04e-2 depending on the run, which coin-
    flipped a spec gate). Largest-first is the best of the orderings measured in
    specs/SPEC_trotter_resolution_floor.md and reproduces the historical green-run values.
    """
    labels = [p.to_label() for p in hamiltonian.paulis]
    coeffs = [complex(c) for c in hamiltonian.coeffs]
    idx = sorted(range(len(labels)), key=lambda i: (-abs(coeffs[i]), labels[i]))
    return SparsePauliOp.from_list([(labels[i], coeffs[i]) for i in idx])


def build_trotter_step(
    hamiltonian: SparsePauliOp, dt: float, order: int = 2, reps: int = 1
) -> QuantumCircuit:
    """One Trotter step exp(-i H dt) as a MATERIALIZED Suzuki-Trotter circuit.

    The synthesized definition is inlined instead of appending the opaque ``PauliEvolutionGate``:
    ``Operator()`` and ``Statevector.evolve()`` evaluate the opaque gate through its EXACT matrix,
    silently ignoring the attached synthesis -- which made every statevector-simulated "Trotter"
    study in this repo exact evolution in disguise (``order``/``reps`` were no-ops, and the old
    "within Trotter error" test passed vacuously; the qDRIFT/QCIVET failure mode in a new coat).
    Materializing guarantees every consumer -- statevector, ``Operator``, Aer, the
    ancilla-controlled hardware path -- sees the same genuinely Trotterized unitary.
    Regression gate: specs/SPEC_trotter_odmd.md G1 (``tests/test_trotter_odmd_spec.py``).

    Terms are canonically ordered first (see :func:`canonical_term_order`) so the circuit is
    deterministic across processes. Gate: specs/SPEC_trotter_resolution_floor.md G1.
    """
    gate = PauliEvolutionGate(
        canonical_term_order(hamiltonian), time=dt,
        synthesis=SuzukiTrotter(order=order, reps=reps)
    )
    qc = QuantumCircuit(hamiltonian.num_qubits)
    qc.compose(gate.definition, inplace=True)
    return qc


def estimate_energy_aer(
    state_circuit: QuantumCircuit,
    hamiltonian: SparsePauliOp,
    energy_offset: float = 0.0,
    noise_model=None,
    shots: Optional[int] = None,
) -> float:
    """Expectation <psi|H|psi> via qiskit-aer, returned in the total-energy frame.

    ``noise_model=None, shots=None`` -> exact statevector. ``shots`` sets finite-sampling
    precision (1/sqrt(shots)); ``noise_model`` injects device noise.
    """
    from qiskit_aer.primitives import EstimatorV2

    options = {}
    if noise_model is not None:
        options["backend_options"] = {"noise_model": noise_model}
    if shots is not None:
        options["default_precision"] = 1.0 / np.sqrt(float(shots))
    estimator = EstimatorV2(options=options) if options else EstimatorV2()
    ev = estimator.run([(state_circuit, hamiltonian)]).result()[0].data.evs
    return float(ev) + float(energy_offset)


class TrotterKrylovSolver:
    """Quantum Krylov subspace solver using Trotterised evolution circuits for the basis."""

    def __init__(
        self,
        molecular_hamiltonian: MolecularHamiltonian,
        dt: Optional[float] = None,
        trotter_order: int = 2,
        trotter_reps: int = 1,
        threshold: float = 1e-10,
    ):
        self.mh = molecular_hamiltonian
        self.offset = molecular_hamiltonian.energy_offset
        self.threshold = threshold
        self._Hop = molecular_hamiltonian.qubit_hamiltonian
        self._H = self._Hop.to_matrix(sparse=True).tocsc()
        self.dt = float(dt) if dt is not None else self._default_dt()
        self.step_circuit = build_trotter_step(self._Hop, self.dt, trotter_order, trotter_reps)
        self._basis: List[Statevector] = [Statevector(molecular_hamiltonian.hf_circuit)]

    def _default_dt(self) -> float:
        try:
            e_max = float(eigsh(self._H, k=1, which="LA", return_eigenvectors=False)[0])
            e_min = float(eigsh(self._H, k=1, which="SA", return_eigenvectors=False)[0])
            width = e_max - e_min
        except Exception:
            width = 2.0 * float(np.sum(np.abs(self._Hop.coeffs)))
        if not np.isfinite(width) or width <= 0:
            width = 2.0 * float(np.sum(np.abs(self._Hop.coeffs)))
        return np.pi / width

    def _ensure_basis(self, dim: int) -> None:
        while len(self._basis) < dim:
            self._basis.append(self._basis[-1].evolve(self.step_circuit))

    def _rows(self, dim: int) -> np.ndarray:
        """The first ``dim`` Krylov vectors as rows, shape (M, N)."""
        return np.array([s.data for s in self._basis[:dim]])

    def _steps(self, dims, certificate_dir, names) -> List[KrylovStep]:
        """Solve every dimension in ``dims``, then certify the energies in one call
        (``certification.certify_energies``: always, lowest energy first)."""
        ritz = []
        for m in dims:
            if m < 1:
                raise ValueError("krylov_dim must be >= 1")
            self._ensure_basis(m)
            B = self._rows(m)
            S = B.conj() @ B.T
            H = B.conj() @ self._H.dot(B.T)
            ritz.append(ritz_pairs(0.5 * (H + H.conj().T), 0.5 * (S + S.conj().T),
                                   self.threshold))
        certs = certify_energies(
            self.mh, [(vals[0], lambda C=C, m=m: C[:, 0] @ self._rows(m))
                      for m, (vals, C, _) in zip(dims, ritz)], certificate_dir, names)
        return [KrylovStep(dim=m, rank=rank, unverified_energy=float(vals[0]) + self.offset,
                           certificate=cert)
                for m, (vals, _, rank), cert in zip(dims, ritz, certs)]

    # Certified on every call, as in ``QuantumKrylovSolver``: read ``.energy`` for the
    # certified number and ``.unverified_energy`` for the raw estimate.
    def solve(self, krylov_dim: int, certificate_dir: str | Path | None = None) -> KrylovStep:
        """Ground-state Ritz energy with the checker's verdict on it."""
        return self._steps([krylov_dim], certificate_dir, None)[0]

    def convergence(self, max_dim: int,
                    certificate_dir: str | Path | None = None) -> List[KrylovStep]:
        """Steps for M = 1 .. max_dim; one checker run covers the sweep."""
        dims = range(1, max_dim + 1)
        return self._steps(dims, certificate_dir, [f"M{m}" for m in dims])


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    fci = mh.ground_state_energy()
    solver = TrotterKrylovSolver(mh, trotter_order=2, trotter_reps=2)
    print(f"H2 FCI={fci:.6f}  dt={solver.dt:.4f}  (2nd-order Trotter, 2 reps/step)")
    for s in solver.convergence(6):
        print(f"  M={s.dim} rank={s.rank} E={s.unverified_energy:.6f} "
              f"err={(s.unverified_energy - fci) * 1e3:+.4f} mHa  {s.certificate.line.split()[0]}")
    e_exact = estimate_energy_aer(mh.hf_circuit, mh.qubit_hamiltonian, mh.energy_offset)
    print(f"Aer exact HF expectation = {e_exact:.6f} (RHF={mh.hf_energy:.6f})")
