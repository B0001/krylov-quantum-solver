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

from typing import List, Optional

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.synthesis import SuzukiTrotter
from scipy.sparse.linalg import eigsh

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import KrylovStep, solve_generalized_eig


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

    def solve(self, krylov_dim: int) -> KrylovStep:
        if krylov_dim < 1:
            raise ValueError("krylov_dim must be >= 1")
        self._ensure_basis(krylov_dim)
        B = np.array([s.data for s in self._basis[:krylov_dim]])
        S = B.conj() @ B.T
        H = B.conj() @ self._H.dot(B.T)
        energy, rank = solve_generalized_eig(
            0.5 * (H + H.conj().T), 0.5 * (S + S.conj().T), self.threshold
        )
        return KrylovStep(dim=krylov_dim, rank=rank, energy=energy + self.offset)

    def convergence(self, max_dim: int) -> List[KrylovStep]:
        self._ensure_basis(max_dim)
        return [self.solve(m) for m in range(1, max_dim + 1)]


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    fci = mh.ground_state_energy()
    solver = TrotterKrylovSolver(mh, trotter_order=2, trotter_reps=2)
    print(f"H2 FCI={fci:.6f}  dt={solver.dt:.4f}  (2nd-order Trotter, 2 reps/step)")
    for s in solver.convergence(6):
        print(f"  M={s.dim} rank={s.rank} E={s.energy:.6f} err={(s.energy - fci) * 1e3:+.4f} mHa")
    e_exact = estimate_energy_aer(mh.hf_circuit, mh.qubit_hamiltonian, mh.energy_offset)
    print(f"Aer exact HF expectation = {e_exact:.6f} (RHF={mh.hf_energy:.6f})")
