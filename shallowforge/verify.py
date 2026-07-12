"""Compile + verification harness (task 3) and baseline CX counting (task 4 headline).

Correctness gate (ADR, spec §7): a compiled Trotter step's unitary must match exact
``exp(-iH·dt)`` to fidelity >= 1 - 1e-8 per step on the simulable golden systems (T0-T2) BEFORE
any gate count is quoted. ``baseline_cx_per_step`` reproduces the first-order-Trotter naive
baseline (the ~6,500 CX/step N2 figure) that every future rung is measured against.
"""

from __future__ import annotations

from typing import Any

# Same basis set the repo's HardwareKrylovSolver.resource_report transpiles to, so CX counts
# are comparable to benchmark_resources.py.
_BASIS_GATES = ["cx", "u", "rz", "h", "x", "sx", "p", "s", "sdg"]


def compile_ir_to_circuit(ir, dt: float, *, order: int = 1, reps: int = 1):
    """Compile a TermStreamIR into a materialized Suzuki-Trotter step circuit exp(-iH·dt)."""
    from hybrid_quantum_solver.trotter_krylov import build_trotter_step

    return build_trotter_step(ir.to_sparse_pauli_op(), dt, order=order, reps=reps)


def count_cx(circuit) -> int:
    """Two-qubit (CX) gate count after transpiling to the reference basis."""
    from qiskit import transpile

    transpiled = transpile(circuit, basis_gates=_BASIS_GATES, optimization_level=0)
    return int(transpiled.count_ops().get("cx", 0))


def step_fidelity(
    ir, dt: float, *, order: int = 1, reps: int = 1, n_states: int = 6, seed: int = 0
) -> float:
    """State fidelity of the compiled step vs exact exp(-iH·dt), averaged over random states.

    Applies the compiled circuit to ``n_states`` seeded Haar-random states and compares against
    ``expm(-iH·dt)|ψ⟩``; returns the mean ``|⟨ψ_exact|ψ_compiled⟩|²``. This uses the exact
    Hamiltonian matrix once (``expm``) but only statevector *applications* of the circuit
    (O(gates·dim)), so it stays tractable at the T2 (LiH, 1024-dim) size where building the full
    compiled ``Operator`` would not. 1.0 iff the compiled step matches exact evolution (up to
    global phase) on the sampled subspace.
    """
    import numpy as np
    from qiskit.quantum_info import Statevector
    from scipy.linalg import expm

    op = ir.to_sparse_pauli_op()
    u_exact = expm(-1j * op.to_matrix() * dt)
    qc = compile_ir_to_circuit(ir, dt, order=order, reps=reps)
    dim = u_exact.shape[0]
    rng = np.random.default_rng(seed)

    total = 0.0
    for _ in range(n_states):
        psi = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        psi /= np.linalg.norm(psi)
        exact = u_exact @ psi
        compiled = Statevector(psi).evolve(qc).data
        total += float(abs(np.vdot(exact, compiled)) ** 2)
    return total / n_states


def baseline_cx_per_step(ir, dt: float) -> int:
    """The frozen baseline: first-order Trotter, one rep, naive term order. CX per step."""
    return count_cx(compile_ir_to_circuit(ir, dt, order=1, reps=1))


def compiled_totals(circuit) -> dict[str, Any]:
    """CX / depth / ancilla totals for a compiled circuit (feeds the manifest)."""
    from qiskit import transpile

    transpiled = transpile(circuit, basis_gates=_BASIS_GATES, optimization_level=0)
    return {
        "cx_per_step": int(transpiled.count_ops().get("cx", 0)),
        "depth": int(transpiled.depth()),
        "ancillas": 0,  # Trotter path uses no ancillas; qubitization (R4, M3+) will.
    }
