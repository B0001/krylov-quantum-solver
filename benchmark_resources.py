#!/usr/bin/env python3
"""
Resource accounting for the on-hardware quantum Krylov solver.

Tabulates the honest quantum cost of the Hadamard-test solver (hardware_krylov.py) as a function
of Krylov dimension M: qubit count, Trotter-step depth and 2-qubit-gate count, the deepest
Hadamard-test circuit, the number of distinct circuits (M(M+1)/2), observable evaluations, and the
total shot budget. This replaces the previous code's fictitious wall-clock "benchmarks" (e.g. a
19-hour 16-qubit job) with reproducible circuit-level resource estimates.

Run:  python benchmark_resources.py   ->  prints a resource table and writes data/krylov_resources.csv
"""
import csv

from hybrid_quantum_solver.hardware_krylov import HardwareKrylovSolver
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

SHOTS = 8192
KRYLOV_DIMS = [2, 4, 6]
SYSTEMS = {
    "H2 (STO-3G)": dict(atom="H 0 0 0; H 0 0 0.74"),
    "N2 CAS(6,6)": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
FIELDS = ["system", "krylov_dim", "qubits", "hamiltonian_pauli_terms",
          "trotter_step_depth", "trotter_step_cx", "deepest_circuit_depth", "deepest_circuit_cx",
          "distinct_pair_circuits", "observable_evaluations", "total_shots"]


def main():
    rows = []
    header = (f"{'system':14s} {'M':>2} {'qb':>3} {'Hterms':>6} {'step_d':>7} {'step_cx':>7} "
              f"{'deep_d':>7} {'deep_cx':>7} {'circs':>6} {'obs_ev':>7} {'shots':>9}")
    print(header)
    print("-" * len(header))
    for name, spec in SYSTEMS.items():
        mh = build_molecular_hamiltonian(basis="sto3g", **spec)
        solver = HardwareKrylovSolver(mh)
        for m in KRYLOV_DIMS:
            r = solver.resource_report(m, shots=SHOTS)
            print(f"{name:14s} {m:2d} {r['qubits']:3d} {r['hamiltonian_pauli_terms']:6d} "
                  f"{r['trotter_step_depth']:7d} {r['trotter_step_cx']:7d} "
                  f"{r['deepest_circuit_depth']:7d} {r['deepest_circuit_cx']:7d} "
                  f"{r['distinct_pair_circuits']:6d} {r['observable_evaluations']:7d} {r['total_shots']:9d}")
            rows.append({k: r[k] for k in FIELDS if k != "system"} | {"system": name})

    with open("data/krylov_resources.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print("\nstep_d/step_cx = single Trotter-step depth & CX count;  deep_* = deepest Hadamard-test "
          "circuit (i=0,j=M-1).\nobs_ev = 4 observables x M(M+1)/2 circuits;  shots = obs_ev x "
          f"{SHOTS}.  -> data/krylov_resources.csv")


if __name__ == "__main__":
    main()
