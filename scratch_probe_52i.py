import numpy as np
from certified_gaps import gap_bracket, reachable_gap
from gap_selfcheck import self_checked_gap_from
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

GEOMS = {
    "linear H4 R=3.0": "H 0 0 0; H 0 0 3.0; H 0 0 6.0; H 0 0 9.0",
    "linear H4 R=2.0 (stretched)": "H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0",
    "asym H4 (1.0/2.0/4.0)": "H 0 0 0; H 0 0 1.0; H 0 0 3.0; H 0 0 7.0",
    "square H4 a=1.4": "H 0 0 0; H 1.4 0 0; H 1.4 1.4 0; H 0 1.4 0",
}

for name, atom in GEOMS.items():
    mh = build_molecular_hamiltonian(atom=atom)
    solver = QuantumKrylovSolver(mh)
    gap = reachable_gap(mh)
    print("=" * 80)
    print(f"{name}: exact reachable gap = {gap:.6f} Ha")
    for m in (6, 8, 12, 16):
        br = gap_bracket(mh, m, solver=solver)
        escaped = not (br.gap_lower - 1e-9 <= gap <= br.gap_upper + 1e-9)
        print(f"  M={m:2d}  Delta_lo={br.gap_lower:.6f}  Delta_hi={br.gap_upper:.6f}  "
              f"escaped={escaped}")

print("### extra asymmetric search ###")
GEOMS2 = {
    "asym B (0.9/2.0/5.0)": "H 0 0 0; H 0 0 0.9; H 0 0 2.9; H 0 0 7.9",
    "asym C (1.0/2.5/6.0)": "H 0 0 0; H 0 0 1.0; H 0 0 3.5; H 0 0 9.5",
    "asym D (1.0/2.0/5.0)": "H 0 0 0; H 0 0 1.0; H 0 0 3.0; H 0 0 8.0",
    "asym E (0.9/1.8/4.5) trimer-like": "H 0 0 0; H 0 0 0.9; H 0 0 2.7; H 0 0 7.2",
}
for name, atom in GEOMS2.items():
    mh = build_molecular_hamiltonian(atom=atom)
    solver = QuantumKrylovSolver(mh)
    gap = reachable_gap(mh)
    print("=" * 80)
    print(f"{name}: exact reachable gap = {gap:.6f} Ha")
    for m in (6, 8, 12, 16):
        br = gap_bracket(mh, m, solver=solver)
        escaped = not (br.gap_lower - 1e-9 <= gap <= br.gap_upper + 1e-9)
        print(f"  M={m:2d}  Delta_lo={br.gap_lower:.6f}  Delta_hi={br.gap_upper:.6f}  "
              f"escaped={escaped}")

print("### gap_selfcheck corroboration on escaping cases ###")
CHECK_GEOMS = {
    "linear H4 R=3.0": "H 0 0 0; H 0 0 3.0; H 0 0 6.0; H 0 0 9.0",
    "asym C (1.0/2.5/6.0)": "H 0 0 0; H 0 0 1.0; H 0 0 3.5; H 0 0 9.5",
}
for name, atom in CHECK_GEOMS.items():
    mh = build_molecular_hamiltonian(atom=atom)
    solver = QuantumKrylovSolver(mh)
    gap = reachable_gap(mh)
    dims = (6, 8, 12, 16)
    (lo, hi), flags = self_checked_gap_from(mh, dims, solver=solver)
    print(f"{name}: gap={gap:.6f}  self-checked=[{lo:.6f},{hi:.6f}] contains={lo<=gap<=hi}  flags={list(zip(dims, flags))}")
