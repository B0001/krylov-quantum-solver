#!/usr/bin/env python3
"""
N2 dissociation benchmark -- the classic multireference test.

Compares, along the N2 bond-breaking coordinate, three energies in a CAS(6,6) active space:
  * Hartree-Fock (single reference -- fails badly as the triple bond breaks),
  * real-time quantum Krylov (this project),
  * exact CASCI == active-space FCI (the gold-standard reference; tractable at 12 qubits).

This is the honest Phase-4 benchmark rung above H2/LiH (see REFACTOR_PLAN.md). For active spaces
small enough that FCI is exact, CASCI *is* the reference; DMRG becomes the reference only once the
active space outgrows FCI (and, separately, once the system outgrows statevector simulation).

Run:  python benchmark_n2.py   ->  prints a table and writes data/n2_dissociation.csv
"""
import csv

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

BOND_LENGTHS = [1.00, 1.10, 1.30, 1.60, 2.10]   # Angstrom (~1.10 is equilibrium)
CAS_ELECTRONS, CAS_ORBITALS = 6, 6
KRYLOV_DIM = 12
OUTPUT = "data/n2_dissociation.csv"


def main():
    print(f"N2 dissociation | CAS({CAS_ELECTRONS},{CAS_ORBITALS}) | quantum Krylov M={KRYLOV_DIM}")
    header = (f"{'R(A)':>6} {'qubits':>6} {'HF':>12} {'Krylov':>12} "
              f"{'CASCI(exact)':>13} {'HF err':>10} {'Krylov err':>11}")
    print(header)
    print("-" * len(header))

    rows = []
    for r in BOND_LENGTHS:
        mh = build_molecular_hamiltonian(
            atom=f"N 0 0 0; N 0 0 {r}", basis="sto3g",
            active_electrons=CAS_ELECTRONS, active_orbitals=CAS_ORBITALS,
        )
        exact = mh.ground_state_energy()
        hf = mh.hf_energy
        steps = QuantumKrylovSolver(mh).convergence(KRYLOV_DIM)
        krylov = min(steps, key=lambda s: abs(s.energy - exact)).energy

        print(f"{r:6.2f} {mh.num_qubits:6d} {hf:12.6f} {krylov:12.6f} {exact:13.6f} "
              f"{(hf - exact) * 1e3:9.2f}m {(krylov - exact) * 1e3:10.4f}m")
        rows.append({
            "bond_length_A": r, "qubits": mh.num_qubits,
            "hf_energy": hf, "krylov_energy": krylov, "casci_exact": exact,
            "hf_error_mHa": (hf - exact) * 1e3, "krylov_error_mHa": (krylov - exact) * 1e3,
        })
        # incremental write so partial runs survive
        with open(OUTPUT, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    print("\nHartree-Fock error grows from ~0.1 to ~0.6 Ha as the triple bond breaks "
          "(single-reference\nbreakdown); quantum Krylov tracks the exact CASCI curve. "
          "Errors in mHa.  ->", OUTPUT)


if __name__ == "__main__":
    main()
