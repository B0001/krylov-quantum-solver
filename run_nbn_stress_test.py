#!/usr/bin/env python3
"""
NbN subspace-scaling study (rewired onto the corrected pipeline).

Sweeps Krylov dimension M and shot budget, comparing the quantum Krylov estimate against
the classical CASCI total energy (the active-space FCI target). Replaces the old 2D
subspace x "noise_variance" sweep that drove the broken orchestrator and reported energies
swinging from +0.003 to -48.7 Ha. See REFACTOR_PLAN.md.

SCIENTIFIC CAVEAT: load_and_compute_integrals builds a finite molecular cluster from the
CIF unit cell (no periodic boundary conditions). This is a code-correctness / convergence
study of the solver, NOT a validated electronic-structure result for solid NbN (Phase 4).
"""
import time

import polars as pl

from hybrid_quantum_solver.chemistry_gateway import load_and_compute_integrals
from hybrid_quantum_solver.pipeline import run_from_integrals

TARGET = "data/nb_structures/NbN_mp-2634.cif"
KRYLOV_DIMS = [4, 8, 12]
SHOTS = [None, 8192, 1024]      # None = exact statevector


def main():
    print("[CLASSICAL] CASCI active-space integrals for NbN ...")
    h1, eri, n_orb, casci_total, e_core, nelecas = load_and_compute_integrals(
        TARGET, cas_electrons=8, cas_orbitals=8
    )
    print(f"  CASCI total = {casci_total:.6f} Ha | qubits = {n_orb * 2} | nelecas = {nelecas}")

    rows = []
    for shots in SHOTS:
        for m in KRYLOV_DIMS:
            t0 = time.time()
            r = run_from_integrals(
                h1, eri, num_particles=nelecas, e_core=float(e_core),
                krylov_dim=m, shots=shots, seed=0,
                reference_energy=float(casci_total), track_convergence=False,
            )
            dt = time.time() - t0
            print(f"  M={m:<3} shots={str(shots):<6}  E={r.computed_energy:.6f}  "
                  f"|err|={r.error_vs_reference:.2e} Ha  ({dt:.1f}s)")
            rows.append({
                "krylov_dim": m, "shots": shots, "qubits": r.n_qubits, "rank": r.rank,
                "computed_energy": r.computed_energy, "casci_total": casci_total,
                "error_vs_casci": r.error_vs_reference, "time_s": round(dt, 2),
            })

    pl.DataFrame(rows).write_csv("data/nbn_krylov_scaling.csv")
    print("\n[DONE] -> data/nbn_krylov_scaling.csv")


if __name__ == "__main__":
    main()
