#!/usr/bin/env python3
"""
Honest benchmark: real-time quantum Krylov vs exact FCI on small molecules.

Replaces the retracted README "Empirical Benchmarks" table (which reported machine-epsilon
collapse as a feature and energies that swung from +0.003 to -48.7 Ha). Every number here is
the solver's error against an exact diagonalisation reference in the same energy frame.

Run:  python benchmark_krylov.py   ->  prints a table and writes data/krylov_benchmark.csv
"""
import csv

from hybrid_quantum_solver.pipeline import run_geometry

SYSTEMS = {
    "H2 (0.74 A, equilibrium)": dict(atom="H 0 0 0; H 0 0 0.74"),
    "H2 (2.0 A, stretched)":    dict(atom="H 0 0 0; H 0 0 2.0"),
    "H4 chain (1.0 A)":         dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
    "LiH (1.6 A)":              dict(atom="Li 0 0 0; H 0 0 1.6"),
}
KRYLOV_DIM = 10


def _pyscf_fci(atom, basis="sto3g"):
    """Independent, fast FCI reference (no dense qubit-matrix diagonalisation)."""
    from pyscf import gto, scf, fci
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    e, _ = fci.FCI(mf).kernel()
    return float(e)


def main():
    header = f"{'system':26s} {'qubits':>6} {'M':>3} {'HF err':>11} {'Krylov err':>12} {'>= FCI':>7}"
    print(header)
    print("-" * len(header))
    rows = []
    for name, spec in SYSTEMS.items():
        fci = _pyscf_fci(spec["atom"])
        r = run_geometry(**spec, krylov_dim=KRYLOV_DIM, reference="none")
        hf_err_mha = (r.hf_energy - fci) * 1e3
        kr_err_mha = (r.computed_energy - fci) * 1e3
        respects_floor = r.computed_energy >= fci - 1e-7
        r.reference_energy = fci
        print(f"{name:26s} {r.n_qubits:6d} {r.krylov_dim:3d} "
              f"{hf_err_mha:9.3f}m {kr_err_mha:10.4f}m {str(respects_floor):>7}")
        rows.append({
            "system": name, "qubits": r.n_qubits, "krylov_dim": r.krylov_dim,
            "fci_energy": r.reference_energy, "hf_energy": r.hf_energy,
            "krylov_energy": r.computed_energy,
            "hf_error_mHa": hf_err_mha, "krylov_error_mHa": kr_err_mha,
            "respects_variational_floor": respects_floor,
        })

    with open("data/krylov_benchmark.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print("\nErrors in milli-Hartree (mHa). 1 kcal/mol = 1.594 mHa ('chemical accuracy').")
    print("Wrote data/krylov_benchmark.csv")


if __name__ == "__main__":
    main()
