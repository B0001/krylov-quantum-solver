#!/usr/bin/env python3
"""
Per-sector validate-and-cost sweep -- each spin sector self-reports to telemetry.

This wires validate_and_cost into the spin-sector machinery of run_nbn_sqd_sweep.py: for every
physically reachable spin sector of a structure it tapers, cross-checks the near-term energy
across the available solvers, prices the FT run, and writes one flat CSV row per sector. Run it
alongside the SQD sweep so each Nb3 sector carries a validated energy and an FT cost tag in your
telemetry, not just a raw SQD number.

  from_cif("data/nb_structures/NbN_mp-2634.cif", cas_electrons=8, cas_orbitals=8,
           output_csv="nb_sector_verdicts.csv")

Each sector is guarded, so a non-converging SCF or out-of-regime stage logs a status instead of
killing the sweep. Reuses run_nbn_sqd_sweep.py and validate_and_cost.py (run in the chem-ft env
for the full three-stage verdict).
"""

import os
import sys
import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_nbn_sqd_sweep import (
    load_geometry, valid_spin_sectors, multiplicity_name, build_scf, integrals_for_spin,
)
from validate_and_cost import validate_and_cost


def _flatten(mol_spin, name, nelec, norb, report):
    tap, xc, ft = report["taper"], report["cross_check"], report["ft_cost"]
    row = {
        "mol_spin": mol_spin, "multiplicity": name,
        "n_alpha": nelec[0], "n_beta": nelec[1], "norb": norb,
        "reference_energy": xc.get("reference"),
        "cross_check_agree": xc.get("agree"),
        "max_dev_mHa": xc.get("max_dev_mHa"),
        "methods_run": "+".join(xc.get("results", {}).keys()),
        "methods_skipped": "+".join(xc.get("skipped", [])),
        "n_qubits_original": tap.get("n_qubits_original"),
        "n_qubits_tapered": tap.get("n_qubits_tapered"),
        "ft_threshold": None, "ft_lambda_DF": None, "ft_ccsd_t_err_mHa": None,
        "ft_toffoli": None, "ft_logical_qubits": None, "ft_status": "n/a",
    }
    if ft is None:
        row["ft_status"] = "no_openfermion"
    elif isinstance(ft, dict) and ft.get("recommended"):
        rec = ft["recommended"]
        row.update(ft_threshold=rec["thresh"], ft_lambda_DF=rec["lambda_DF"],
                   ft_ccsd_t_err_mHa=rec["ccsd_t_err_mHa"], ft_toffoli=rec["toffoli_total"],
                   ft_logical_qubits=rec["logical_qubits"], ft_status="costed")
    elif isinstance(ft, dict):
        row["ft_status"] = "out_of_regime" if "error" in ft else "no_threshold_met"
    return row


def _to_frame(rows):
    """Build a polars frame from possibly-ragged rows (a guarded sector may emit only a
    status row). Polars needs a uniform schema, so union the keys and null-fill the gaps --
    matching pandas' default ragged-dict behaviour while preserving first-seen column order."""
    columns = list(dict.fromkeys(k for row in rows for k in row))
    return pl.DataFrame([{c: row.get(c) for c in columns} for row in rows],
                        infer_schema_length=None)


def validate_sweep(atom_str, basis_set, ecp_dict, cas_electrons, cas_orbitals,
                   output_csv="sector_verdicts.csv"):
    """Run validate_and_cost over all physical spin sectors and write the verdict CSV."""
    sectors = valid_spin_sectors(cas_electrons, cas_orbitals)
    print(f"[VERDICT SWEEP] spin sectors (mol.spin = 2S): {sectors}")
    rows = []
    for mol_spin in sectors:
        name = multiplicity_name(mol_spin)
        try:
            _, mf = build_scf(atom_str, basis_set, ecp_dict, mol_spin)
            h1, eri, e_core, nelec, _ = integrals_for_spin(mf, cas_orbitals, cas_electrons, mol_spin)
            report = validate_and_cost(h1, eri, e_core, nelec, h1.shape[0])
            row = _flatten(mol_spin, name, nelec, h1.shape[0], report)
            row["status"] = "OK"
            print(f"  [{name:>8}] (na,nb)=({nelec[0]},{nelec[1]})  E={row['reference_energy']:.6f}  "
                  f"agree={row['cross_check_agree']}  ft={row['ft_status']}")
        except Exception as exc:  # noqa: BLE001
            row = {"mol_spin": mol_spin, "multiplicity": name, "status": f"FAILED: {exc}"}
            print(f"  [{name:>8}] FAILED: {exc}")
        rows.append(row)

    _to_frame(rows).write_csv(output_csv)
    print(f"[COMPLETE] {len(rows)} sector verdicts -> {output_csv}")
    return rows


def from_cif(cif_filepath, cas_electrons=8, cas_orbitals=8, output_csv="sector_verdicts.csv"):
    atom_str, basis_set, ecp_dict = load_geometry(cif_filepath)
    return validate_sweep(atom_str, basis_set, ecp_dict, cas_electrons, cas_orbitals, output_csv)


if __name__ == "__main__":
    # No CIF on hand here: drive the same code path from a directly built molecule so the
    # per-sector verdict and CSV emission are exercised end to end.
    rows = validate_sweep("Li 0 0 0; H 0 0 1.6", "sto-3g", {},
                          cas_electrons=2, cas_orbitals=2, output_csv="sector_verdicts_demo.csv")
    print("\nCSV columns:", _to_frame(rows).columns)
