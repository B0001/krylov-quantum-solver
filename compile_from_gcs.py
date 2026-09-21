#!/usr/bin/env python3
"""
CERTCHEM DATASET COMPILER: compile_from_gcs.py

Downloads completed result blobs from GCS and writes one CSV per molecular
family, reading each point's identity from the provenance the worker recorded
at compute time.

An earlier version reconstructed the geometry by computing local CASCI
energies for every candidate coordinate and taking the argmin of
|E_solver - E_CASCI|, asserting in a comment that the nearest match was
"mathematically guaranteed to be correct". It was circular: the coordinate
label came out of the reference, so the resulting dataset agreed with CASCI to
<1e-4 Ha by construction and could never falsify anything. A reference that
cannot disagree is not a reference.

Provenance is now stated by the producer. Blobs written before the worker
stamped provenance cannot be labelled honestly and are reported as skipped,
not guessed.

Usage:
    python compile_from_gcs.py [--bucket NAME] [--verify]

    --verify  recompute CASCI from each recorded geometry and check that the
              certified bracket actually contains it. Slow, and the only part
              of this script entitled to call itself a check.
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict

from google.cloud import storage

# Thread pinning must precede any BLAS-backed import (pyscf, under --verify).
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

DEFAULT_BUCKET = "tlsmbl-compute-certchem-results"

BASE_FIELDS = [
    "family",
    "coordinate_var",
    "molecule",
    "basis",
    "active_space",
    "best_estimate_hartree",
    "lower_bound_hartree",
    "upper_bound_hartree",
    "bracket_width_hartree",
    "wall_time_s",
    "job_id",
]
VERIFY_FIELDS = [
    "casci_reference_hartree",
    "error_vs_casci_hartree",
    "bracket_contains_reference",
]


def row_from_blob(data: dict) -> dict | None:
    """Build a dataset row from a result blob, or None if it has no provenance."""
    prov = data.get("provenance")
    if not prov or not prov.get("molecule"):
        return None

    meta = prov.get("metadata") or {}
    return {
        "family": meta.get("family", "unlabelled"),
        "coordinate_var": meta.get("coordinate"),
        "molecule": prov["molecule"],
        "basis": prov.get("basis"),
        "active_space": ",".join(str(n) for n in prov.get("active_space", [])),
        "best_estimate_hartree": data.get("best_estimate_hartree"),
        "lower_bound_hartree": data.get("lower_bound_hartree"),
        "upper_bound_hartree": data.get("upper_bound_hartree"),
        "bracket_width_hartree": data.get("bracket_width_hartree"),
        "wall_time_s": data.get("wall_time_s"),
        "job_id": prov.get("job_id"),
    }


def casci_energy(molecule: str, basis: str, active_space: str) -> float | None:
    """Exact CASCI energy for a recorded geometry -- an independent reference.

    Independent because the geometry comes from the blob's provenance, not
    from this calculation. It is free to disagree with the solver, which is
    the entire point of computing it.
    """
    from pyscf import gto, mcscf, scf

    try:
        nelec, norb = (int(n) for n in active_space.split(","))
        mol = gto.M(atom=molecule, basis=basis, spin=0, verbose=0)
        mc = mcscf.CASCI(scf.RHF(mol).run(), norb, nelec)
        mc.kernel()
        return float(mc.e_tot)
    except Exception as err:  # noqa: BLE001 -- one bad point must not kill the batch
        print(f"     [WARN] CASCI reference failed for {molecule[:40]}...: {err}")
        return None


def verify(row: dict) -> dict:
    """Attach the CASCI reference and whether the certified bracket contains it."""
    ref = casci_energy(row["molecule"], row["basis"], row["active_space"])
    row["casci_reference_hartree"] = ref
    if ref is None:
        row["error_vs_casci_hartree"] = None
        row["bracket_contains_reference"] = None
        return row

    est = row.get("best_estimate_hartree")
    row["error_vs_casci_hartree"] = None if est is None else est - ref

    lo, hi = row.get("lower_bound_hartree"), row.get("upper_bound_hartree")
    if lo is None or hi is None:
        row["bracket_contains_reference"] = None
    else:
        row["bracket_contains_reference"] = bool(lo <= ref <= hi)
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bucket", default=DEFAULT_BUCKET)
    ap.add_argument("--outdir", default="data")
    ap.add_argument(
        "--verify",
        action="store_true",
        help="recompute CASCI per point and check the certified bracket contains it",
    )
    args = ap.parse_args()

    print("=" * 80)
    print("CERTCHEM DATASET COMPILATION (provenance-based)")
    print("=" * 80)
    print(f"Bucket: gs://{args.bucket}")

    try:
        bucket = storage.Client().bucket(args.bucket)
        blobs = [b for b in bucket.list_blobs(prefix="results/") if b.name.endswith(".json")]
    except Exception as err:
        print(f"[FATAL] Could not read the results bucket: {err}")
        print("Authenticate with: gcloud auth application-default login")
        return 1

    print(f"Found {len(blobs)} result blobs.\n")
    if not blobs:
        print("[WARNING] Nothing to compile.")
        return 0

    families = defaultdict(list)
    skipped_no_provenance = 0
    unreadable = 0

    for blob in blobs:
        try:
            data = json.loads(blob.download_as_text())
        except Exception as err:
            print(f"  [WARN] Unreadable blob {blob.name}: {err}")
            unreadable += 1
            continue

        row = row_from_blob(data)
        if row is None:
            skipped_no_provenance += 1
            continue
        families[row["family"]].append(row)

    labelled = sum(len(v) for v in families.values())
    print(f"Labelled from provenance: {labelled}")
    if skipped_no_provenance:
        print(
            f"Skipped (no provenance -- predates the stamping worker): "
            f"{skipped_no_provenance}. These cannot be labelled honestly; "
            f"recompute them rather than inferring their geometry."
        )
    if unreadable:
        print(f"Unreadable blobs: {unreadable}")

    if args.verify:
        print(f"\n[VERIFY] Recomputing CASCI references for {labelled} points...")
        for rows in families.values():
            for row in rows:
                verify(row)

    print("\n" + "=" * 80)
    os.makedirs(args.outdir, exist_ok=True)
    fields = BASE_FIELDS + (VERIFY_FIELDS if args.verify else [])
    violations = 0

    for family, rows in sorted(families.items()):
        rows.sort(key=lambda r: (r["coordinate_var"] is None, r["coordinate_var"]))
        path = os.path.join(args.outdir, f"cloud_sweep_{family}.csv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

        note = ""
        if args.verify:
            bad = [r for r in rows if r["bracket_contains_reference"] is False]
            violations += len(bad)
            note = f"  [{len(bad)} bracket violations]" if bad else "  [brackets hold]"
        print(f"  -> {path}  ({len(rows)} points){note}")

    if args.verify and violations:
        print(
            f"\n[FAIL] {violations} certified brackets do NOT contain the CASCI "
            f"reference. A certified bracket that excludes the true energy is a "
            f"broken certificate, not a small error."
        )
        return 1

    print(f"\nCompiled {len(families)} families.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
