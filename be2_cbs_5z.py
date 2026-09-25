#!/usr/bin/env python3
"""
be2_cbs_5z.py -- bead chem-mom: does cc-pV5Z move the Be2 CBS well depth?

SPEC_be2_cbs.md attributes the remaining ~461 cm^-1 underbinding (De=469 vs 929.7 cm^-1) to the
METHOD (CASCI(4,8)+NEVPT2), not the basis -- but its "CBS" is only a TZ/QZ two-point Helgaker
extrapolation. This driver runs casci_nevpt2_point on be2_cbs.py's 13-point grid at TZ, QZ and 5Z,
then compares the TZ/QZ and QZ/5Z CBS wells with the SAME wide quartic fit be2_cbs.__main__ uses.

Pre-registered criterion (bead chem-mom, fixed before any 5Z number was seen):
  |De_CBS(QZ/5Z) - De_CBS(TZ/QZ)| <~ 20 cm^-1  -> method attribution CONFIRMED;
  a material move toward 929.7                 -> attribution WRONG, correct the spec.

Resumable: every point is appended to a TRACKED CSV (results/be2_cbs_5z/points.csv, data/ is
gitignored) as soon as it finishes; a restart skips points already on disk.
"""
from __future__ import annotations

import csv
import os
import platform
import time

import numpy as np

from be2_cbs import HA2CM, casci_nevpt2_point, cbs_extrapolate_correlation

RS = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 4.0, 6.0, 8.0]  # be2_cbs.__main__ grid
BASES = ["ccpvtz", "ccpvqz", "ccpv5z"]
CARDINAL = {"ccpvtz": 3, "ccpvqz": 4, "ccpv5z": 5}
OUT_DIR = "results/be2_cbs_5z"
POINTS_CSV = os.path.join(OUT_DIR, "points.csv")
FIELDS = ["R", "basis", "e_casci", "e_corr", "e_tot", "wall_s"]
EXPERIMENT_DE = 929.7


def load_points():
    pts = {}
    if os.path.exists(POINTS_CSV):
        with open(POINTS_CSV) as f:
            for row in csv.DictReader(f):
                pts[(row["basis"], float(row["R"]))] = {k: float(row[k]) for k in
                                                         ("e_casci", "e_corr", "e_tot", "wall_s")}
    return pts


def compute_missing():
    os.makedirs(OUT_DIR, exist_ok=True)
    pts = load_points()
    new_file = not os.path.exists(POINTS_CSV)
    with open(POINTS_CSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()
        for basis in BASES:
            for R in RS:
                if (basis, R) in pts:
                    continue
                t0 = time.time()
                p = casci_nevpt2_point(R, basis)
                wall = time.time() - t0
                row = dict(R=R, basis=basis, e_casci=p.e_casci, e_corr=p.e_corr, e_tot=p.e_tot,
                           wall_s=round(wall, 1))
                w.writerow(row)
                f.flush()
                pts[(basis, R)] = row
                print(f"{basis} R={R:4.2f} E_tot={p.e_tot:.8f} E_corr={p.e_corr:.8f} "
                      f"wall={wall:.1f}s", flush=True)
    return load_points()


def quartic_well(curve):
    """be2_cbs.__main__'s wide quartic fit over 2.1<=R<=3.0; De vs the R=8.0 asymptote."""
    win = [R for R in RS if 2.1 <= R <= 3.0]
    p = np.poly1d(np.polyfit(win, [curve[R] for R in win], 4))
    roots = [x.real for x in p.deriv().r if abs(x.imag) < 1e-6 and win[0] <= x.real <= win[-1]]
    Re = min(roots, key=lambda x: p(x))
    return float(Re), float((curve[8.0] - p(Re)) * HA2CM)


def cbs_curve(pts, lo, hi):
    return {R: pts[(hi, R)]["e_casci"] + cbs_extrapolate_correlation(
        CARDINAL[lo], pts[(lo, R)]["e_corr"], CARDINAL[hi], pts[(hi, R)]["e_corr"]) for R in RS}


def summarize(pts):
    curves = {b: {R: pts[(b, R)]["e_tot"] for R in RS} for b in BASES}
    curves["CBS(TZ/QZ)"] = cbs_curve(pts, "ccpvtz", "ccpvqz")
    curves["CBS(QZ/5Z)"] = cbs_curve(pts, "ccpvqz", "ccpv5z")
    lines = ["label,Re_A,De_cm-1,De_minus_expt_cm-1"]
    wells = {}
    for label, c in curves.items():
        Re, De = quartic_well(c)
        wells[label] = (Re, De)
        lines.append(f"{label},{Re:.4f},{De:.1f},{De - EXPERIMENT_DE:.1f}")
    shift = wells["CBS(QZ/5Z)"][1] - wells["CBS(TZ/QZ)"][1]
    verdict = "CONFIRMED" if abs(shift) <= 20.0 else "NOT CONFIRMED"
    lines.append(f"# De shift CBS(QZ/5Z) - CBS(TZ/QZ) = {shift:+.1f} cm^-1 -> method attribution "
                 f"{verdict} (pre-registered threshold ~20 cm^-1)")
    walls = {b: np.mean([pts[(b, R)]["wall_s"] for R in RS]) for b in BASES}
    lines.append("# mean wall s/point: " + ", ".join(f"{b}={w:.1f}" for b, w in walls.items())
                 + f"; host {platform.machine()} {os.cpu_count()} cpu")
    text = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "wells.csv"), "w") as f:
        f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    summarize(compute_missing())
