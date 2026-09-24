#!/usr/bin/env python3
"""
Pre-registered thermodynamic-limit analysis for the localized H_n table (bead chem-oxm,
specs/SPEC_hchain_largen2.md §11). Written and committed BEFORE the large-n numbers existed; the
fit forms, windows and bar definition below must not be changed after seeing the data.

  * 8 fits: e(n) = a + b/n  and  e(n) = a + b/n + c/n^2, each over n_min in {8, 12, 16, 20} and
    n_max = the largest n in the table.
  * Systematic-inclusive bar: the envelope of `a` over the 8 fits.
    headline = midpoint of [min a, max a];  bar = (max a - min a)/2 + max single-fit stderr.
  * Leave-one-out (LOO): the driver's definition -- the linear 1/n fit over ALL n, minus the same
    fit with the largest n dropped.
  * Bulk estimator: bulk_per_site_energy over the two largest n; compared with the all-n linear fit.

Pure numpy + csv: no pyscf, no block2, so the spec gate that reads the vendored table is instant.

Run:  uv run python hchain_tdl_analysis.py [table.csv]   (default: the vendored spec table)
"""
from __future__ import annotations

import csv
import sys

import numpy as np

from hybrid_quantum_solver.dmrg_reference import bulk_per_site_energy

VENDORED_TABLE = "specs/hchain_tdl_localized_table.csv"
MOTTA_E_INF = -0.540493      # Motta et al., PRX 7, 031059 (2017), STO-6G, R = 1.8 bohr, Ha/atom
N_MINS = (8, 12, 16, 20)
ORDERS = (1, 2)              # 1: a + b/n ; 2: a + b/n + c/n^2


def load_table(path=VENDORED_TABLE):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: int(r["n"]))
    return rows


def fit_1_over_n(ns, e_per_atom, order):
    """Least squares e = sum_k c_k n^-k (k=0..order). Returns (a, stderr_a); stderr is nan when
    there are no residual degrees of freedom."""
    x = 1.0 / np.asarray(ns, dtype=float)
    y = np.asarray(e_per_atom, dtype=float)
    X = np.vander(x, order + 1, increasing=True)          # columns 1, 1/n, 1/n^2
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(y) - (order + 1)
    if dof <= 0:
        return float(coef[0]), float("nan")
    s2 = float(np.sum((y - X @ coef) ** 2)) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    return float(coef[0]), float(np.sqrt(max(cov[0, 0], 0.0)))


def analyse(rows):
    ns = [int(r["n"]) for r in rows]
    epa = [float(r["e_per_atom"]) for r in rows]
    tot = [float(r["e_dmrg_extrap"]) for r in rows]
    fits = []
    for order in ORDERS:
        for n_min in N_MINS:
            sel = [(n, e) for n, e in zip(ns, epa) if n >= n_min]
            a, se = fit_1_over_n([s[0] for s in sel], [s[1] for s in sel], order)
            fits.append({"form": "a+b/n" if order == 1 else "a+b/n+c/n^2", "n_min": n_min,
                         "n_max": ns[-1], "npts": len(sel), "a": a, "stderr": se})
    a_vals = [f["a"] for f in fits]
    ses = [f["stderr"] for f in fits if np.isfinite(f["stderr"])]
    lo, hi = min(a_vals), max(a_vals)
    headline = 0.5 * (lo + hi)
    bar = 0.5 * (hi - lo) + (max(ses) if ses else 0.0)
    e_all, _ = fit_1_over_n(ns, epa, 1)
    e_loo, _ = fit_1_over_n(ns[:-1], epa[:-1], 1)
    e_bulk = bulk_per_site_energy(ns, tot)
    gap = headline - MOTTA_E_INF
    return {
        "fits": fits, "envelope": (lo, hi), "max_stderr": max(ses) if ses else 0.0,
        "headline": headline, "bar": bar,
        "e_inf_linear_all": e_all, "e_inf_loo": e_loo, "loo_shift": abs(e_all - e_loo),
        "e_bulk": e_bulk, "bulk_vs_fit": abs(e_bulk - e_all),
        "motta": MOTTA_E_INF, "gap": gap, "inside": abs(gap) <= bar,
        "gap_in_bars": abs(gap) / bar if bar > 0 else float("inf"),
    }


def report(res):
    print(f"{'form':<14}{'n_min':>6}{'n_max':>6}{'pts':>5}{'a (Ha/atom)':>15}{'stderr':>11}")
    for f in res["fits"]:
        print(f"{f['form']:<14}{f['n_min']:>6}{f['n_max']:>6}{f['npts']:>5}{f['a']:>15.6f}"
              f"{f['stderr']:>11.2e}")
    lo, hi = res["envelope"]
    print(f"\nenvelope of a: [{lo:.6f}, {hi:.6f}]  width {1e3 * (hi - lo):.3f} mHa; "
          f"max stderr {1e3 * res['max_stderr']:.3f} mHa")
    print(f"HEADLINE e_inf = {res['headline']:.6f} +/- {res['bar']:.6f} Ha/atom")
    print(f"linear all-n fit {res['e_inf_linear_all']:.6f}; LOO {res['e_inf_loo']:.6f} "
          f"(shift {1e3 * res['loo_shift']:.3f} mHa/atom)")
    print(f"bulk per-site {res['e_bulk']:.6f} (vs linear all-n fit {1e3 * res['bulk_vs_fit']:.3f} "
          f"mHa/atom)")
    print(f"Motta {res['motta']:.6f}: gap {1e3 * res['gap']:+.3f} mHa/atom = "
          f"{res['gap_in_bars']:.2f} bar-widths -> {'INSIDE' if res['inside'] else 'OUTSIDE'}")




# ---- Vendoring (not analysis): merge the driver CSVs into the tracked table the gate reads. ----
VENDOR_FIELDS = ["n", "e_dmrg_extrap", "e_per_atom", "stderr", "dw_dmax", "regime", "bond_dims",
                 "dw_per_D", "e_per_D", "extrap_method", "fci_energy", "threads", "wall_s", "source"]


def vendor(primary="data/hchain_tdl_localized.csv",
           rerun="data/hchain_tdl_localized_rerun.csv", out=VENDORED_TABLE):
    """A rerun row (larger D ladder) replaces the primary row for the same n."""
    by_n = {}
    for path, tag in ((primary, "ladder 100/200/400"), (rerun, "ladder 200/400/800")):
        with open(path) as f:
            for r in csv.DictReader(f):
                r["source"] = tag
                r["dw_dmax"] = r["dw_per_D"].split("/")[-1]
                by_n[int(r["n"])] = r
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=VENDOR_FIELDS, extrasaction="ignore")
        w.writeheader()
        for n in sorted(by_n):
            w.writerow(by_n[n])


if __name__ == "__main__":
    if sys.argv[1:2] == ["--vendor"]:
        vendor()
        sys.argv.pop(1)
    report(analyse(load_table(sys.argv[1] if len(sys.argv) > 1 else VENDORED_TABLE)))
