#!/usr/bin/env python3
"""
Pre-registered thermodynamic-limit analysis for the DMRG Hubbard-chain table (bead chem-tjr,
specs/SPEC_hubbard_bethe.md §10). Written and committed BEFORE the production numbers existed; the
point checks, fit forms and bar definitions below must not be changed after seeing the data.

Per point (one DMRG bond-dimension ladder at one (route, U, L)):
  * OK iff  regime in {converged, truncation}
        and E(D_max) - E_extrap <= spread       (spread = max_D E(D) - min_D E(D), the ladder's own)
        and E_extrap - E(D_max) <= ENERGY_NOISE (an extrapolation above its best point is wrong)
        and every stage's last-sweep |dE| <= STAGE_DE  (stall detector; regime cannot see a stall
                                                        on a 1/D-axis ladder, see §10.1)
  * sigma = max(stderr, |E_extrap - E(D_max)|)   -- the point's extrapolation uncertainty.

Route (a), open chains:  e_a = bulk_per_site_energy over the two largest L (surface term cancels);
  bar_a = |e_a - q(next-lower pair)| + (sigma(L_n) + sigma(L_{n-1})) / (L_n - L_{n-1}),
  q(next-lower pair) = the same difference quotient one step down.
Route (b), closed-shell rings:  e_b = intercept of the least-squares fit e(L) = e_b + a / L^2 over
  every ring L;  bar_b = |e_b - (same fit with the smallest L dropped)| + intercept stderr
  + max_L sigma(L) / L.
Acceptance:  (1) |e_a - lieb_wu(U)| < 1e-3 Ha/site for U in {2,4,8} (route (a) is the headline);
             (2) |e_a - e_b| < bar_a + bar_b.

Pure numpy + csv (+ scipy quadrature for lieb_wu): no pyscf, no block2.

Run:  uv run python hubbard_tdl_analysis.py [table.csv]   (default: the vendored spec table)
"""
from __future__ import annotations

import csv
import sys

import numpy as np

from hybrid_quantum_solver.dmrg_reference import bulk_per_site_energy
from hybrid_quantum_solver.model_hamiltonians import lieb_wu_energy

VENDORED_TABLE = "specs/hubbard_lieb_wu_table.csv"
FIELDS = ["route", "U", "L", "bond_dims", "e_per_D", "dw_per_D", "stage_dE", "e_extrap", "stderr",
          "method", "regime", "threads", "wall_s"]
ENERGY_NOISE = 1e-6      # Ha; same 1 uHa as dmrg_reference.ENERGY_NOISE
STAGE_DE = 1e-6          # Ha; last-sweep energy change allowed for a converged stage
TDL_GATE = 1e-3          # Ha/site, acceptance (1)


def load_table(path=VENDORED_TABLE):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: (r["route"], float(r["U"]), int(r["L"])))
    return rows


def _floats(s):
    return [float(x) for x in s.split("/")]


def point_check(r):
    """Returns (ok, sigma, reasons) for one table row."""
    es = _floats(r["e_per_D"])
    e_x, e_dmax = float(r["e_extrap"]), es[-1]
    spread = max(es) - min(es)
    reasons = []
    if r["regime"] not in ("converged", "truncation"):
        reasons.append(f"regime={r['regime']}")
    if e_dmax - e_x > spread:
        reasons.append(f"undershoot {e_dmax - e_x:.2e} > spread {spread:.2e}")
    if e_x - e_dmax > ENERGY_NOISE:
        reasons.append(f"extrap above E(Dmax) by {e_x - e_dmax:.2e}")
    stage = _floats(r["stage_dE"])
    if not all(np.isfinite(stage)) or max(stage) > STAGE_DE:
        reasons.append(f"stage not converged: max|dE|={max(stage):.2e}")
    sigma = max(float(r["stderr"]), abs(e_x - e_dmax))
    return not reasons, sigma, reasons


def route_a(rows):
    """Open chains -> (e_a, bar_a, detail) from the two largest L; None if < 3 lengths."""
    pts = sorted((int(r["L"]), float(r["e_extrap"]), point_check(r)[1]) for r in rows)
    if len(pts) < 3:
        return None
    Ls, Es, sig = zip(*pts)
    e_a = bulk_per_site_energy(Ls, Es)
    q_lower = bulk_per_site_energy(Ls[:-1], Es[:-1])
    step = Ls[-1] - Ls[-2]
    bar = abs(e_a - q_lower) + (sig[-1] + sig[-2]) / step
    return e_a, bar, {"Ls": list(Ls), "q_lower": q_lower}


def _fit_inv_L2(Ls, eps):
    x = 1.0 / np.asarray(Ls, dtype=float) ** 2
    y = np.asarray(eps, dtype=float)
    X = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(y) - 2
    if dof <= 0:
        return float(coef[0]), 0.0
    s2 = float(np.sum((y - X @ coef) ** 2)) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    return float(coef[0]), float(np.sqrt(max(cov[0, 0], 0.0)))


def route_b(rows):
    """Closed-shell rings -> (e_b, bar_b, detail); None if < 3 lengths."""
    pts = sorted((int(r["L"]), float(r["e_extrap"]), point_check(r)[1]) for r in rows)
    if len(pts) < 3:
        return None
    Ls, Es, sig = zip(*pts)
    eps = [e / L for L, e in zip(Ls, Es)]
    e_b, se = _fit_inv_L2(Ls, eps)
    e_drop, _ = _fit_inv_L2(Ls[1:], eps[1:])
    bar = abs(e_b - e_drop) + se + max(s / L for s, L in zip(sig, Ls))
    return e_b, bar, {"Ls": list(Ls), "drop_smallest": e_drop, "stderr": se}


def analyse(rows):
    out = {"points": [], "U": {}}
    for r in rows:
        ok, sigma, reasons = point_check(r)
        out["points"].append((r["route"], float(r["U"]), int(r["L"]), ok, sigma, reasons))
    for U in sorted({float(r["U"]) for r in rows}):
        sel = [r for r in rows if float(r["U"]) == U]
        a = route_a([r for r in sel if r["route"] == "open"])
        b = route_b([r for r in sel if r["route"] == "ring"])
        lw = lieb_wu_energy(U)
        res = {"lieb_wu": lw, "a": a, "b": b}
        if a:
            res["resid_a"] = a[0] - lw
            res["pass_1"] = abs(a[0] - lw) < TDL_GATE
        if b:
            res["resid_b"] = b[0] - lw
        if a and b:
            res["diff_ab"] = abs(a[0] - b[0])
            res["pass_2"] = res["diff_ab"] < a[1] + b[1]
        out["U"][U] = res
    return out


def report(res):
    print("route  U     L    ok   sigma(Ha)   reasons")
    for route, U, L, ok, sigma, reasons in res["points"]:
        print(f"{route:5s} {U:4.1f} {L:5d}  {'Y' if ok else 'N'}  {sigma:9.2e}   {'; '.join(reasons)}")
    for U, r in res["U"].items():
        print(f"\nU={U}:  lieb_wu = {r['lieb_wu']:.8f}")
        if r["a"]:
            e, bar, d = r["a"]
            print(f"  (a) open  e_a = {e:.8f} +/- {bar:.2e}   Ls={d['Ls']}  q_lower={d['q_lower']:.8f}"
                  f"   resid = {r['resid_a'] * 1e3:+.4f} mHa/site  -> {'PASS' if r['pass_1'] else 'FAIL'} (1)")
        if r["b"]:
            e, bar, d = r["b"]
            print(f"  (b) ring  e_b = {e:.8f} +/- {bar:.2e}   Ls={d['Ls']}"
                  f"   resid = {r['resid_b'] * 1e3:+.4f} mHa/site")
        if "diff_ab" in r:
            print(f"  |e_a - e_b| = {r['diff_ab']:.2e}  vs bar_a+bar_b = {r['a'][1] + r['b'][1]:.2e}"
                  f"  -> {'PASS' if r['pass_2'] else 'FAIL'} (2)")


if __name__ == "__main__":
    report(analyse(load_table(sys.argv[1] if len(sys.argv) > 1 else VENDORED_TABLE)))
