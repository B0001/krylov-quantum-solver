#!/usr/bin/env python3
"""
chem-5fs: does the floor-guard mechanism finding (BACKLOG "The floor-guard mechanism survives")
hold for the asymmetric-H6 witness and under an irrep filter, not just linear H6?

That BACKLOG entry re-ran SPEC_subspace_floor_resolvability's guard-PASSING-but-invalid witnesses
at conv_tol=1e-13 for linear H6 R=1.0/1.1/1.2 ONLY (no asym-H6, no irrep filter), and concluded the
~1e-4-amplitude reachable level is PHYSICAL (populations in the 1e-6..1e-3 band identical to 6 sig
figs at default and tight SCF). The third recorded escape witness -- asym H6
[0.9,0.9,2.2,0.9,0.9], d=3, M=6 (SPEC_subspace_floor_resolvability.md) -- was never re-checked, nor
was the irrep filter from SPEC_symmetry_reachability.md ever applied to any of them.

This script runs the SAME check (population band stability under conv_tol, plus the irrep-filter
diagnostic) on all four geometries so the asym-H6 and irrep-filter gaps get recorded next to the
linear-H6 numbers.

Run: uv run python scripts/chem_5fs_floor_guard_mechanism_check.py
"""
from __future__ import annotations

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from reachability import (
    TIGHT_SCF_CONV_TOL,
    _dense_hf_projection,
    hf_orbital_irreps,
    scf_symmetry_status,
    symmetry_allowed,
    symmetry_filter_available,
)

BAND_LO, BAND_HI = 1e-6, 1e-3
DEFAULT_CONV_TOL = 1e-9


def _chain(n, r):
    return "; ".join(f"H 0 0 {i * r}" for i in range(n))


def _asym_chain(gaps):
    z = [0.0]
    for g in gaps:
        z.append(z[-1] + g)
    return "; ".join(f"H 0 0 {zi}" for zi in z)


GEOMETRIES = [
    ("linear H6 R=1.0", _chain(6, 1.0)),
    ("linear H6 R=1.1", _chain(6, 1.1)),
    ("linear H6 R=1.2", _chain(6, 1.2)),
    ("asym H6 [0.9,0.9,2.2,0.9,0.9]", _asym_chain([0.9, 0.9, 2.2, 0.9, 0.9])),
]


def band_populations(atom, conv_tol):
    mh = build_molecular_hamiltonian(atom=atom, conv_tol=conv_tol)
    w, vecs, pops = _dense_hf_projection(mh)
    idx = np.flatnonzero((pops > BAND_LO) & (pops < BAND_HI))
    return mh, w, vecs, pops, idx


def main():
    print(f"{'geometry':34s} {'conv_tol':>9s} {'#in-band':>8s} {'max pop':>10s} "
          f"{'SCF':>10s} {'irrep filter':>13s} {'#flagged forbidden':>18s}")
    results = []
    for name, atom in GEOMETRIES:
        broken, de = scf_symmetry_status(atom)
        avail = symmetry_filter_available(atom)
        row_default = band_populations(atom, DEFAULT_CONV_TOL)
        row_tight = band_populations(atom, TIGHT_SCF_CONV_TOL)

        for (conv_tol_label, conv_tol_val), (mh, w, vecs, pops, idx) in zip(
            [("default(1e-9)", DEFAULT_CONV_TOL), ("tight(1e-13)", TIGHT_SCF_CONV_TOL)],
            [row_default, row_tight],
        ):
            n_flag = "n/a"
            if avail:
                allowed = symmetry_allowed(mh, w, vecs)
                n_flag = int(np.sum(~allowed[idx]))
            maxpop = float(pops[idx].max()) if len(idx) else float("nan")
            print(f"{name:34s} {conv_tol_label:>9s} {len(idx):8d} {maxpop:10.3e} "
                  f"{'broken' if broken else 'symmetric':>10s} "
                  f"{'available' if avail else 'unavailable':>13s} {str(n_flag):>18s}")

        mh_d, w_d, vecs_d, pops_d, idx_d = row_default
        mh_t, w_t, vecs_t, pops_t, idx_t = row_tight
        same_count = len(idx_d) == len(idx_t)
        # 6-sig-fig match, same methodology as the linear-H6 BACKLOG entry: compare sorted
        # in-band populations elementwise (only meaningful if the counts match).
        matches_6sf = same_count and bool(
            np.allclose(np.sort(pops_d[idx_d]), np.sort(pops_t[idx_t]), rtol=1e-6, atol=0.0)
        )
        results.append(dict(
            name=name, broken=broken, dE=de, avail=avail,
            n_default=len(idx_d), n_tight=len(idx_t),
            matches_6sf=matches_6sf,
        ))

    print()
    print("Summary (mirrors the linear-H6 BACKLOG claim: 'populations in the 1e-6..1e-3 band "
          "identical to 6 sig figs at default and conv_tol=1e-13'):")
    for r in results:
        verdict = "PHYSICAL (matches linear-H6 behavior)" if r["matches_6sf"] else \
            "DIFFERS from linear-H6 behavior -- investigate"
        print(f"  {r['name']:34s} SCF={'broken' if r['broken'] else 'symmetric':9s} "
              f"dE={r['dE']:+.3e}  irrep-filter={'available' if r['avail'] else 'UNAVAILABLE':12s} "
              f"n_band(default/tight)={r['n_default']}/{r['n_tight']}  -> {verdict}")


if __name__ == "__main__":
    main()
