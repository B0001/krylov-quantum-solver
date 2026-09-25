#!/usr/bin/env python3
"""
DMRG half-filled 1D Hubbard chains vs the exact Lieb-Wu energy (bead chem-tjr,
specs/SPEC_hubbard_bethe.md §10). Composes `hubbard_chain_integrals` with
`dmrg_energy_extrapolated` and appends one row per (route, U, L) to the TRACKED table
specs/hubbard_lieb_wu_table.csv, which `hubbard_tdl_analysis.py` and the spec gate read.

  route "open": open chain (on-site interaction only, cheap); bulk energy by difference of lengths.
  route "ring": closed-shell ring (hubbard_chain_integrals' default boundary phase), sites put in the
                folded order 0, L-1, 1, L-2, ... so every hopping spans <= 2 MPS sites.

Every ladder starts from the uniform-filling MPS (init_occs = 1 per site): block2's default random
start stalls in the Mott insulator (§10.1). Resumable: rows already in the table are skipped.
block2 only -- never import pyscf/qiskit-aer here (CLAUDE.md segfault note).

Run:  nohup uv run python benchmark_hubbard_lieb_wu.py --route open --U 4 --L 20,40,60,80,100 \
          --bond-dims 200,400,800 > logs/hubbard_open_U4.log 2>&1 &
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import time

import numpy as np

from hubbard_tdl_analysis import FIELDS, VENDORED_TABLE
from hybrid_quantum_solver.dmrg_reference import dmrg_energy_extrapolated
from hybrid_quantum_solver.model_hamiltonians import hubbard_chain_integrals


def fold_order(L):
    """0, L-1, 1, L-2, ...: a ring's wrap bond becomes a range-2 hop on the MPS."""
    return np.array([i // 2 if i % 2 == 0 else L - 1 - i // 2 for i in range(L)])


def integrals(route, L, U):
    m = hubbard_chain_integrals(L, U, open_chain=(route == "open"))
    h1, eri = m.h1, m.eri
    if route == "ring":
        p = fold_order(L)
        h1, eri = h1[np.ix_(p, p)], eri[np.ix_(p, p, p, p)]
    return h1, eri, m.nelec, m.e_core


def done_keys(path):
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        return {(r["route"], float(r["U"]), int(r["L"])) for r in csv.DictReader(f)}


def append_row(path, row):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerow(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--route", choices=("open", "ring"), required=True)
    ap.add_argument("--U", type=float, required=True)
    ap.add_argument("--L", required=True, help="comma-separated even lengths")
    ap.add_argument("--bond-dims", required=True)
    ap.add_argument("--sweeps", type=int, default=8)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--stack-mem-gb", type=float, default=None)
    ap.add_argument("--table", default=VENDORED_TABLE)
    a = ap.parse_args()
    Ds = tuple(int(x) for x in a.bond_dims.split(","))
    for L in (int(x) for x in a.L.split(",")):
        if (a.route, a.U, L) in done_keys(a.table):
            print(f"skip {a.route} U={a.U} L={L} (in table)", flush=True)
            continue
        h1, eri, ne, ec = integrals(a.route, L, a.U)
        scratch = f"./.dmrg_tmp/hub_{a.route}_{a.U}_{L}_{os.getpid()}"
        t0 = time.time()
        res = dmrg_energy_extrapolated(
            h1, eri, ne, ec, bond_dims=Ds, n_sweeps_per=a.sweeps, protocol="perD",
            scratch=scratch, n_threads=a.threads, init_occs=np.ones(L),
            stack_mem=None if a.stack_mem_gb is None else int(a.stack_mem_gb * 2**30),
        )
        wall = time.time() - t0
        shutil.rmtree(scratch, ignore_errors=True)
        row = {
            "route": a.route, "U": a.U, "L": L,
            "bond_dims": "/".join(str(p[0]) for p in res.per_D),
            "e_per_D": "/".join(repr(p[2]) for p in res.per_D),
            "dw_per_D": "/".join(f"{p[1]:.3e}" for p in res.per_D),
            "stage_dE": "/".join(f"{x:.3e}" for x in res.stage_dE),
            "e_extrap": repr(res.energy), "stderr": f"{res.stderr:.3e}",
            "method": res.method, "regime": res.regime,
            "threads": a.threads, "wall_s": f"{wall:.0f}",
        }
        append_row(a.table, row)
        print(f"{a.route} U={a.U} L={L}  E={res.energy:.10f}  e/L={res.energy / L:.8f}  "
              f"regime={res.regime}  E(Dmax)-E={res.per_D[-1][2] - res.energy:.2e}  "
              f"stage_dE={max(res.stage_dE):.1e}  wall={wall:.0f}s", flush=True)


if __name__ == "__main__":
    main()
