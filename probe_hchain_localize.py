#!/usr/bin/env python3
"""
Cost probe for bead chem-rhw (specs/SPEC_hchain_largen2.md §10): H_n DMRG in canonical RHF orbitals
vs Loewdin site orbitals, IDENTICAL sweep schedule, threads and memory pool.

Each run is its own subprocess, so it gets its own wall-clock cap, its own peak RSS (wait4
rusage), and its own block2 runtime. block2's per-sweep log (iprint=1) streams to
data/probe_<tag>.log; when a run hits the cap, the last completed sweep in that log is what it
reached.

Run:  uv run python probe_hchain_localize.py [--n 20] [--D 400] [--cap-min 90]
      -> data/hchain_localize_probe.json (untracked; the numbers are transcribed into the spec)
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time

# Same schedule for both bases: D held fixed; noise annealed to zero; block2's default early stop
# (|dE| < 1e-8 after a noise-free sweep) applies identically to both.
N_SWEEPS = 20
NOISES = [1e-4] * 4 + [1e-5] * 4 + [1e-6] * 4 + [0.0] * 8
THRDS = [1e-10] * N_SWEEPS
SEED = 1234


def worker(n, D, localize, threads, stack_gb, scratch):
    from benchmark_hchain_tdl import integrals
    h1, eri, ne, ec, e_hf = integrals(n, localize=localize)
    import block2
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    block2.Random.rand_seed(SEED)
    drv = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=threads,
                     stack_mem=int(stack_gb * 1024**3))
    drv.initialize_system(n_sites=n, n_elec=sum(ne), spin=0, orb_sym=None)
    t0 = time.time()
    mpo = drv.get_qc_mpo(h1e=h1, g2e=eri, ecore=ec, iprint=0)
    ket = drv.get_random_mps(tag="KET", bond_dim=D, nroots=1)
    e = drv.dmrg(mpo, ket, n_sweeps=N_SWEEPS, bond_dims=[D] * N_SWEEPS, noises=NOISES,
                 thrds=THRDS, iprint=1)
    dws = list(drv._dmrg.discarded_weights)
    out = {"energy": float(e), "discarded_weight": float(dws[-1]), "sweeps": len(dws),
           "e_hf": float(e_hf), "dmrg_seconds": time.time() - t0}
    print("RESULT " + json.dumps(out), flush=True)


_SWEEP_E = re.compile(r"E\s*=\s*(-?\d+\.\d+)")
_SWEEP_DW = re.compile(r"DW\s*=\s*([0-9.eE+-]+)")


def last_sweep(log_path):
    """(energy, discarded weight, completed sweeps) from block2's iprint=1 per-sweep lines."""
    e = dw = None
    k = 0
    with open(log_path) as f:
        for line in f:
            if line.startswith("Time elapsed"):
                m, w = _SWEEP_E.search(line), _SWEEP_DW.search(line)
                if m:
                    e, k = float(m.group(1)), k + 1
                if w:
                    dw = float(w.group(1))
    return e, dw, k


def run_one(tag, n, D, localize, threads, stack_gb, cap_s):
    os.makedirs("data", exist_ok=True)
    log = f"data/probe_{tag}.log"
    cmd = [sys.executable, "-u", __file__, "--worker", "--n", str(n), "--D", str(D),
           "--threads", str(threads), "--stack-mem-gb", str(stack_gb),
           "--scratch", f".dmrg_tmp/probe_{tag}"] + (["--localize"] if localize else [])
    t0 = time.time()
    with open(log, "w") as fh:
        p = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT)
        status, ru = None, None
        while True:
            pid, st, ru = os.wait4(p.pid, os.WNOHANG)
            if pid:
                status = os.waitstatus_to_exitcode(st)
                break
            if time.time() - t0 > cap_s:
                p.kill()
                _, st, ru = os.wait4(p.pid, 0)
                status = "capped"
                break
            time.sleep(2)
    wall = time.time() - t0
    rec = {"tag": tag, "n": n, "D": D, "localize": localize, "threads": threads,
           "stack_mem_gb": stack_gb, "wall_seconds": wall, "peak_rss_gb": ru.ru_maxrss / 1024**2,
           "status": status}
    result = [ln for ln in open(log) if ln.startswith("RESULT ")]
    if result:
        rec.update(json.loads(result[-1][7:]))
    else:
        rec["energy"], rec["discarded_weight"], rec["sweeps"] = last_sweep(log)
    print(json.dumps(rec), flush=True)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--D", type=int, default=400)
    ap.add_argument("--localize", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--stack-mem-gb", type=float, default=6.0)
    ap.add_argument("--scratch", default=".dmrg_tmp/probe")
    ap.add_argument("--cap-min", type=float, default=90.0)
    ap.add_argument("--runs", default="loc,can",
                    help="comma list of loc|can, optionally with @D, e.g. loc@100,loc,can")
    ap.add_argument("--out", default="data/hchain_localize_probe.json")
    a = ap.parse_args()
    if a.worker:
        worker(a.n, a.D, a.localize, a.threads, a.stack_mem_gb, a.scratch)
        return
    recs = json.load(open(a.out)) if os.path.exists(a.out) else []
    for spec in a.runs.split(","):
        basis, _, d = spec.partition("@")
        D = int(d) if d else a.D
        tag = f"n{a.n}_{basis}_D{D}"
        rec = run_one(tag, a.n, D, basis == "loc", a.threads, a.stack_mem_gb, a.cap_min * 60)
        recs = [r for r in recs if r["tag"] != tag] + [rec]
        with open(a.out, "w") as f:
            json.dump(recs, f, indent=1)


if __name__ == "__main__":
    main()
