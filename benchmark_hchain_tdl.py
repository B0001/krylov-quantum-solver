#!/usr/bin/env python3
"""
H_n thermodynamic-limit study (specs/SPEC_hchain_tdl.md).

Per-atom ground-state energy of the minimal-basis hydrogen chain, with TWO controlled
extrapolations and honest error bars:

  1. bond dimension D -> infinity, per n, via the discarded-weight rule
     (hybrid_quantum_solver.dmrg_reference.dmrg_energy_extrapolated), and
  2. system size n -> infinity, via E(n)/n ~ e_inf + a/n  (thermodynamic_limit_fit).

Validated against exact FCI where tractable. Resumable: completed n are skipped on rerun.

HONEST SCOPE (see SPEC §2, §7): minimal-basis MODEL at a fixed geometry, open boundaries; this
reproduces benchmark physics (cf. Motta et al., PRX 7, 031059, 2017), it does not extend it, and
the quantum/statevector solver plays no role at these n.

Run:  python benchmark_hchain_tdl.py [--threads N]   ->  data/hchain_tdl.csv + e_inf estimate
"""
import argparse
import csv
import math
import os

from pyscf import gto, scf, ao2mo

from hybrid_quantum_solver.dmrg_reference import (
    dmrg_energy_extrapolated,
    thermodynamic_limit_fit,
    fci_energy,
    dmrg_available,
)

R_ANG = 1.8 * 0.529177210903          # 1.8 Bohr ~ 0.9525 A (near the cohesive minimum)
# Fast, laptop-tractable default (~minutes) giving a stable e_inf. The per-D protocol in
# dmrg_energy_extrapolated runs a SEPARATE converged DMRG per bond dimension (clean truncation
# points but ~3x the cost of a single ramp), so large n is slow: n=20 takes many minutes, n=30
# ~an hour. Add larger n on a faster/GPU node (the driver is resumable). A cheaper variant could
# read block2's per-stage get_dmrg_results() from one ramping run -- a documented follow-up.
CHAIN_LENGTHS = [8, 10, 12, 16]
BOND_DIMS = (200, 400, 800)           # D -> inf schedule (near-exact for 1D chains; raise on big nodes)
FCI_DET_CUTOFF = 5_000_000
OUTPUT = "data/hchain_tdl.csv"
FIELDS = ["n", "qubits", "ndet", "hf_energy", "e_dmrg_extrap", "stderr",
          "e_per_atom", "fci_energy", "dmrg_vs_fci", "extrap_method"]


def integrals(n):
    atom = "; ".join(f"H 0 0 {i * R_ANG:.6f}" for i in range(n))
    mol = gto.M(atom=atom, basis="sto-6g", verbose=0)
    mf = scf.RHF(mol).run()
    mo = mf.mo_coeff
    h1 = mo.T @ mf.get_hcore() @ mo
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mo), n)
    ne = (mol.nelectron // 2, mol.nelectron // 2)
    return h1, eri, ne, float(mol.energy_nuc()), mf.e_tot


def load_done(path):
    done = {}
    if os.path.exists(path):
        with open(path) as f:
            for row in csv.DictReader(f):
                done[int(row["n"])] = (float(row["e_per_atom"]), float(row["e_dmrg_extrap"]))
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    if not dmrg_available():
        print("[FATAL] block2 required: pip install block2.")
        return

    print(f"H_n thermodynamic-limit study | R={R_ANG:.4f} A | sto-6g | D-schedule {BOND_DIMS}")
    hdr = (f"{'n':>4} {'qubits':>6} {'ndet':>16} {'HF':>12} {'E_extrap':>13} {'stderr':>9} "
           f"{'E/atom':>11} {'FCI':>13} {'|D-FCI|':>9}")
    print(hdr)
    print("-" * len(hdr))

    os.makedirs("data", exist_ok=True)
    done = load_done(OUTPUT)
    file_exists = os.path.exists(OUTPUT)
    per_atom = dict((n, v[0]) for n, v in done.items())

    for n in CHAIN_LENGTHS:
        if n in done:
            print(f"{n:>4}  (cached)")
            continue
        h1, eri, ne, ec, e_hf = integrals(n)
        ndet = math.comb(n, ne[0]) * math.comb(n, ne[1])
        res = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=BOND_DIMS, n_threads=args.threads)
        e_fci = fci_energy(h1, eri, ne, ec) if ndet <= FCI_DET_CUTOFF else None
        d_fci = abs(res.energy - e_fci) if e_fci is not None else None
        epa = res.energy / n
        per_atom[n] = epa

        def c(x, w=13, p=6): return f"{x:{w}.{p}f}" if x is not None else f"{'--':>{w}}"
        cf = f"{d_fci:9.1e}" if d_fci is not None else f"{'--':>9}"
        print(f"{n:>4} {2 * n:6d} {ndet:16,d} {c(e_hf, 12)} {c(res.energy)} {res.stderr:9.1e} "
              f"{epa:11.6f} {c(e_fci)} {cf}")

        row = {"n": n, "qubits": 2 * n, "ndet": ndet, "hf_energy": e_hf,
               "e_dmrg_extrap": res.energy, "stderr": res.stderr, "e_per_atom": epa,
               "fci_energy": e_fci, "dmrg_vs_fci": d_fci, "extrap_method": res.method}
        with open(OUTPUT, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            if not file_exists:
                w.writeheader()
                file_exists = True
            w.writerow(row)

    ns = sorted(per_atom)
    if len(ns) >= 3:
        e_inf, stderr = thermodynamic_limit_fit(ns, [per_atom[n] for n in ns])
        e_inf_lo, _ = thermodynamic_limit_fit(ns[:-1], [per_atom[n] for n in ns[:-1]])
        print(f"\nThermodynamic limit  e_inf = {e_inf:.6f} +/- {stderr:.6f} Ha/atom")
        print(f"  leave-one-out (drop n={ns[-1]}): {e_inf_lo:.6f}  "
              f"(shift {abs(e_inf - e_inf_lo) * 1e3:.2f} mHa/atom)")
    print(f"\nMinimal-basis model, open chain -- reproduces benchmark physics, not a new result.  ->  {OUTPUT}")


if __name__ == "__main__":
    main()
