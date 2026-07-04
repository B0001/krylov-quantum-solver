#!/usr/bin/env python3
"""
NbN CAS(14,14) DMRG reference -- the "DMRG-referenced transition-metal active space" backlog item.

The 2-atom NbN cluster's CAS(14,14) holds ~1.18e7 determinants -- beyond the repo's 5e6 FCI
cutoff -- and its spin-scanned ground state (cached in data/nbn_scf.chk by benchmark_nbn.py)
sits in the HIGH-SPIN nelec=(10,4) (2S=6) sector. Two genuinely independent sweep schedules
(perD: separate converged runs per bond dimension; ramp: one ramping run) referee each other:

    headline (2026-07-04, 16 GB laptop):
      A  perD 400/800/1200:  E = -110.04602843 Ha
      B  ramp 300/600/1200:  E = -110.04602846 Ha      |E_A - E_B| = 3e-8 Ha

HONEST FINDING (gated): this CAS is a SOFT DMRG target. The high-spin sector is
low-entanglement -- D=400 converges to sub-nHa (discarded weight ~1e-9), and even D<=300 sits
within 1 uHa with usable discarded weights. "FCI-intractable by determinant count" did not mean
"strongly correlated"; a hard TM benchmark needs a low-spin/multireference sector or a larger
cluster (out of scope, recorded in specs/SPEC_nbn_dmrg_reference.md). A reference number, not a
materials claim (finite cluster, LANL2DZ ECP, fixed geometry).

pyscf + block2 only -- no qiskit imports (block2 OpenMP isolation; see CLAUDE.md).
"""
from __future__ import annotations

from pyscf import ao2mo, lib, mcscf, scf

from hybrid_quantum_solver.dmrg_reference import dmrg_energy_extrapolated

SCHEDULES = {
    # headline (driver-level, ~20 min + ~3 min on a 16 GB laptop)
    "A": dict(bond_dims=(400, 800, 1200), protocol="perD", n_sweeps_per=10,
              scratch="./.dmrg_tmp/nbn_A", stack_mem=2 * 1024 ** 3, seed=1),
    "B": dict(bond_dims=(300, 600, 1200), protocol="ramp", sweeps_per_stage=6,
              scratch="./.dmrg_tmp/nbn_B", stack_mem=2 * 1024 ** 3, seed=2),
    # cheap CI-gate variants (~2 min total, still in the discarded-weight regime)
    "A'": dict(bond_dims=(100, 200, 300), protocol="perD", n_sweeps_per=8,
               scratch="./.dmrg_tmp/nbn_Ac", stack_mem=1 * 1024 ** 3, seed=11),
    "B'": dict(bond_dims=(80, 160, 300), protocol="ramp", sweeps_per_stage=4,
               scratch="./.dmrg_tmp/nbn_Bc", stack_mem=1 * 1024 ** 3, seed=12),
}


def load_nbn_cas(norb: int = 14, nelec_cas: int = 14, chk: str = "data/nbn_scf.chk"):
    """(h1, eri, nelec, e_core) of the NbN CAS from the cached spin-scanned SCF (no SCF run)."""
    mol = lib.chkfile.load_mol(chk)
    res = lib.chkfile.load(chk, "scf")
    mf = scf.UHF(mol) if mol.nelec[0] != mol.nelec[1] else scf.RHF(mol)
    mf.mo_coeff, mf.mo_occ = res["mo_coeff"], res["mo_occ"]
    mf.mo_energy, mf.e_tot = res["mo_energy"], res["e_tot"]
    cas = mcscf.CASCI(mf, norb, nelec_cas)
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return h1, eri, (int(cas.nelecas[0]), int(cas.nelecas[1])), float(e_core)


def run_schedule(name: str, n_threads: int = 2):
    """One independent DMRG schedule by name ('A'/'B' headline, "A'"/"B'" cheap CI variants)."""
    h1, eri, nelec, e_core = load_nbn_cas()
    kw = dict(SCHEDULES[name])
    return dmrg_energy_extrapolated(h1, eri, nelec, e_core, n_threads=n_threads, **kw)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--schedule", choices=sorted(SCHEDULES), default="A'")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    r = run_schedule(args.schedule, n_threads=args.threads)
    print(f"SCHEDULE {args.schedule}: E={r.energy:.8f} Ha  stderr={r.stderr * 1e3:.4f} mHa  "
          f"method={r.method}")
    print("per_D:", r.per_D)
