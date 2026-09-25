#!/usr/bin/env python3
"""
NbN CAS(14,14): is the committed DMRG reference in the right spin sector? (bead chem-dc7)

SPEC_nbn_dmrg_reference.md records E(10,4) = -110.046028 (high-spin, the gated reference) and a
never-gated low-spin E(7,7) = -110.042500 "3.5 mHa above". This driver re-derives both:

  cif      reconstruct the gitignored data/nb_structures/NbN_mp-2634.cif and regenerate
           data/nbn_scf.chk through benchmark_nbn.ground_state_mf (the committed SCF path)
  fci      EXACT active-space FCI per spin sector with <S^2> reported; with --twos, a spin
           penalty pins S. The unconstrained Ms=0 (7,7) root is the global CAS ground of ANY S
           (spin-restricted CASCI integrals => spin-free H, every S has an Ms=0 component).
  dmrg     one named schedule of nbn_dmrg_reference.SCHEDULES in a chosen sector. block2 runs in
           SU(2) mode with spin = na - nb, so nelec=(7,7) is the S=0 SINGLET, (10,4) is S=3.
           --orbitals no rotates the CAS into the target state's own natural orbitals (from the
           exact FCI 1-RDM): the energy must be invariant; the discarded weight need not be.

Every run appends one JSON line (energy, per-D weights, wall time, host) to
results/nbn_low_spin/runs.jsonl -- data/ is gitignored, the record is not.

pyscf-only for `fci`/`cif`; `dmrg` imports block2 -- one mode per process (CLAUDE.md isolation).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import time

import numpy as np

RESULTS = "results/nbn_low_spin/runs.jsonl"
CIF = "data/nb_structures/NbN_mp-2634.cif"
NO_CACHE = "data/nbn_no_{na}_{nb}.npz"

# mp-2634 is WC-type NbN, P-6m2: Nb (0,0,0), N (1/3,2/3,1/2). The 2-atom "cluster" the pipeline
# builds is therefore a diatomic whose only geometric parameter is d(Nb-N) = sqrt(a^2/3 + c^2/4).
# The original CIF was never committed and materialsproject.org is unreachable from this session;
# MP's own description gives "All Nb-N bond lengths are 2.25 A" (2 decimals). RECONSTRUCTION, not
# the original file: see SPEC_nbn_low_spin.md for the energy fingerprint showing it does NOT
# reproduce the committed -110.046028.
D_NBN = 2.25
A_HEX = 2.95


def write_cif(path: str = CIF, d: float = D_NBN, a: float = A_HEX) -> None:
    c = 2 * math.sqrt(d * d - a * a / 3)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f"""# RECONSTRUCTED WC-type NbN (mp-2634 geometry class), d(Nb-N) = {d} A -- see nbn_low_spin.py
data_NbN
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   {a:.10f}
_cell_length_b   {a:.10f}
_cell_length_c   {c:.10f}
_cell_angle_alpha   90.0
_cell_angle_beta   90.0
_cell_angle_gamma   120.0
_symmetry_Int_Tables_number   1
loop_
 _symmetry_equiv_pos_as_xyz
  'x, y, z'
loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Nb0  Nb  0.0  0.0  0.0  1
  N1  N  0.3333333333333333  0.6666666666666667  0.5  1
""")


def _record(row: dict) -> None:
    row.update(host=platform.node(), cpu=platform.processor() or platform.machine(),
               ncpu=os.cpu_count(), stamp=time.strftime("%Y-%m-%dT%H:%M:%S"))
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "a") as fh:
        fh.write(json.dumps(row) + "\n")
    print("RECORD", json.dumps(row), flush=True)


def spin_fci(h1, eri, nelec, e_core, twos=None, nroots=1, shift=0.01, max_cycle=300):
    """Exact CAS FCI in the Ms=(na-nb)/2 sector; with ``twos`` a spin penalty pins S=twos/2.
    Returns [(E, <S^2>)...] and the CI vectors.

    The penalty is shift*(S^2 - ss)^2. A large shift (0.5 put the septet +72 Ha up) wrecks the
    Davidson preconditioner -- the first singlet attempt ran 300 cycles to a nonsense -108.53 Ha.
    For an even-S target in a balanced sector, direct_spin0 (spin-symmetric CI: even S only)
    removes the odd-S states outright, so a small shift only has to lift S=2,4,6."""
    from pyscf import fci

    norb = h1.shape[0]
    even = twos is not None and twos % 2 == 0 and nelec[0] == nelec[1]
    solver = (fci.direct_spin0 if even else fci.direct_spin1).FCI()
    solver.max_cycle, solver.conv_tol, solver.nroots = max_cycle, 1e-10, nroots
    if twos is not None:
        s = twos / 2
        solver = fci.addons.fix_spin_(solver, shift=shift, ss=s * (s + 1))
    e, ci = solver.kernel(h1, eri, norb, nelec, ecore=e_core)
    if not np.all(solver.converged):
        raise RuntimeError(f"FCI did not converge: {e}")
    e, ci = np.atleast_1d(e), (ci if nroots > 1 else [ci])
    out = [(float(ei), float(fci.spin_op.spin_square0(c, norb, nelec)[0])) for ei, c in zip(e, ci)]
    return out, ci


def natural_orbitals(ci, norb, nelec):
    from pyscf import fci

    dm1 = fci.direct_spin1.make_rdm1(ci, norb, nelec)
    occ, u = np.linalg.eigh(dm1)
    order = np.argsort(-occ)
    return occ[order], u[:, order]


def rotate(h1, eri, u):
    h1n = u.T @ h1 @ u
    erin = np.einsum("pqrs,pi,qj,rk,sl->ijkl", eri, u, u, u, u, optimize=True)
    return h1n, erin


def cmd_cif(_):
    write_cif()
    import benchmark_nbn as bn

    if os.path.exists(bn.CHKFILE):
        raise SystemExit(f"{bn.CHKFILE} exists; remove it to regenerate")
    t = time.time()
    mf = bn.ground_state_mf(CIF)
    _record(dict(kind="scf", d_ang=D_NBN, spin=int(mf.mol.spin), e_scf=float(mf.e_tot),
                 wall_s=round(time.time() - t, 1)))


def cmd_scan(args):
    """Reconstruction-independence check: at each d(Nb-N), regenerate the spin-scanned SCF (own
    chk under .dmrg_tmp/) and compare exact FCI in the (10,4) sector (lowest S>=3) with (9,5)
    (lowest S>=2). If (9,5) is lower at every d, the committed S=3 sector is not the CAS ground
    whatever the lost CIF's exact bond length was."""
    import benchmark_nbn as bn
    from pyscf import ao2mo, mcscf

    for d in args.d:
        tag = f".dmrg_tmp/scan_{d:.4f}"
        os.makedirs(".dmrg_tmp", exist_ok=True)
        write_cif(tag + ".cif", d)
        bn.CHKFILE = tag + ".chk"
        if os.path.exists(bn.CHKFILE):
            os.remove(bn.CHKFILE)
        t = time.time()
        mf = bn.ground_state_mf(tag + ".cif")
        row = dict(kind="scan", d_ang=d, scf_spin=int(mf.mol.spin), e_scf=float(mf.e_tot))
        for nelec in ((10, 4), (9, 5)):
            cas = mcscf.CASCI(mf, 14, nelec)
            h1, e_core = cas.get_h1eff()
            eri = ao2mo.restore(1, cas.get_h2eff(), 14)
            (e, ss), = spin_fci(h1, eri, nelec, e_core)[0]
            row[f"E_{nelec[0]}{nelec[1]}"], row[f"s2_{nelec[0]}{nelec[1]}"] = e, ss
        row["wall_s"] = round(time.time() - t, 1)
        _record(row)


def cmd_fci(args):
    from nbn_dmrg_reference import load_nbn_cas

    nelec = tuple(args.nelec)
    h1, eri, nelec, e_core = load_nbn_cas(nelec=nelec)
    t = time.time()
    roots, ci = spin_fci(h1, eri, nelec, e_core, twos=args.twos, nroots=args.nroots)
    wall = round(time.time() - t, 1)
    for i, (e, ss) in enumerate(roots):
        _record(dict(kind="fci", nelec=list(nelec), twos_pinned=args.twos, root=i, energy=e,
                     s2=ss, ndet=math.comb(14, nelec[0]) * math.comb(14, nelec[1]), wall_s=wall))
    if args.save_no:
        occ, u = natural_orbitals(ci[0], 14, nelec)
        os.makedirs("data", exist_ok=True)
        np.savez(NO_CACHE.format(na=nelec[0], nb=nelec[1]), occ=occ, u=u)
        _record(dict(kind="natocc", nelec=list(nelec), twos_pinned=args.twos,
                     occ=[round(float(x), 6) for x in occ]))


def cmd_dmrg(args):
    from hybrid_quantum_solver.dmrg_reference import dmrg_energy_extrapolated
    from nbn_dmrg_reference import SCHEDULES, load_nbn_cas

    h1, eri, nelec, e_core = load_nbn_cas(nelec=tuple(args.nelec), chk=args.chk)
    if args.orbitals == "no":
        u = np.load(NO_CACHE.format(na=nelec[0], nb=nelec[1]))["u"]
        h1, eri = rotate(h1, eri, u)
    kw = dict(SCHEDULES[args.schedule])
    if args.dims:
        kw["bond_dims"] = tuple(args.dims)
    kw["scratch"] = f"{kw['scratch']}_{nelec[0]}{nelec[1]}_{args.orbitals}_{os.path.basename(args.chk)}"
    os.makedirs(kw["scratch"], exist_ok=True)
    t = time.time()
    r = dmrg_energy_extrapolated(h1, eri, nelec, e_core, n_threads=args.threads, **kw)
    _record(dict(kind="dmrg", chk=args.chk, schedule=args.schedule, nelec=list(nelec), orbitals=args.orbitals,
                 bond_dims=list(kw["bond_dims"]), protocol=kw["protocol"], energy=r.energy,
                 stderr=r.stderr, method=r.method, regime=getattr(r, "regime", None),
                 per_D=[[int(d), float(w), float(e)] for d, w, e in r.per_D],
                 threads=args.threads, wall_s=round(time.time() - t, 1)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("cif").set_defaults(fn=cmd_cif)
    sc = sub.add_parser("scan")
    sc.add_argument("--d", type=float, nargs="+", required=True)
    sc.set_defaults(fn=cmd_scan)
    f = sub.add_parser("fci")
    f.add_argument("--nelec", type=int, nargs=2, default=[7, 7])
    f.add_argument("--twos", type=int, default=None)
    f.add_argument("--nroots", type=int, default=1)
    f.add_argument("--save-no", action="store_true")
    f.set_defaults(fn=cmd_fci)
    d = sub.add_parser("dmrg")
    d.add_argument("--schedule", required=True)
    d.add_argument("--nelec", type=int, nargs=2, default=[7, 7])
    d.add_argument("--orbitals", choices=("scf", "no"), default="scf")
    d.add_argument("--dims", type=int, nargs="*")
    d.add_argument("--threads", type=int, default=2)
    d.add_argument("--chk", default="data/nbn_scf.chk")
    d.set_defaults(fn=cmd_dmrg)
    a = ap.parse_args()
    a.fn(a)
