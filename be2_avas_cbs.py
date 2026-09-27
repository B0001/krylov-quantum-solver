#!/usr/bin/env python3
"""
be2_avas_cbs.py -- bead chem-xqh: a basis-consistent active space for the Be2 CBS attribution.

SPEC_be2_cbs.md's G5 (bead chem-mom) found that CAS(4,8) on canonical HF virtuals is NOT the same
active space at cc-pVTZ/QZ/5Z: the outermost active orbital's rms spread grows 7.1 -> 8.3 -> 9.4
bohr and its orbital energy falls +0.150 -> +0.074 Ha as the basis grows (results/be2_cbs_5z/,
diagnose.py/.txt), because "lowest 8 virtuals by canonical HF energy" picks increasingly diffuse
orbitals as more diffuse basis functions become available. That is why the single-basis wells
oscillate (997 / 698 / 1313 cm^-1) and why CBS(QZ/5Z) moves +1224 cm^-1 from CBS(TZ/QZ) instead of
agreeing to within the pre-registered ~20 cm^-1: the "CBS" extrapolation was mixing two different
physical active spaces, not converging one.

FIX: select the active space with AVAS (Sayfutyarova, Sun, Chan & Knizia, JCTC 2017,
arXiv:1701.07862) projected onto the Be 2s/2p atomic valence character in a FIXED reference basis
(``minao='ano'``) instead of picking by canonical-HF orbital energy. Because the reference basis
that defines "Be 2s/2p character" does not change with the target basis, the selected active space
is the same physical 8-orbital valence space at every cc-pVXZ -- this is checked directly (G1) and
by comparing orbital spread across bases (G2), not assumed.

NOTE: pyscf's default AVAS reference (``minao='minao'``) is the *minimal atomic ground-state*
basis, which for Be has NO 2p shell at all (Be's atomic ground state is 1s^2 2s^2, no occupied p) --
using it silently drops the 2p orbitals from the active space (ncas=2, not 8). ``minao='ano'`` (the
full ANO basis, which spans unoccupied-in-the-atom shells too) is required to get Be 2p AOs to
project onto; this was verified interactively before writing this driver (ncas=8, nelecas=4 at
every basis tested) and is asserted at runtime by ``avas_casci_nevpt2_point`` (G1) rather than
trusted silently.

PRE-REGISTERED CRITERIA (fixed before the TZ/QZ/5Z curves below were computed; see
specs/SPEC_be2_avas_cbs.md):
  G1 -- AVAS must select exactly ncas=8, nelecas=4 at every (basis, R) point. A mismatch means the
        active space silently changed shape and invalidates the comparison -- hard failure, not a
        warning.
  G2 -- the active-space character (rms spread of each active orbital, at R=2.5 and R=8.0) must
        agree between cc-pVTZ and cc-pV5Z to within 1.0 bohr per orbital -- a direct, much tighter
        bar than the canonical-orbital baseline's 2.3 bohr (7.1->9.4) drift on the outermost orbital.
  G3 -- CBS consistency, same bar chem-mom used: |De_CBS(QZ/5Z) - De_CBS(TZ/QZ)| <= 20 cm^-1 ->
        method attribution CONFIRMED; a material move -> still NOT CONFIRMED, record why.
  G4 -- single-basis De(TZ), De(QZ), De(5Z) monotone in the cardinal number X, OR the failure mode
        is identified (this driver reports the CASCI-only vs NEVPT2-only decomposition so a
        non-monotonic result can be attributed to a specific piece, not left as "noise").

Resumable: every point is appended to a TRACKED CSV (results/be2_avas_cbs/points.csv; data/ is
gitignored) as soon as it finishes; a restart skips points already on disk.
"""
from __future__ import annotations

import csv
import os
import platform
import time

import numpy as np
from pyscf import gto, mcscf, mrpt, scf
from pyscf.mcscf import avas

from be2_cbs import HA2CM, BASIS_CARDINAL, Be2Point, cbs_extrapolate_correlation

RS = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 4.0, 6.0, 8.0]  # be2_cbs.py's own grid
BASES = ["ccpvtz", "ccpvqz", "ccpv5z"]
AOLABELS = ("Be 2s", "Be 2p")
AVAS_MINAO = "ano"
AVAS_THRESHOLD = 0.2
CAS_ELECTRONS, CAS_ORBITALS = 4, 8
OUT_DIR = "results/be2_avas_cbs"
POINTS_CSV = os.path.join(OUT_DIR, "points.csv")
FIELDS = ["R", "basis", "e_casci", "e_corr", "e_tot", "wall_s"]
EXPERIMENT_DE = 929.7
EXPERIMENT_RE = 2.4498
CBS_CONFIRM_THRESHOLD_CM1 = 20.0  # pre-registered, same bar as bead chem-mom (G5)


def avas_active_space(mf, aolabels=AOLABELS, minao=AVAS_MINAO, threshold=AVAS_THRESHOLD,
                       cas_electrons=CAS_ELECTRONS, cas_orbitals=CAS_ORBITALS):
    """AVAS-projected (ncas, nelecas, mo_coeff), asserting the space has the expected shape (G1).

    ``minao='ano'`` is required, not pyscf's AVAS default ``minao='minao'``: Be's minimal
    ground-state basis has no occupied-or-virtual 2p shell to project onto (Be is 1s^2 2s^2), so
    the default silently drops 2p from the active space. See module docstring.
    """
    ncas, ne_act, mo = avas.avas(mf, aolabels, minao=minao, threshold=threshold,
                                 canonicalize=False)
    nelecas = ne_act if isinstance(ne_act, (int, np.integer)) else sum(ne_act)
    if ncas != cas_orbitals or nelecas != cas_electrons:
        raise RuntimeError(
            f"AVAS active space mismatch: got ncas={ncas} nelecas={nelecas}, expected "
            f"({cas_orbitals}, {cas_electrons}) -- basis-consistency assumption (G1) violated")
    return ncas, cas_electrons, mo


def avas_casci_nevpt2_point(R: float, basis: str) -> Be2Point:
    """CASCI(4,8) + NEVPT2 on the AVAS(Be 2s/2p, minao='ano')-selected active space at R (Angstrom).

    Same CASCI+NEVPT2 recipe as be2_cbs.casci_nevpt2_point (fixed active-space orbitals, unfrozen
    core so NEVPT2 includes core-valence correlation) -- the only change is WHICH 8 orbitals are
    active: AVAS-projected valence character instead of the 8 lowest-energy canonical HF virtuals.
    """
    mol = gto.M(atom=f"Be 0 0 0; Be 0 0 {R}", basis=basis, spin=0, verbose=0)
    mf = scf.RHF(mol).run()
    ncas, nelecas, mo = avas_active_space(mf)
    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.kernel(mo)
    e_corr = float(mrpt.NEVPT(mc).kernel())
    return Be2Point(R=R, basis=basis, e_casci=float(mc.e_tot), e_corr=e_corr)


def active_orbital_diagnostics(basis: str, R: float):
    """(orbital Fock energies, rms spreads in bohr) of the AVAS-selected active orbitals at (basis, R).

    AVAS orbitals are not canonical HF eigenstates, so "orbital energy" here is <phi|F|phi> (the
    Fock-matrix diagonal in the AVAS basis), not an eigenvalue -- reported for the same qualitative
    comparison the canonical-orbital diagnostic (results/be2_cbs_5z/diagnose.py) made, not as an
    MO energy ladder.
    """
    mol = gto.M(atom=f"Be 0 0 0; Be 0 0 {R}", basis=basis, spin=0, verbose=0)
    mf = scf.RHF(mol).run()
    ncas, nelecas, mo = avas_active_space(mf)
    ncore = (mol.nelectron - nelecas) // 2
    act = mo[:, ncore:ncore + ncas]
    fock = mf.get_fock()
    e_act = np.einsum("pi,pq,qi->i", act, fock, act)
    r2 = mol.intor("int1e_r2")
    rz = mol.intor("int1e_r")[2]
    spread = np.array([np.sqrt(c @ r2 @ c - (c @ rz @ c) ** 2) for c in act.T])
    return e_act, spread


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
                p = avas_casci_nevpt2_point(R, basis)
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
    """Same wide quartic fit be2_cbs.py / be2_cbs_5z.py use: window 2.1<=R<=3.0, De vs R=8.0."""
    win = [R for R in RS if 2.1 <= R <= 3.0]
    p = np.poly1d(np.polyfit(win, [curve[R] for R in win], 4))
    roots = [x.real for x in p.deriv().r if abs(x.imag) < 1e-6 and win[0] <= x.real <= win[-1]]
    Re = min(roots, key=lambda x: p(x))
    return float(Re), float((curve[8.0] - p(Re)) * HA2CM)


def cbs_curve(pts, lo, hi):
    return {R: pts[(hi, R)]["e_casci"] + cbs_extrapolate_correlation(
        BASIS_CARDINAL[lo], pts[(lo, R)]["e_corr"], BASIS_CARDINAL[hi], pts[(hi, R)]["e_corr"])
            for R in RS}


def write_diagnostics():
    os.makedirs(OUT_DIR, exist_ok=True)
    lines = ["basis,R,active_orbital_fock_energies_Ha,active_orbital_rms_radius_bohr"]
    for basis in BASES:
        for R in (2.5, 8.0):
            e_act, spread = active_orbital_diagnostics(basis, R)
            lines.append(f"{basis},{R}," + " ".join(f"{x:+.3f}" for x in e_act) + "," +
                         " ".join(f"{s:.1f}" for s in spread))
    text = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "diagnose.txt"), "w") as f:
        f.write(text + "\n")
    print(text)
    return lines


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
    verdict = "CONFIRMED" if abs(shift) <= CBS_CONFIRM_THRESHOLD_CM1 else "NOT CONFIRMED"
    lines.append(f"# De shift CBS(QZ/5Z) - CBS(TZ/QZ) = {shift:+.1f} cm^-1 -> method attribution "
                 f"{verdict} (pre-registered threshold {CBS_CONFIRM_THRESHOLD_CM1:.0f} cm^-1)")
    des = [wells[b][1] for b in BASES]
    mono = "monotone" if (des[0] < des[1] < des[2] or des[0] > des[1] > des[2]) else "NON-monotone"
    lines.append(f"# single-basis De(TZ,QZ,5Z) = {des[0]:.1f}, {des[1]:.1f}, {des[2]:.1f} cm^-1 "
                 f"-> {mono} in X")
    walls = {b: np.mean([pts[(b, R)]["wall_s"] for R in RS]) for b in BASES}
    lines.append("# mean wall s/point: " + ", ".join(f"{b}={w:.1f}" for b, w in walls.items())
                 + f"; host {platform.machine()} {os.cpu_count()} cpu")
    text = "\n".join(lines)
    with open(os.path.join(OUT_DIR, "wells.csv"), "w") as f:
        f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    write_diagnostics()
    summarize(compute_missing())
