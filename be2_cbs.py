#!/usr/bin/env python3
"""
be2_cbs.py -- Be2 well depth toward experiment: core-valence dynamic correlation on top of the
validated CAS(4,8) reference, extrapolated across cc-pVTZ/QZ to the complete basis set (CBS) limit.

BACKLOG CLAIM (specs/BACKLOG.md): "core-valence correlation + a cc-pVxZ->CBS extrapolation moves
the FCI/DMRG well depth from ~305 cm^-1 toward the experimental 929.7." Gate as originally stated:
`|D_e - 930| < 100 cm^-1` at CBS+CV, `R_e` within 0.1 A of 2.45.

METHOD: study_be2.py's active space (CASCI(4,8), the valence 2s/2p of both Be atoms, frozen-core
FCI) is the validated static-correlation reference -- it is NOT reoptimized here (see R1 in
specs/SPEC_be2_cbs.md: CASSCF orbital optimization is numerically unstable for this near-degenerate
system, confirmed by conv=False at several bond lengths/bases during development; CASCI on fixed
canonical HF orbitals converges cleanly everywhere it's used). NEVPT2 supplies the dynamic
correlation on top of that CASCI reference, and -- because the two Be 1s orbitals are NOT frozen
for the perturbative step -- it includes core-valence and core-core correlation automatically, no
separate all-electron active space needed. The NEVPT2 correlation energy is CBS-extrapolated via
the standard two-point Helgaker form E_corr(X) = E_CBS + B/X^3 (X = cardinal number, TZ=3/QZ=4);
the CASCI reference energy is taken at the larger (QZ) basis, since it is far closer to basis
saturation than the correlation energy (DZ->TZ->QZ CASCI-alone shifts are ~2 mHa then ~0.7 mHa at
R=2.45, an order of magnitude below the correlation-energy basis dependence).

HONEST FINDING (see specs/SPEC_be2_cbs.md G4): this composition does NOT clear the backlog's
original tolerance. It reproduces a genuine bound well at the right general location -- unlike the
frozen-core CAS(4,8)/cc-pVDZ FCI baseline in study_be2.py, whose ~305 cm^-1 "well" is a spurious
long-range artifact at R~4.5 A with no real minimum near R_e -- but underbinds by roughly half
(D_e ~ 460-470 cm^-1 vs 929.7) and overshoots R_e by ~0.15 A (~2.6 A vs 2.45 A). The residual is
attributed to two known limitations of this composition, not tuned away: (i) NEVPT2 is a
second-order perturbative correction, weaker than the CCSD(T)-F12/large-active-space MRCI+Q
treatments the literature needed to match experiment; (ii) the 8-orbital valence active space
excludes higher Rydberg-like virtuals that contribute non-negligibly to Be2's unusually
correlation-driven bond. Both are recorded as follow-up scope, not silently absorbed into a passing
gate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Sequence, Tuple

import numpy as np
from pyscf import gto, mcscf, mrpt, scf

HA2CM = 219474.6313702
BASIS_CARDINAL = {"ccpvdz": 2, "ccpvtz": 3, "ccpvqz": 4}


@dataclass
class Be2Point:
    R: float
    basis: str
    e_casci: float
    e_corr: float

    @property
    def e_tot(self) -> float:
        return self.e_casci + self.e_corr


def casci_nevpt2_point(R: float, basis: str, cas_electrons: int = 4,
                       cas_orbitals: int = 8) -> Be2Point:
    """CASCI(cas_electrons, cas_orbitals) + NEVPT2 total energy at bond length R (Angstrom).

    Fixed canonical HF orbitals (no CASSCF reorientation -- see module docstring R1). The two Be
    1s orbitals stay doubly occupied outside the active space and unfrozen for NEVPT2, so its
    perturbative correction includes core-valence/core-core dynamic correlation automatically.
    """
    mol = gto.M(atom=f"Be 0 0 0; Be 0 0 {R}", basis=basis, spin=0, verbose=0)
    mf = scf.RHF(mol).run()
    mc = mcscf.CASCI(mf, cas_orbitals, cas_electrons)
    mc.kernel()
    e_corr = float(mrpt.NEVPT(mc).kernel())
    return Be2Point(R=R, basis=basis, e_casci=float(mc.e_tot), e_corr=e_corr)


def cbs_extrapolate_correlation(x_lo: int, e_lo: float, x_hi: int, e_hi: float) -> float:
    """Two-point Helgaker CBS extrapolation of a correlation energy: E(X) = E_CBS + B/X^3."""
    return float((e_hi * x_hi ** 3 - e_lo * x_lo ** 3) / (x_hi ** 3 - x_lo ** 3))


def cbs_point(R: float, lo_basis: str = "ccpvtz", hi_basis: str = "ccpvqz",
             cas_electrons: int = 4, cas_orbitals: int = 8) -> float:
    """CBS(lo_basis/hi_basis)-extrapolated NEVPT2 correlation + CASCI reference at hi_basis."""
    lo = casci_nevpt2_point(R, lo_basis, cas_electrons, cas_orbitals)
    hi = casci_nevpt2_point(R, hi_basis, cas_electrons, cas_orbitals)
    e_corr_cbs = cbs_extrapolate_correlation(BASIS_CARDINAL[lo_basis], lo.e_corr,
                                             BASIS_CARDINAL[hi_basis], hi.e_corr)
    return hi.e_casci + e_corr_cbs


def quadratic_well(Rs: Sequence[float], Es: Sequence[float], r_asymptote: float,
                   asymptote_energy: float) -> Tuple[float, float]:
    """Vertex (Re, De) in (Angstrom, cm^-1) of a quadratic fit through the (Rs, Es) points.

    ``Es`` must bracket a local minimum (the well); ``asymptote_energy`` is the separately
    computed dissociated-limit energy (Hartree) that defines De = E(asymptote) - E(Re).
    """
    Rs = np.asarray(Rs, dtype=float)
    Es = np.asarray(Es, dtype=float)
    a, b, c = np.polyfit(Rs, Es, 2)
    Re = float(-b / (2 * a))
    Emin = a * Re ** 2 + b * Re + c
    De = float((asymptote_energy - Emin) * HA2CM)
    return Re, De


if __name__ == "__main__":
    import csv
    import os

    # The full record: 13-point curve x {ccpvdz, ccpvtz, ccpvqz}, quartic-fit well from a wide
    # window -- the careful version of the numbers the gates pin cheaply from 4 points.
    Rs = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 4.0, 6.0, 8.0]
    bases = ["ccpvdz", "ccpvtz", "ccpvqz"]
    os.makedirs("data", exist_ok=True)
    rows = []
    points: Dict[str, Dict[float, Be2Point]] = {b: {} for b in bases}
    for basis in bases:
        for R in Rs:
            p = casci_nevpt2_point(R, basis)
            points[basis][R] = p
            rows.append(dict(R=R, basis=basis, e_casci=p.e_casci, e_corr=p.e_corr,
                             e_tot=p.e_tot))
            print(f"{basis} R={R:4.2f}  E_casci={p.e_casci:.6f}  E_corr={p.e_corr:.6f}")
    with open("data/be2_cbs_curve.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    cbs_curve = {}
    for R in Rs:
        lo, hi = points["ccpvtz"][R], points["ccpvqz"][R]
        e_corr_cbs = cbs_extrapolate_correlation(3, lo.e_corr, 4, hi.e_corr)
        cbs_curve[R] = hi.e_casci + e_corr_cbs
    win = [R for R in Rs if 2.1 <= R <= 3.0]
    coeffs = np.polyfit(win, [cbs_curve[R] for R in win], 4)
    p = np.poly1d(coeffs)
    roots = [x.real for x in p.deriv().r if abs(x.imag) < 1e-6 and win[0] <= x.real <= win[-1]]
    Re = min(roots, key=lambda x: p(x))
    De = (cbs_curve[8.0] - p(Re)) * HA2CM
    print(f"\nCBS(TZ/QZ) + CV, wide quartic fit: Re={Re:.4f} A  De={De:.1f} cm^-1")
    print("Experiment (Merritt, Bondybey & Heaven, Science 2009): Re=2.45 A  De=929.7 cm^-1")
