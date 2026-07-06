#!/usr/bin/env python3
"""
Strain / pressure response of the Nb3X8 downfolded dimers -- Grueneisen parameters in the hopping.

Uniaxial compression of a breathing-kagome Nb3X8 layer increases the inter-layer hopping |t| (the
dimerization overlap) at roughly fixed on-site/inter-site Coulomb U0, Us. So d ln X / d ln |t| is
the leading strain response of any cluster observable X -- a Grueneisen parameter with |t| as the
pressure knob. Strain-tunable magnetism in Nb3Cl8 is an active experimental target (e.g. the
strain-tuned spin-liquid-candidate work, arXiv:2601.14524), which makes these signs and magnitudes
falsifiable predictions, not just internal numbers.

From the exact dimer (cRPA parameters of arXiv:2501.10320, via nb3x8_gaps / odmd_spin):

THE FINDINGS (specs/SPEC_nb3x8_strain.md):
  * SPIN gap J stiffens under compression EVERYWHERE (gamma_J > 0), and gamma_J runs monotonically
    from the atomic-limit value 2 (Nb3F8, where J ~ 4t^2/(U0-Us) so d ln J/d ln t = 2) toward the
    strong-hopping limit 1 (Nb3I8: 2.00 -> 1.89 -> 1.78 -> 1.52) -- the halide series traces the
    correlation crossover, and the closed-form dJ/dt matches central finite differences exactly.
  * CHARGE gap is NON-MONOTONIC in |t| (a minimum at |t*|), and the halide family STRADDLES it:
    Nb3F8/Nb3Cl8 sit below their minima (compression softens the charge gap, gamma_gap < 0),
    Nb3Br8/Nb3I8 sit above (compression stiffens it, gamma_gap > 0). So spin and charge respond in
    OPPOSITE directions only for the light halides; |t*| falls F->I (271->152->122->76 meV) as U0
    drops.
  * THE SHARP PREDICTION: Nb3Cl8 -- the most-studied member -- sits almost exactly AT its charge-gap
    minimum (|t| = 136 vs |t*| ~ 152 meV), so its strain response is SPIN-CHARGE DECOUPLED: strong
    spin-gap stiffening (gamma_J ~ 1.9) with a near-vanishing charge-gap response (|gamma_gap| <
    0.05), a > 30x split. Straining Nb3Cl8 should move its singlet-triplet gap (and, scaling with J,
    its chi(T) maximum and Schottky peak) while barely touching its Mott gap.

HONEST SCOPE: |t| is the sole strain proxy -- real strain also shifts U0, Us and the geometry, and
changes in-plane couplings this isolated dimer omits; a linear-response (Grueneisen) statement, not
a strain phase diagram; density-density only. Nb3F8's spin gap (J ~ 0.05 meV) is below the model's
neglected terms -- its gamma_J = 2 is the analytic atomic-limit exponent, not a magnitude claim.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar

from nb3x8_gaps import exact_charge_gap
from odmd_spin import dimer_exchange_analytic


def spin_gap_gruneisen(U0: float, t: float, Us: float) -> float:
    """gamma_J = d ln J / d ln |t| of the exact singlet-triplet gap, closed form. Runs from 2
    (atomic, J ~ 4t^2/(U0-Us)) to 1 (strong hopping, J ~ 2|t|)."""
    D = U0 - Us
    root = np.sqrt(0.25 * D * D + 4.0 * t * t)
    J = root - 0.5 * D
    dJdt = 4.0 * t / root
    return float(t * dJdt / J)


def charge_gap_gruneisen(U0: float, t: float, Us: float, rel: float = 1e-4) -> float:
    """gamma_gap = d ln(Delta_c) / d ln |t| of the exact Mott gap, central finite difference."""
    h = rel * abs(t)
    gp = exact_charge_gap(U0, t + h, Us)
    gm = exact_charge_gap(U0, t - h, Us)
    g0 = exact_charge_gap(U0, t, Us)
    return float(t * (gp - gm) / (2.0 * h) / g0)


def charge_gap_min_hopping(U0: float, Us: float, t_hi: float = 400.0) -> float:
    """|t*| (meV) minimizing the exact charge gap at fixed U0, Us -- the strain-response sign flip
    (below: gamma_gap < 0; above: gamma_gap > 0)."""
    res = minimize_scalar(lambda tt: exact_charge_gap(U0, -tt, Us),
                          bounds=(1.0, t_hi), method="bounded")
    return float(res.x)


if __name__ == "__main__":
    from nb3x8_gaps import NB3X8_LT_BULK
    from nb3x8_magnetometry import chi_max_temperature
    from nb3x8_thermo import schottky_peak_temperature

    def fd_gruneisen(f, p, rel=1e-5):
        h = rel * abs(p["t"])
        pu, pd = dict(p, t=p["t"] + h), dict(p, t=p["t"] - h)
        return p["t"] * (f(**pu) - f(**pd)) / (2 * h) / f(**p)

    print("Nb3X8 strain response -- Grueneisen parameters gamma = d ln X / d ln |t|")
    print(f"{'set':>7} | {'|t|':>6} | {'|t*|':>6} | {'gam_J':>6} | {'gam_J(fd)':>9} | "
          f"{'gam_gap':>8} | {'gam_chiTmax':>11} | {'gam_Schottky':>12}")
    for name, p in NB3X8_LT_BULK.items():
        gJ = spin_gap_gruneisen(**p)
        gJ_fd = fd_gruneisen(dimer_exchange_analytic, p)
        gg = charge_gap_gruneisen(**p)
        tstar = charge_gap_min_hopping(p["U0"], p["Us"])
        gchi = fd_gruneisen(chi_max_temperature, p) if abs(p["t"]) > 1 else float("nan")
        gsch = fd_gruneisen(schottky_peak_temperature, p) if abs(p["t"]) > 1 else float("nan")
        print(f"{name:>7} | {abs(p['t']):6.1f} | {tstar:6.1f} | {gJ:6.3f} | {gJ_fd:9.3f} | "
              f"{gg:8.3f} | {gchi:11.3f} | {gsch:12.3f}")
    print("\nFinding: spin gap stiffens everywhere (gamma_J: 2->1 across the series); charge gap is")
    print("non-monotonic (family straddles |t*|); Nb3Cl8 sits at its charge-gap minimum -> spin-")
    print("charge-decoupled strain response. See specs/SPEC_nb3x8_strain.md.")
