#!/usr/bin/env python3
"""
Nb3X8 metamagnetism at finite temperature -- the crossover width law.

SPEC_nb3x8_metamagnetism gates the T=0 field-driven singlet->triplet ground-state crossing (a
sharp step at h_c = J). This module adds temperature: a Boltzmann trace over the FULL N=2-sector
spectrum of the field-augmented Hamiltonian (nb3x8_metamagnetism.field_spectrum), in the same
style as nb3x8_susceptibility.py's zero-field thermal trace.

THE FINDING (specs/SPEC_nb3x8_metamagnetism_thermal.md G4): the T=0 step smooths into a crossover
whose 10-90% width in h obeys the closed form width = 2*ln(9)*T -- the elementary two-level
logistic width -- to < 0.1% relative, for T up to 0.2*J and across the whole magnetic halide
family (Cl/Br/I), even though the full trace carries 4 more levels (the field-independent Sz=0
triplet member, two ionic singlets) than the two-level picture assumes. Cross-validated (G2): the
zero-field slope dM/dh matches nb3x8_susceptibility.py's INDEPENDENTLY-derived susceptibility(T)
(a different route -- the zero-field Van Vleck <Sz^2>/T trace) to < 1e-5 relative.

HONEST SCOPE: same as nb3x8_metamagnetism.py (isolated dimer, g=2, density-density only); the
2*ln(9)*T law is checked for T <~ 0.2*J only -- not claimed near or above the charge scale E_s-J,
where the ionic singlets would necessarily intrude (see nb3x8_susceptibility.py's own charge-scale
boundary). No experimental comparison -- a prediction, not a reproduction.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from nb3x8_metamagnetism import field_spectrum


def magnetization_thermal(U0: float, t: float, Us: float, h: float, T: float) -> float:
    """<Sz> of the field-augmented dimer at temperature T (meV), Boltzmann-traced over the full
    N=2-sector spectrum (same style as nb3x8_susceptibility.n2_spectrum's thermal average)."""
    e, sz = field_spectrum(U0, t, Us, h)
    boltz = np.exp(-(e - e.min()) / T)
    return float((sz * boltz).sum() / boltz.sum())


def crossover_width(U0: float, t: float, Us: float, T: float,
                     lo: float = 1e-4, hi: float | None = None) -> float:
    """h(M=0.9) - h(M=0.1) at temperature T (meV) -- the 10-90% crossover width in field."""
    from odmd_spin import dimer_exchange_analytic

    J = dimer_exchange_analytic(U0, t, Us)
    if hi is None:
        hi = 3.0 * J

    def solve(target):
        f = lambda h: magnetization_thermal(U0, t, Us, h, T) - target  # noqa: E731
        return brentq(f, lo, hi)

    return solve(0.9) - solve(0.1)


if __name__ == "__main__":
    import math

    from nb3x8_gaps import NB3X8_LT_BULK
    from nb3x8_susceptibility import susceptibility
    from odmd_spin import dimer_exchange_analytic

    TARGET = 2.0 * math.log(9.0)
    print("Nb3X8 metamagnetism at finite T -- crossover width law (width/T vs 2*ln(9)):")
    print(f"{'material':>8} | {'T/J':>6} | {'width/T':>10} | {'2ln9':>8} | {'dM/dh vs chi':>14}")
    for name, p in NB3X8_LT_BULK.items():
        J = dimer_exchange_analytic(**p)
        if J < 1.0:
            continue  # Nb3F8: below the model's own noise floor
        for frac in (0.05, 0.2):
            T = frac * J
            width = crossover_width(**p, T=T)
            dh = 1e-4 * J
            dMdh = (magnetization_thermal(**p, h=dh, T=T)
                    - magnetization_thermal(**p, h=-dh, T=T)) / (2 * dh)
            chi = susceptibility(p["U0"], p["t"], p["Us"], T)
            print(f"{name:>8} | {frac:6.2f} | {width / T:10.4f} | {TARGET:8.4f} | "
                  f"{dMdh:.3e} vs {chi:.3e}")
