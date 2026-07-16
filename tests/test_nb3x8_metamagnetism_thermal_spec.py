"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_metamagnetism_thermal.md (finite-T magnetization,
the crossover width law).

Test-first: ``nb3x8_metamagnetism_thermal`` does not exist yet, so this file is RED until the
spec is implemented. Claim: the full N=2-sector Boltzmann trace under the Zeeman term (a) recovers
SPEC_nb3x8_metamagnetism's sharp T=0 step as T->0, (b) has its zero-field slope dM/dh match the
INDEPENDENTLY-BUILT nb3x8_susceptibility.susceptibility(T) -- a real cross-check between two
separately-derived modules, not the same code twice -- and (c) shows the step smoothing into a
crossover whose 10-90% width in h obeys the closed form width = 2*ln(9)*T, material-independent.

PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import math

import numpy as np

from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_metamagnetism import magnetization as magnetization_zero_T
from nb3x8_metamagnetism_thermal import crossover_width, magnetization_thermal
from nb3x8_susceptibility import susceptibility
from odmd_spin import dimer_exchange_analytic

MAGNETIC = ("Nb3Cl8", "Nb3Br8", "Nb3I8")  # Nb3F8 excluded: J below the model's noise floor


def test_G1_low_T_recovers_the_sharp_step():
    """At T = 1e-4*J, the thermal trace matches the T=0 step (SPEC_nb3x8_metamagnetism) closely."""
    for name in MAGNETIC:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        T = 1e-4 * J
        m_below = magnetization_thermal(**p, h=0.9 * J, T=T)
        m_above = magnetization_thermal(**p, h=1.1 * J, T=T)
        assert m_below < 1e-3, (name, m_below)
        assert m_above > 1 - 1e-3, (name, m_above)
        assert abs(m_below - magnetization_zero_T(**p, h=0.9 * J)) < 1e-3, name
        assert abs(m_above - magnetization_zero_T(**p, h=1.1 * J)) < 1e-3, name


def test_G2_matches_independently_built_susceptibility():
    """DEFINITION OF DONE: central finite-difference dM/dh|_{h=0} at T=0.1*J matches
    nb3x8_susceptibility.susceptibility(T) -- built from a different derivation route (the
    zero-field Van Vleck <Sz^2>/T trace) -- to < 1e-5 relative."""
    for name in MAGNETIC:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        T = 0.1 * J
        dh = 1e-4 * J
        dMdh = (magnetization_thermal(**p, h=dh, T=T)
                - magnetization_thermal(**p, h=-dh, T=T)) / (2 * dh)
        chi = susceptibility(p["U0"], p["t"], p["Us"], T)
        assert abs(dMdh / chi - 1.0) < 1e-5, (name, dMdh, chi)


def test_G3_monotonic_in_field():
    """magnetization_thermal is non-decreasing in h across the crossing, at two temperatures."""
    for name in MAGNETIC:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        for T in (0.02 * J, 0.15 * J):
            hs = np.linspace(-0.5 * J, 2.0 * J, 25)
            ms = [magnetization_thermal(**p, h=h, T=T) for h in hs]
            assert all(b >= a - 1e-9 for a, b in zip(ms, ms[1:])), (name, T, ms)


def test_G4_crossover_width_law_and_its_breakdown():
    """width(10-90%) / T == 2*ln(9) (~4.394) to < 1e-4 relative for T <~ 0.1*J, across materials --
    the full 6-level trace reduces to the elementary two-level logistic width in this regime.
    THE BOUNDARY: by T = 0.2*J the deviation has already crossed 0.1%, and by T = 0.3*J it exceeds
    1% -- the third level (the field-independent Sz=0 triplet member) measurably contaminates the
    two-level picture once T is no longer small compared to J. Gated as a finding, not hidden."""
    target = 2.0 * math.log(9.0)
    for name in MAGNETIC:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        for frac in (0.01, 0.05, 0.1):
            T = frac * J
            width = crossover_width(p["U0"], p["t"], p["Us"], T)
            assert abs(width / T / target - 1.0) < 1e-4, (name, frac, width / T, target)

        # the boundary: the law measurably breaks down well before T ~ J
        width_02 = crossover_width(p["U0"], p["t"], p["Us"], 0.2 * J)
        width_03 = crossover_width(p["U0"], p["t"], p["Us"], 0.3 * J)
        dev_02 = abs(width_02 / (0.2 * J) / target - 1.0)
        dev_03 = abs(width_03 / (0.3 * J) / target - 1.0)
        assert dev_02 > 1e-3, (name, dev_02)
        assert dev_03 > 1e-2, (name, dev_03)
        assert dev_03 > dev_02, (name, dev_02, dev_03)  # deviation grows monotonically with T
