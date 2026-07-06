"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_thermo.md (magnetic heat capacity & entropy).

Claim: the exact-spectrum magnetic C_m(T) and S_m(T) of the downfolded Nb3X8 dimer show (a) a
Schottky anomaly pinned at T ~ 0.352 J by the analytic two-level (singlet, triplet) result -- a
J-scale fingerprint whose ratio to the chi(T) peak (0.625 J) is universal -- and (b) a localized-
moment entropy plateau R ln 4/dimer (R ln 2/cluster) that is clean ONLY when the charge scale
E_s >> J, degrading monotonically with E_s/J and vanishing for the iodide (the same charge-scale
boundary as the Bleaney-Bowers regime in nb3x8_susceptibility).

REPRODUCTION of textbook magnetothermodynamics with ab-initio parameters. Isolated single dimer,
density-density only. PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np

from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_magnetometry import chi_max_temperature
from nb3x8_susceptibility import ionic_singlet_energy
from nb3x8_thermo import (
    entropy,
    entropy_plateau,
    heat_capacity,
    schottky_peak_temperature,
    two_level_schottky_ratio,
)
from odmd_spin import dimer_exchange_analytic

LN4 = np.log(4.0)
LN6 = np.log(6.0)
# J-resolvable members (Nb3F8 J~0.05 meV has no spin feature -- excluded), ordered by E_s/J.
_RESOLVED = ("Nb3Cl8", "Nb3Br8", "Nb3I8")


def test_G1_schottky_peak_pinned_by_two_level():
    """The analytic (1,3) two-level Schottky peak is at T/J ~ 0.3515, and each dimer's C_m(T) peaks
    there within 3%. C >= 0 everywhere; third-law S(T->0)=0 (nondegenerate singlet ground)."""
    r = two_level_schottky_ratio(3.0)
    assert abs(r - 0.3515) < 0.01, r
    for name in _RESOLVED:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        Tpk = schottky_peak_temperature(**p)
        assert abs(Tpk / J - r) < 0.03 * r + 0.01, (name, Tpk / J, r)
        grid = np.linspace(0.05 * J, 30.0 * J, 500)
        assert np.all(heat_capacity(**p, T=grid) > -1e-12), name          # C >= 0
        assert entropy(**p, T=1e-4 * J) < 1e-6, name                      # S(0) -> 0


def test_G2_localized_moment_entropy_plateau():
    """The well-separated Nb3Cl8 (E_s/J ~ 17) has a genuine R ln 4/dimer (= R ln 2/cluster) entropy
    plateau: flat (min |dS/dlnT| < 0.10) and within 3% of ln 4; the sector entropy saturates at
    ln 6 at high T."""
    p = NB3X8_LT_BULK["Nb3Cl8"]
    Splat, flat = entropy_plateau(**p)
    assert flat < 0.10, flat
    assert abs(Splat - LN4) < 0.03 * LN4, (Splat, LN4)
    E_s = ionic_singlet_energy(**p)
    assert abs(entropy(**p, T=50.0 * E_s) - LN6) < 1e-3                    # full-sector saturation
    assert abs(0.5 * Splat - np.log(2.0)) < 0.03 * np.log(2.0)            # R ln 2 per cluster


def test_G3_plateau_cleanliness_is_charge_scale_set():
    """DEFINITION OF DONE (the finding): the entropy plateau degrades monotonically as E_s/J
    shrinks -- both its flatness metric and its deviation from ln 4 increase strictly Cl -> Br -> I
    -- and the iodide (E_s/J ~ 3) has NO clean plateau (deviation > 10% of ln 4). Same charge-scale
    boundary as chi(T)'s Bleaney-Bowers breakdown."""
    es_over_j, flats, devs = [], [], []
    for name in _RESOLVED:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        Splat, flat = entropy_plateau(**p)
        es_over_j.append(ionic_singlet_energy(**p) / J)
        flats.append(flat)
        devs.append(abs(Splat - LN4) / LN4)
    # E_s/J decreases Cl -> Br -> I; flatness and |dev| increase strictly (worse plateau)
    assert es_over_j[0] > es_over_j[1] > es_over_j[2], es_over_j
    assert flats[0] < flats[1] < flats[2], flats
    assert devs[0] < devs[1] < devs[2], devs
    assert devs[2] > 0.10, devs[2]                                        # iodide: no clean plateau
    assert devs[0] < 0.05, devs[0]                                        # chloride: clean plateau


def test_G4_cross_tie_to_chi_and_scope():
    """The heat-capacity peak and the chi(T) maximum are both exact two-level features, so their
    ratio is universal and material-independent: C-peak/chi-peak ~ 0.56 across the family (within
    5%). Nb3F8 (J~0) carries no resolvable spin feature."""
    r_expected = two_level_schottky_ratio(3.0) / 0.625                    # ~0.562
    for name in _RESOLVED:
        p = NB3X8_LT_BULK[name]
        ratio = schottky_peak_temperature(**p) / chi_max_temperature(**p)
        assert abs(ratio / r_expected - 1.0) < 0.05, (name, ratio, r_expected)
    assert dimer_exchange_analytic(**NB3X8_LT_BULK["Nb3F8"]) < 1.0        # F excluded (J~0)
