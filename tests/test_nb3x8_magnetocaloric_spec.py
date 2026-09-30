"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_magnetocaloric.md (Nb3X8 magnetocaloric S(T,B)).

Test-first: ``nb3x8_magnetocaloric`` does not exist yet, so this file is RED until the spec is
implemented. Two independent claims: (i) MACHINERY -- the Maxwell relation (dS/dB)_T = (dM/dT)_B,
computed from two independently-derived Boltzmann traces (entropy's ln Z + <E>/T route here vs
magnetization_thermal's <Sz> route in nb3x8_metamagnetism_thermal.py), agrees to machine precision;
(ii) VERDICT -- |Delta S_M| stays far under 10% of R ln2 per formula unit for B<=100T, T>=2K, and
far under GGG's theoretical entropy ceiling in matching (J/(kg K)) units.

PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import math

import numpy as np
import scipy.constants as sc

from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_magnetocaloric import (
    GGG_ENTROPY_CEILING_J_PER_KG_K,
    GGG_MOLAR_MASS_G_PER_MOL,
    MAGNETIC_HALIDES,
    NB3X8_FORMULA_MASS_G_PER_MOL,
    delta_entropy_isothermal,
    delta_s_m_fraction_of_r_ln2,
    delta_s_m_j_per_kg_k,
    dm_dT_numeric,
    ds_dh_numeric,
    scan_max_delta_s,
)
from nb3x8_metamagnetism import G_MU_B, critical_field_tesla
from odmd_spin import dimer_exchange_analytic

# The swept spin window for G1 -- see SPEC_nb3x8_magnetocaloric.md G1 scope note: points near
# T/J < 0.1 or h/J very close to an exact crossing/edge (INCLUDING h/J = 0 exactly, where both
# (dS/dh)_T and (dM/dT)_h vanish identically by the h -> -h symmetry, making a relative comparison
# a 0/0) were empirically found to be ill-conditioned for a finite-difference relative comparison,
# not a machinery bug -- excluded here, not hidden.
T_FRACS = (0.1, 0.15, 0.2, 0.3, 0.5, 0.8)
H_FRACS = (-0.3, 0.05, 0.3, 0.6, 0.9, 1.1, 1.4, 1.8)


def test_G1_maxwell_relation_definition_of_done():
    """DEFINITION OF DONE: (dS/dh)_T (entropy route) and (dM/dT)_h (magnetization_thermal route,
    an independently-derived module) agree to < 1e-6 relative across the documented spin window."""
    worst = 0.0
    for name in MAGNETIC_HALIDES:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        for Tfrac in T_FRACS:
            T = Tfrac * J
            dT = 1e-4 * T
            for hfrac in H_FRACS:
                h = hfrac * J
                dh = 1e-4 * J
                dSdh = ds_dh_numeric(**p, h=h, T=T, dh=dh)
                dMdT = dm_dT_numeric(**p, h=h, T=T, dT=dT)
                rel = abs(dSdh / dMdT - 1.0)
                worst = max(worst, rel)
                assert rel < 1e-6, (name, Tfrac, hfrac, dSdh, dMdT, rel)
    assert worst < 1e-5, worst  # sanity: the whole window is comfortably inside the gate


def test_G2_peak_is_near_Bc():
    """Supports the verdict's premise: |Delta S_M| is non-decreasing in B up to 100 T (the search
    max legitimately sits at the boundary), and is far larger at B_c than at 100 T."""
    for name in MAGNETIC_HALIDES:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        T = 0.3 * J
        Bs = np.linspace(0.0, 100.0, 15)
        vals = [abs(delta_entropy_isothermal(**p, T=T, h=B * G_MU_B)) for B in Bs]
        assert all(b >= a - 1e-12 for a, b in zip(vals, vals[1:])), (name, vals)

        T_low = 0.05 * J
        Bc = critical_field_tesla(**p)
        dS_100T = abs(delta_entropy_isothermal(**p, T=T_low, h=100.0 * G_MU_B))
        dS_Bc = abs(delta_entropy_isothermal(**p, T=T_low, h=Bc * G_MU_B))
        assert dS_Bc > 10.0 * dS_100T, (name, dS_100T, dS_Bc)


def test_G3_ten_percent_kill_criterion_verdict():
    """VERDICT -- take seriously if this fails: for each magnetic halide, the max |Delta S_M| over
    B in [0,100]T, T in [2K, 5*J/k_B] stays under 10% of R ln2 per formula unit."""
    for name in MAGNETIC_HALIDES:
        p = NB3X8_LT_BULK[name]
        max_dS, B_star, T_star = scan_max_delta_s(**p)
        frac = delta_s_m_fraction_of_r_ln2(**p, T_kelvin=T_star, B_tesla=B_star)
        assert 0.0 <= B_star <= 100.0 and T_star >= 2.0, (name, B_star, T_star)
        assert frac < 0.10, (name, frac, B_star, T_star, max_dS)
        # record how far under the threshold we actually are (2+ orders of magnitude expected)
        assert frac < 0.05, (name, frac)  # tighter check than the bare kill criterion


def test_G4_sourced_unit_matched_ggg_comparison():
    """Molar masses match hand-checked values (guards silent drift); the GGG ceiling matches its
    closed form; Nb3Cl8's best-case Delta S_M is a small fraction of that ceiling."""
    # hand-checked against IUPAC standard atomic weights (Nb 92.906, Cl 35.45, Br 79.904,
    # I 126.904, Gd 157.25, Ga 69.723, O 15.999 g/mol)
    assert abs(NB3X8_FORMULA_MASS_G_PER_MOL["Nb3Cl8"] - 562.32) < 0.1
    assert abs(NB3X8_FORMULA_MASS_G_PER_MOL["Nb3Br8"] - 917.95) < 0.1
    assert abs(NB3X8_FORMULA_MASS_G_PER_MOL["Nb3I8"] - 1293.95) < 0.1
    assert abs(GGG_MOLAR_MASS_G_PER_MOL - 1012.35) < 0.1

    expected_ceiling = 3.0 * sc.R * math.log(8.0) * 1000.0 / GGG_MOLAR_MASS_G_PER_MOL
    assert abs(GGG_ENTROPY_CEILING_J_PER_KG_K / expected_ceiling - 1.0) < 1e-9
    assert abs(GGG_ENTROPY_CEILING_J_PER_KG_K - 51.2) < 0.5

    p = NB3X8_LT_BULK["Nb3Cl8"]
    max_dS, B_star, T_star = scan_max_delta_s(**p)
    dS_phys = delta_s_m_j_per_kg_k(**p, T_kelvin=T_star, B_tesla=B_star,
                                   formula_mass_g_per_mol=NB3X8_FORMULA_MASS_G_PER_MOL["Nb3Cl8"])
    ratio = abs(dS_phys) / GGG_ENTROPY_CEILING_J_PER_KG_K
    assert ratio < 0.01, (dS_phys, GGG_ENTROPY_CEILING_J_PER_KG_K, ratio)
