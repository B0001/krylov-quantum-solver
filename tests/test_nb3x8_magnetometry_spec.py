"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_magnetometry.md.

Claim: the ab-initio interlayer singlet-triplet gap J of the Nb3X8 bilayer dimers (odmd_spin /
nb3x8_susceptibility, from the cRPA parameters of arXiv:2501.10320), with NO fit parameters,
predicts the SCALE and ORDERING of the measured magnetic-singlet transitions of Nb3Cl8 (~90 K,
Sheckelton 2017) and Nb3Br8 (~382 K, Haraguchi 2017), and its quantified MISS localizes the
missing physics -- the isolated-cluster -> cooperative-lattice renormalization, and the fact that
the interlayer J is a different coupling from the weak in-plane exchange that sets theta_W.

Comparison against *measured* references (experiment is the ground truth here), not a fit. Isolated
bilayer dimer, density-density only. PySCF/qiskit, no block2; `make gates` runs it in its own
process.
"""
import numpy as np
import pytest

from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_magnetometry import (
    EXPERIMENT,
    MEV_PER_K,
    chi_max_temperature,
    chi_ht_curie_ratio,
    overprediction_factor,
    predicted_transition_K,
    theta_over_measured,
    z_eff_ht,
)
from nb3x8_susceptibility import curie_weiss_theta, ionic_singlet_energy, susceptibility
from odmd_spin import dimer_exchange_analytic

# Members with an established interlayer-singlet transition (the falsifiable comparison set).
_MEASURED = ("Nb3Cl8", "Nb3Br8")


def test_G1_prediction_machinery_is_anchored():
    """The predictor rests on validated dimer thermodynamics: chi*T -> 1/2 (S=1/2 pair Curie
    constant) deep in the spin window, theta_CW = -J/4 exactly, and the exact chi maximum sits at
    ~0.625 J (Bleaney-Bowers), so chi_max_temperature is a faithful transition-scale estimator."""
    for name, p in NB3X8_LT_BULK.items():
        J = dimer_exchange_analytic(**p)
        assert abs(curie_weiss_theta(J) - (-J / 4.0)) < 1e-12, name
        if J < 1.0:  # Nb3F8 (J ~ 0.05 meV) has no resolvable spin window -- excluded, not gated
            continue
        # deep-spin Curie constant of the coupled S=1/2 pair -- only where a clean window
        # J << T << E_s exists (guarded exactly as SPEC_nb3x8_susceptibility G2): the strongly
        # hybridized halides have E_s/J ~ 3-10, so most have no such window.
        E_s = ionic_singlet_energy(**p)
        T_plat = 100.0 * J  # J/T = 0.01
        if T_plat <= E_s / 20.0:
            assert abs(susceptibility(**p, T=T_plat) * T_plat - 0.5) < 2e-3, name
        # exact chi maximum tracks the Bleaney-Bowers 0.625 J to a few percent
        Tmax = chi_max_temperature(**p)
        assert abs(Tmax / (0.625 * J) - 1.0) < 0.05, (name, Tmax, 0.625 * J)


def test_G2_reproduces_measured_scale_and_ordering():
    """THE WIN: the parameter-free prediction lands within an order of magnitude of every measured
    Tc, and reproduces the observed Cl < Br ordering."""
    for name in _MEASURED:
        predK = predicted_transition_K(**NB3X8_LT_BULK[name])
        obs = EXPERIMENT[name]["Tc_K"]
        assert 1.0 < predK / obs < 10.0, (name, predK, obs)          # right scale, not a fit
    # ordering preserved: predicted and observed both increase Cl -> Br
    assert (predicted_transition_K(**NB3X8_LT_BULK["Nb3Cl8"])
            < predicted_transition_K(**NB3X8_LT_BULK["Nb3Br8"]))
    assert EXPERIMENT["Nb3Cl8"]["Tc_K"] < EXPERIMENT["Nb3Br8"]["Tc_K"]


def test_G3_the_finding_overprediction_and_two_couplings():
    """DEFINITION OF DONE (the finding). (a) The isolated dimer OVERPREDICTS Tc for both halides
    (factor > 2), and the overprediction WEAKENS monotonically Cl -> Br (the cluster -> lattice
    renormalization). (b) The interlayer J does NOT set theta_W: -J/4 overshoots the measured
    theta_W of Nb3Cl8 by > 3x -- interlayer J sets Tc, a separate weak in-plane exchange sets
    theta_W.

    Threshold note (bd chem-jiy, 2026-09-30): this gate originally asserted > 5x against a
    theta_W = -13.1 K that does not appear anywhere in the cited primary source (Sheckelton et al.,
    arXiv:1701.05528) -- it reached this repo via an unverified search-engine summary. The
    primary-source-verified value is theta_W = -51.2 K (p.4, "T > 140 K" fit), which gives a real
    but smaller ratio of 3.75x (see EXPERIMENT and test_G5 below). The bound is corrected here to
    match the verified number, not loosened to fit a convenient result -- see
    specs/SPEC_nb3x8_magnetometry.md and specs/BACKLOG.md."""
    over_cl = overprediction_factor("Nb3Cl8")
    over_br = overprediction_factor("Nb3Br8")
    assert over_cl > 2.0 and over_br > 2.0, (over_cl, over_br)
    assert over_cl > over_br, (over_cl, over_br)                     # weakens down the series
    # the two-coupling separation: -J/4 is far from the measured theta_W
    assert abs(theta_over_measured("Nb3Cl8")) > 3.0, theta_over_measured("Nb3Cl8")


def test_G4_honest_boundary_iodide_excluded_and_scale_only():
    """The honest boundary: Nb3I8 has no interlayer-singlet transition, so it is NOT in the
    experimental comparison set; and the predictor is a *scale* estimator (kelvin, positive,
    monotone in J), never a claim to reproduce a first-order cooperative transition."""
    assert "Nb3I8" not in EXPERIMENT
    # scale estimator is well-defined and monotone in J across the family
    prevK = -np.inf
    for name in ("Nb3Cl8", "Nb3Br8", "Nb3I8"):  # increasing J
        predK = predicted_transition_K(**NB3X8_LT_BULK[name])
        assert predK > 0.0 and np.isfinite(predK), name
        assert predK > prevK, (name, predK, prevK)                  # monotone in J
        prevK = predK
    assert MEV_PER_K == 0.0861733


def test_G5_phase_correction_does_not_survive_primary_source_verification():
    """bd chem-jiy -- REFINEMENT of G3(b). specs/BACKLOG.md's Nb3X8 entry conjectured that the
    Curie-Weiss miss (interlayer J vs measured theta_W) is a phase-assignment artifact: theta_W is
    fitted in the HT (undimerized) phase, but this module compared it to a LT-phase J; swapping in
    the HT-phase J (nb3x8_gaps.NB3X8_HT_BULK, already unused) was claimed to invert the miss to a
    clean mean-field z_eff = 4.1 (theta = -z*J/4), using theta_W = -13.1 K.

    Primary-source check (arXiv:1701.05528 p.4): "an analysis of the inverse susceptibility data
    for T > 140 K yields a Curie constant C = 0.484 emu K mol^-1 Oe^-1 (p_eff = 1.97, consistent
    with S_eff = 1/2) and a Weiss temperature of theta = -51.2 K." C and mu_eff were correctly
    quoted; theta_W = -13.1 K does NOT appear anywhere in the paper.

    Killable checks (bd chem-jiy, pre-registered): dies if z_eff is unphysical (<1 or >12 for a
    layered stack) or if chi_HT(200K) misses the measured Curie constant by >2x -- either restores
    the original (now-corrected) miss as real, not a phase artifact.

    RESULT: z_eff_ht (recomputed with the VERIFIED theta_W = -51.2 K) is 16.1 -- OUTSIDE the [1,12]
    bound. The phase-correction hypothesis is KILLED: it was an artifact of the wrong, unverified
    -13.1 K datum, not a real phase-assignment fix. chi_HT(200K) alone is within the 2x sanity
    bound (does not independently kill), but the z_eff bound is sufficient and z_eff is reported as
    an extracted diagnostic, never a free-standing prediction."""
    z = z_eff_ht("Nb3Cl8")
    assert EXPERIMENT["Nb3Cl8"]["theta_K"] == -51.2, "theta_W must be the verified primary-source value"
    assert EXPERIMENT["Nb3Cl8"]["C_emu"] == 0.484 and EXPERIMENT["Nb3Cl8"]["mu_eff"] == 1.97
    assert not (1.0 <= z <= 12.0), z                       # UNPHYSICAL: kills the phase-artifact claim
    assert 10.0 < z < 20.0, z                               # pins the actual (unphysical) value found
    chi_ratio = chi_ht_curie_ratio("Nb3Cl8")
    assert 0.5 < chi_ratio < 2.0, chi_ratio                 # this check alone does not kill
    # the verified, corrected miss (LT J vs the verified theta_W) is real but smaller than the
    # originally-published, uncited 15x
    assert 3.0 < abs(theta_over_measured("Nb3Cl8")) < 5.0, theta_over_measured("Nb3Cl8")
    for fn in (z_eff_ht, chi_ht_curie_ratio):   # no primary-source HT datum for Br: ValueError, not TypeError
        with pytest.raises(ValueError, match="Nb3Br8"):
            fn("Nb3Br8")
