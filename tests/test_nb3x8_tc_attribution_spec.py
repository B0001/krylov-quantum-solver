"""
Acceptance gates G1-G7 for specs/SPEC_nb3x8_tc_attribution.md (beads chem-phh, chem-cfl).

chem-g78 killed the "coordination/mean-field reduction" attribution for the Nb3X8 5.3x (Cl) /
2.3x (Br) Tc overprediction (tests/test_nb3x8_gaps_spec.py::test_G8_*), leaving the attribution
open with three named candidates (SPEC_nb3x8_magnetometry.md section 7): the cooperative/first-
order structural transition, in-plane kagome exchange, and a direct lattice renormalization of
t_s_perp. This module pre-registers a killable check for each. chem-phh executed candidate 3 from
data already vendored in this repo; chem-cfl executed candidate 1, sourcing elastic/lattice data
from a second primary reference (arXiv:2412.03360) this sandbox could fetch via WebFetch (unlike
WebSearch).

Gates:
  * G1 -- machinery anchor: at eta=1 (no rescaling) the module's own predicted-Tc machinery
    reproduces the already-published, independently-gated overprediction factors
    (nb3x8_magnetometry.overprediction_factor: 5.3x Cl / 2.3x Br) -- if this drifts, nothing
    below can be trusted.
  * G2 -- required_t_rescale is self-consistent: plugging eta back into predicted_transition_K
    reproduces the measured Tc to <1e-6 relative, for both halides, and eta in (0, 1) (a
    reduction, matching the overprediction's sign).
  * G3 -- candidate 3 (DEFINITION OF DONE): the rescaling required to match measured Tc is more
    than 10x the largest same-structural-motif |t_s_perp| spread this family's own vendored data
    shows (bulk vs bilayer, both LT/dimerized) -- KILLED for both Cl and Br.
  * G4 -- honest boundary: candidate 2 is named as UNTESTED, with a non-empty pre-registered
    check and a stated reason it is blocked (so a silent "tested and passed" claim would fail
    this gate).
  * G5 -- candidate 1 (DEFINITION OF DONE, chem-cfl): the E_lattice needed to make the F_phase(T)
    HT/LT crossing land on the measured Tc is positive and well below the sourced exfoliation-
    energy-per-formula-unit upper bound (arXiv:2412.03360) for both halides -- SURVIVES, not
    killed, recorded honestly either way.
  * G6 -- required_e_lattice is self-consistent: plugging the solved E_lattice back into the
    F_HT(T)-F_LT(T) crossing condition reproduces the measured Tc to high precision.
  * G7 -- the crossing is unique: the E_lattice needed to put the crossing at temperature T is
    non-decreasing in T across a wide range, so a fixed E_lattice cannot cross F_HT-F_LT more than
    once near the measured Tc.

PySCF/NumPy/SciPy only (no block2); `make gates` runs it in its own process.
"""
import numpy as np
import pytest

from nb3x8_gaps import NB3X8_HT_BULK, NB3X8_LT_BULK
from nb3x8_magnetometry import EXPERIMENT, MEV_PER_K, overprediction_factor, predicted_transition_K
from nb3x8_tc_attribution import (
    UNTESTED_CANDIDATES,
    candidate1_verdict,
    candidate3_verdict,
    exfoliation_energy_per_fu,
    ht_lt_t_ratio,
    max_family_same_motif_t_spread,
    mean_energy,
    required_e_lattice,
    required_t_rescale,
    same_motif_t_spread,
)
from nb3x8_thermo import entropy


@pytest.mark.parametrize("name,expected_factor", [("Nb3Cl8", 5.3), ("Nb3Br8", 2.3)])
def test_G1_eta1_reproduces_published_overprediction(name, expected_factor):
    """eta=1 (no rescaling) must give back the already-gated overprediction factor, else the
    machinery this module builds on has drifted."""
    p = NB3X8_LT_BULK[name]
    tc_pred = predicted_transition_K(**p)
    factor = tc_pred / EXPERIMENT[name]["Tc_K"]
    assert abs(factor - overprediction_factor(name)) < 1e-9
    assert abs(factor - expected_factor) / expected_factor < 0.02


@pytest.mark.parametrize("name", ["Nb3Cl8", "Nb3Br8"])
def test_G2_required_rescale_is_self_consistent(name):
    """Plugging eta back in must reproduce the measured Tc to high precision, and eta must be a
    genuine reduction (0 < eta < 1), matching the fact that the isolated dimer OVERpredicts."""
    p = NB3X8_LT_BULK[name]
    eta = required_t_rescale(name)
    assert 0.0 < eta < 1.0
    tc_check = predicted_transition_K(p["U0"], eta * p["t"], p["Us"])
    tc_measured = EXPERIMENT[name]["Tc_K"]
    assert abs(tc_check - tc_measured) / tc_measured < 1e-6


def test_G3_same_motif_anchor_matches_vendored_table():
    """Regression pin on the empirical anchor (bulk vs bilayer |t_s_perp|, both LT/dimerized) --
    pure arithmetic on NB3X8_CLUSTERS, already pinned by test_nb3x8_gaps_spec.py's own table."""
    expected = {"Nb3F8": 2.041, "Nb3Cl8": 0.147, "Nb3Br8": 0.118, "Nb3I8": 0.092}
    for name, pct in expected.items():
        assert abs(same_motif_t_spread(name) * 100 - pct) < 0.01, name
    assert abs(max_family_same_motif_t_spread() * 100 - 2.041) < 0.01


@pytest.mark.parametrize("name,min_needed_pct", [("Nb3Cl8", 50.0), ("Nb3Br8", 30.0)])
def test_G4_candidate3_killed_for_both_halides(name, min_needed_pct):
    """DEFINITION OF DONE: the t_s_perp rescaling needed to match measured Tc (57.7% Cl / 35.8%
    Br) is 10x+ the 2.04% family-wide same-motif anchor -- candidate 3 (a same-structural-motif
    renormalization of t_s_perp alone) is KILLED for both halides, not just asserted."""
    v = candidate3_verdict(name, margin=10.0)
    assert v["needed_pct"] > min_needed_pct
    assert v["bound_pct"] == pytest.approx(10.0 * max_family_same_motif_t_spread() * 100, rel=1e-9)
    assert v["needed_pct"] > v["bound_pct"]
    assert v["killed"] is True
    # context sanity: a genuine PHASE change (HT/LT) moves t_s_perp far more than this same-motif
    # check allows, so killing candidate 3 does not also kill candidate 1 (a distinct mechanism).
    assert 0.0 < ht_lt_t_ratio(name) < 1.0
    assert abs(1.0 - ht_lt_t_ratio(name)) > v["needed_pct"] / 100.0


def test_G4_candidate_2_is_honestly_untested():
    """Candidate 2 (in-plane kagome exchange) must be named as untested, with a concrete
    pre-registered check and a stated reason it is blocked -- guards against this module (or a
    future edit) silently claiming it was tested without the data existing. Candidate 1 has moved
    OUT of this dict (chem-cfl tested it, see test_G5_* below) -- a regression that left it in here
    too would mean the module claims it both tested and untested at once."""
    assert set(UNTESTED_CANDIDATES) == {"in_plane_kagome_exchange"}
    for key, v in UNTESTED_CANDIDATES.items():
        assert v["description"].strip(), key
        assert v["pre_registered_check"].strip(), key
        assert v["blocked_on"].strip(), key
        assert v["follow_up_bead"].strip(), key
        assert "KILL" in v["pre_registered_check"], key


@pytest.mark.parametrize("name,max_frac", [("Nb3Cl8", 0.1), ("Nb3Br8", 0.1)])
def test_G5_candidate1_survives_for_both_halides(name, max_frac):
    """DEFINITION OF DONE (chem-cfl): the lattice energy needed to pull the F_phase(T) HT/LT
    crossing onto the measured Tc is positive (correct sign -- LT must be net lattice-favored at
    T=0, matching the observed LT ground state) and a small fraction (<10%) of the physically
    plausible upper bound sourced from arXiv:2412.03360 (the full monolayer exfoliation energy per
    formula unit) -- candidate 1 SURVIVES this check for both halides, comfortably (not by tuning
    right up against the bound)."""
    v = candidate1_verdict(name)
    assert v["needed_meV"] > 0.0
    assert v["bound_meV"] == pytest.approx(exfoliation_energy_per_fu(name), rel=1e-12)
    assert v["needed_frac_of_bound"] == pytest.approx(v["needed_meV"] / v["bound_meV"], rel=1e-12)
    assert 0.0 < v["needed_frac_of_bound"] < max_frac
    assert v["killed"] is False
    # robustness: even a 2x-too-generous bound (1-vs-2-formula-units-per-cell ambiguity, see
    # exfoliation_energy_per_fu's docstring) does not flip the verdict.
    assert v["needed_meV"] < 0.5 * v["bound_meV"]


@pytest.mark.parametrize("name", ["Nb3Cl8", "Nb3Br8"])
def test_G6_required_e_lattice_is_self_consistent(name):
    """Plugging the solved E_lattice back into F_HT(Tc)-F_LT(Tc) must vanish at the measured Tc to
    high precision -- required_e_lattice solves linear algebra, not a numerical root-find, so this
    should be exact up to float round-off."""
    Tc_meV = EXPERIMENT[name]["Tc_K"] * MEV_PER_K
    ht, lt = NB3X8_HT_BULK[name], NB3X8_LT_BULK[name]
    e_lat = required_e_lattice(name)
    dE = mean_energy(**ht, T=Tc_meV) - mean_energy(**lt, T=Tc_meV)
    dS = entropy(**ht, T=Tc_meV) - entropy(**lt, T=Tc_meV)
    residual = e_lat + dE - Tc_meV * dS
    assert abs(residual) < 1e-8


@pytest.mark.parametrize("name", ["Nb3Cl8", "Nb3Br8"])
def test_G7_required_e_lattice_is_monotone_in_target_temperature(name):
    """Sanity check behind the "unique crossing" claim in the module docstring: the E_lattice
    needed to put the crossing at temperature T (T*dS_mag(T) - dE_mag(T)) must be non-decreasing in
    T across a wide range spanning the measured Tc -- otherwise a fixed E_lattice could cross F_HT-
    F_LT more than once and "the" crossing temperature would be ambiguous."""
    ht, lt = NB3X8_HT_BULK[name], NB3X8_LT_BULK[name]
    T = np.exp(np.linspace(np.log(0.01), np.log(3000.0), 2000))
    dE = mean_energy(**ht, T=T) - mean_energy(**lt, T=T)
    dS = entropy(**ht, T=T) - entropy(**lt, T=T)
    needed = T * dS - dE
    assert np.all(np.diff(needed) >= -1e-9)
