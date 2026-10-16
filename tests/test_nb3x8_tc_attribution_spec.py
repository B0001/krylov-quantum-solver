"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_tc_attribution.md (bead chem-phh).

chem-g78 killed the "coordination/mean-field reduction" attribution for the Nb3X8 5.3x (Cl) /
2.3x (Br) Tc overprediction (tests/test_nb3x8_gaps_spec.py::test_G8_*), leaving the attribution
open with three named candidates (SPEC_nb3x8_magnetometry.md section 7): the cooperative/first-
order structural transition, in-plane kagome exchange, and a direct lattice renormalization of
t_s_perp. This module pre-registers a killable check for each and executes the one that is
computable from data already vendored in this repo (candidate 3).

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
  * G4 -- honest boundary: candidates 1 and 2 are named as UNTESTED, each with a non-empty
    pre-registered check and a stated reason they are blocked (so a silent "tested and passed"
    claim for either would fail this gate).

PySCF/NumPy/SciPy only (no block2); `make gates` runs it in its own process.
"""
import pytest

from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_magnetometry import EXPERIMENT, overprediction_factor, predicted_transition_K
from nb3x8_tc_attribution import (
    UNTESTED_CANDIDATES,
    candidate3_verdict,
    ht_lt_t_ratio,
    max_family_same_motif_t_spread,
    required_t_rescale,
    same_motif_t_spread,
)


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


def test_G5_candidates_1_and_2_are_honestly_untested():
    """Candidates 1 (cooperative/first-order structural transition) and 2 (in-plane kagome
    exchange) must be named as untested, each with a concrete pre-registered check and a stated
    reason it is blocked -- guards against this module (or a future edit) silently claiming they
    were tested without the data existing."""
    assert set(UNTESTED_CANDIDATES) == {
        "cooperative_structural_transition",
        "in_plane_kagome_exchange",
    }
    for key, v in UNTESTED_CANDIDATES.items():
        assert v["description"].strip(), key
        assert v["pre_registered_check"].strip(), key
        assert v["blocked_on"].strip(), key
        assert v["follow_up_bead"].strip(), key
        assert "KILL" in v["pre_registered_check"], key
