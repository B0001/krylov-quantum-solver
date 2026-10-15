"""
Acceptance gates G1-G7 for specs/SPEC_nb3x8_tc_attribution.md (bead chem-phh / chem-b0y).

chem-g78 killed the "coordination/mean-field reduction" attribution for the Nb3X8 5.3x (Cl) /
2.3x (Br) Tc overprediction (tests/test_nb3x8_gaps_spec.py::test_G8_*), leaving the attribution
open with three named candidates (SPEC_nb3x8_magnetometry.md section 7): the cooperative/first-
order structural transition, in-plane kagome exchange, and a direct lattice renormalization of
t_s_perp. This module pre-registers a killable check for each; candidates 2 (chem-b0y) and 3
(chem-phh) are now executed, candidate 1 remains blocked on unvendored elastic/lattice data.

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
  * G5 -- honest boundary: candidate 1 is named as UNTESTED, with a non-empty pre-registered
    check and a stated reason it is blocked (so a silent "tested and passed" claim would fail
    this gate); candidate 2 must NOT appear there any more now that it has been executed.
  * G6 -- candidate 2 sourced data pin: NB3X8_INPLANE_EXCHANGE's J1/J2 values match the primary
    source (arXiv:2508.17352 Table 1) exactly -- a regression pin against accidental edits.
  * G7 -- candidate 2 (DEFINITION OF DONE): the RPA-required effective in-plane coupling exceeds
    the most generous sourced bound (and exceeds J_interlayer itself) for both halides -- KILLED.

PySCF/NumPy/SciPy only (no block2); `make gates` runs it in its own process.
"""
import pytest

from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_magnetometry import EXPERIMENT, overprediction_factor, predicted_transition_K
from nb3x8_tc_attribution import (
    NB3X8_INPLANE_EXCHANGE,
    UNTESTED_CANDIDATES,
    candidate2_verdict,
    candidate3_verdict,
    ht_lt_t_ratio,
    max_family_same_motif_t_spread,
    required_inplane_coupling,
    required_t_rescale,
    same_motif_t_spread,
    sourced_inplane_bound,
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


def test_G5_candidate1_is_honestly_untested_and_candidate2_is_not():
    """Candidate 1 (cooperative/first-order structural transition) must be named as untested,
    with a concrete pre-registered check and a stated reason it is blocked -- guards against this
    module (or a future edit) silently claiming it was tested without the data existing. Candidate
    2 (in-plane kagome exchange) must NOT appear here any more: chem-b0y sourced data for it and
    executed the pre-registered check (see test_G5/test_G6 below)."""
    assert set(UNTESTED_CANDIDATES) == {"cooperative_structural_transition"}
    for key, v in UNTESTED_CANDIDATES.items():
        assert v["description"].strip(), key
        assert v["pre_registered_check"].strip(), key
        assert v["blocked_on"].strip(), key
        assert v["follow_up_bead"].strip(), key
        assert "KILL" in v["pre_registered_check"], key


@pytest.mark.parametrize(
    "name,j1,j2",
    [("Nb3Cl8", 0.12, -0.24), ("Nb3Br8", -2.74, -0.37)],
)
def test_G6_candidate2_sourced_data_matches_primary_source(name, j1, j2):
    """Regression pin: NB3X8_INPLANE_EXCHANGE must match arXiv:2508.17352 Table 1 exactly (meV
    per Nb3 trimer, 'positive = AFM' convention), verified against the downloaded PDF text."""
    c = NB3X8_INPLANE_EXCHANGE[name]
    assert c["J1"] == pytest.approx(j1)
    assert c["J2"] == pytest.approx(j2)
    assert "2508.17352" in c["source"]


@pytest.mark.parametrize("name,min_needed_over_bound", [("Nb3Cl8", 1000.0), ("Nb3Br8", 10.0)])
def test_G7_candidate2_killed_for_both_halides(name, min_needed_over_bound):
    """DEFINITION OF DONE for candidate 2: the RPA-required effective in-plane coupling
    (steelmanning the conjecture with the exact self-consistency closure J'_eff*chi_dimer(Tc)=1)
    exceeds even the most generous sourced bound (full triangular coordination, |J1|+|J2|, no
    cancellation) by orders of magnitude for both halides -- and exceeds J_interlayer itself,
    which is the literal pre-registered kill condition in UNTESTED_CANDIDATES before chem-b0y
    resolved it ('requires J' unphysically close to or exceeding J_interlayer')."""
    v = candidate2_verdict(name)
    assert v["needed_meV"] == pytest.approx(required_inplane_coupling(name))
    assert v["bound_meV"] == pytest.approx(sourced_inplane_bound(name))
    assert v["needed_over_bound"] > min_needed_over_bound
    assert v["needed_over_j_interlayer"] > 1.0
    assert v["killed"] is True
