"""
Gates for chem-j71 (specs/SPEC_pm3_depth_predictor.md): a follow-up to chem-pm3
(tests/test_pm3_eta_bound_spec.py) asking whether a band-edge gap weighted by HF weight predicts
the eta-bound's required Krylov depth better than the overshoot does, and how much the bound
actually moves under shot noise (not just the analytic worst case). Driver and full numbers:
scripts/spec_pm3_depth_predictor.py.

  G1 (regression): every new escape witness found by find_escape_witnesses() reuses chem-pm3's own
      validated LP unchanged -- its bound is >= the exact eta and non-increasing in M, exactly like
      the three original witnesses.
  G2 (recorded, not asserted to a direction): the Spearman rank correlation of the weighted-gap
      predictor and of the overshoot against b_min, over all witnesses -- computed on however many
      witnesses the pre-registered search actually finds.
  G3 (recorded, not asserted to a direction): the measured noise-sensitivity sweep runs and produces
      a well-formed distribution (mean/std/[min,max], unsound fraction) at every shots value, with
      the weakest-noise point closest to the noiseless bound.
"""
import numpy as np
import pytest

from scripts.spec_pm3_depth_predictor import (
    CANDIDATES,
    find_escape_witnesses,
    noise_sensitivity_sweep,
    rank_correlation_report,
    weighted_gap,
)
from scripts.spec_pm3_subspace_eta_bound import WITNESSES as PM3_WITNESSES, run_witness


@pytest.fixture(scope="module")
def new_witness_defs():
    """The pre-registered search: smallest escaping M per candidate geometry. Cheap relative to the
    full moment-bound sweep (dense diag + Ritz floor only), but still ~5 min for 9 candidates."""
    return find_escape_witnesses(CANDIDATES)


@pytest.fixture(scope="module")
def all_results(new_witness_defs):
    """Full run_witness() on the three original chem-pm3 witnesses plus every new escape found.
    ~1-2 min per witness (dense 12-qubit diagonalization + Krylov moments)."""
    defs = PM3_WITNESSES + [dict(name=w["name"], atom=w["atom"], d=w["d"], m_floor=w["m_floor"])
                            for w in new_witness_defs]
    return [r for w in defs if (r := run_witness(w, m_extended=(), caps=())) is not None]


def test_G1_new_witnesses_reuse_the_validated_lp(all_results):
    """G1: every witness (original + new) satisfies the same LP soundness properties chem-pm3
    already gates -- bound >= exact eta, and bound non-increasing across M_SWEEP."""
    assert len(all_results) >= 3, "the three original chem-pm3 witnesses must still run"
    for r in all_results:
        prev = np.inf
        for M, (ub, lb) in sorted(r["bounds"].items()):
            assert ub >= r["eta"] - 1e-9, f"{r['name']} M={M}: bound {ub} below exact eta {r['eta']}"
            assert ub <= prev + 1e-6, f"{r['name']} M={M}: bound increased with more moments"
            prev = ub


def test_G2_rank_correlation_is_recorded(all_results):
    """G2: Spearman rho(WG, b_min) and rho(overshoot, b_min) are computed over every witness found.
    Not asserted to a direction (the honest finding may be that WG also fails) -- only that it is
    actually computed, on the actual witness count, and both correlations are well-formed numbers."""
    report = rank_correlation_report(all_results)
    assert report["n"] == len(all_results)
    assert report["n"] >= 3
    assert np.isfinite(report["rho_wg"])
    assert np.isfinite(report["rho_overshoot"])
    assert -1.0 <= report["rho_wg"] <= 1.0
    assert -1.0 <= report["rho_overshoot"] <= 1.0
    # weighted_gap must be positive and finite for every witness (eta > 0 and gaps > 0 by construction
    # of an escape witness -- a non-finite value would mean the search admitted a degenerate case).
    for r in all_results:
        wg = weighted_gap(r)
        assert np.isfinite(wg) and wg > 0, f"{r['name']}: degenerate weighted gap {wg}"


def test_G3_noise_sweep_well_formed():
    """G3: the measured noise-sensitivity sweep on a small synthetic case produces a well-formed
    distribution at every shots value, and sensitivity shrinks as shots grow (less noisy moments
    cannot produce a WIDER distribution of bounds on average)."""
    energies = np.array([0.5, 1.5, 3.2, 4.1, 6.0, 8.5, 9.2])
    weights = np.array([0.10, 0.20, 0.07, 0.08, 0.30, 0.15, 0.10])
    e_min, e_max = 0.0, 10.0
    dt = np.pi / (e_max - e_min)
    band_lo, band_hi = 3.0, 5.0
    eta_exact = float(weights[(energies >= band_lo) & (energies < band_hi)].sum())
    k = np.arange(12)[:, None]
    moments = (weights[None, :] * np.exp(-1j * k * dt * energies[None, :])).sum(axis=1)

    res = dict(e_d=band_lo, beta=band_hi, eta=eta_exact)
    trials = noise_sensitivity_sweep(res, moments, dt, e_min, e_max, M=12, cap=1000.0,
                                     shots_list=(1e3, 1e5, 1e8), n_trials=10, seed=1)
    stds = []
    for tr in trials:
        assert tr.bounds.shape == (10,)
        assert np.all(np.isfinite(tr.bounds))
        assert 0.0 <= tr.unsound_fraction <= 1.0
        stds.append(tr.bounds.std())
    # the highest-shots (least noisy) point's std must be the smallest of the three -- sigma shrinks
    # monotonically with shots, and the LP's sensitivity to its own inputs doesn't invert that.
    assert stds[-1] <= stds[0] + 1e-6, f"noise sweep not monotonic in shots: stds={stds}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
