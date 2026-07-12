"""Spec gate for ScreenLoop (specs/SPEC_screenloop.md, tasks 1-6). Mostly pure/fast."""

import subprocess
import sys

import pytest

from screenloop import (
    Interval,
    SyntheticOracle,
    Verdict,
    classify,
    point_estimate_screen,
    screen,
)

TARGET = Interval(4.0, 6.0)
# 50 candidates with fixed ground truth 0.0, 0.2, ..., 9.8. Hits (true in [4,6]) = ids 20..30 (11).
TRUTHS = {i: round(i * 0.2, 3) for i in range(50)}
TRUE_HITS = {i for i, v in TRUTHS.items() if TARGET.contains(v)}


def _oracle(**kw):
    return SyntheticOracle(TRUTHS, **kw)


def _exhaustive(candidates, target, oracle, precision):
    """Evaluate everything at full precision once; survivors = non-eliminated. Returns (set, cost)."""
    survivors, cost = set(), 0.0
    for c in candidates:
        cost += oracle.cost(precision)
        if classify(oracle.bracket(c, precision), target) is not Verdict.ELIMINATED:
            survivors.add(c)
    return survivors, cost


# --- Gate 1: the dominance rule ----------------------------------------------------------


def test_classify_hand_cases():
    assert classify(Interval(1.0, 2.0), TARGET) is Verdict.ELIMINATED   # fully below
    assert classify(Interval(7.0, 8.0), TARGET) is Verdict.ELIMINATED   # fully above
    assert classify(Interval(4.5, 5.5), TARGET) is Verdict.CONFIRMED    # inside
    assert classify(Interval(3.0, 5.0), TARGET) is Verdict.UNDECIDED    # straddles lower
    # boundary contact is conservative: touching, not disjoint -> kept (undecided), never eliminated
    assert classify(Interval(2.0, 4.0), TARGET) is Verdict.UNDECIDED
    assert classify(Interval(6.0, 9.0), TARGET) is Verdict.UNDECIDED


def test_pruning_imports_without_solver_deps():
    code = (
        "import screenloop.pruning, screenloop.screen, sys;"
        "bad=[m for m in ('pyscf','qiskit','numpy','scipy') if m in sys.modules];"
        "assert not bad, bad; print('clean')"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert "clean" in proc.stdout


# --- Gate 2: zero false eliminations (THE invariant) ------------------------------------


@pytest.mark.parametrize("acquisition", [False, True])
def test_zero_false_eliminations(acquisition):
    candidates = list(TRUTHS)
    result = screen(candidates, TARGET, _oracle(), acquisition=acquisition)
    eliminated = set(result.eliminated)
    # No genuine hit is ever eliminated — the whole point.
    assert eliminated.isdisjoint(TRUE_HITS), f"false eliminations: {eliminated & TRUE_HITS}"
    # And every true hit survives.
    assert TRUE_HITS <= set(result.survivors)


def test_screen_matches_exhaustive_hit_set():
    candidates = list(TRUTHS)
    exhaustive_hits, _ = _exhaustive(candidates, TARGET, _oracle(), precision=4)
    for acquisition in (False, True):
        result = screen(candidates, TARGET, _oracle(), max_precision=4, acquisition=acquisition)
        assert set(result.survivors) == exhaustive_hits


# --- Gate 3: acquisition/loop is cheaper than exhaustive-at-full ------------------------


def test_loop_cheaper_than_exhaustive():
    candidates = list(TRUTHS)
    _, exhaustive_cost = _exhaustive(candidates, TARGET, _oracle(), precision=4)
    result = screen(candidates, TARGET, _oracle(), max_precision=4, acquisition=True)
    assert result.total_cost < 0.5 * exhaustive_cost  # spec: < half the exhaustive cost


# --- Gate 4: the point-estimate baseline is unsound ------------------------------------


def test_baseline_makes_false_eliminations():
    candidates = list(TRUTHS)
    baseline = point_estimate_screen(candidates, TARGET, _oracle(), precision=0)
    baseline_false = set(baseline.eliminated) & TRUE_HITS
    assert len(baseline_false) > 0, "baseline should falsely eliminate some true hits"

    # Same space, same oracle: the bracket-aware loop makes zero.
    loop = screen(candidates, TARGET, _oracle())
    assert (set(loop.eliminated) & TRUE_HITS) == set()


# --- Gate 5: any BoundedOracle plugs in unchanged --------------------------------------


def test_second_oracle_plugs_in_unchanged():
    # A "conformal-ML mock": different (wider) interval widths, different bias schedule.
    ml_oracle = _oracle(w0=2.0, b0=1.2, ratio=0.5, seed=7)
    result = screen(list(TRUTHS), TARGET, ml_oracle, max_precision=5)
    assert (set(result.eliminated) & TRUE_HITS) == set()
    assert TRUE_HITS <= set(result.survivors)


# --- certchem-backed oracle smoke test -------------------------------------------------


def test_certchem_oracle_smoke():
    from screenloop import CertchemEnergyOracle

    # Three H2 bond lengths; screen for "ground-state energy in a window that only ~0.9 A hits".
    candidates = [
        ("H 0 0 0; H 0 0 0.7", "sto-3g", (2, 2)),
        ("H 0 0 0; H 0 0 0.9", "sto-3g", (2, 2)),
        ("H 0 0 0; H 0 0 1.5", "sto-3g", (2, 2)),
    ]
    oracle = CertchemEnergyOracle()
    # Whole-range target: nothing eliminated, loop runs end to end.
    result = screen(candidates, Interval(-1.2, -1.0), oracle, max_precision=1)
    assert result.n_bracket_calls >= len(candidates)
    assert len(result.survivors) + len(result.eliminated) == len(candidates)
