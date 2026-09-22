"""Spec gate for CADE Milestone D2 (specs/SPEC_robot_workspace_d2.md)."""

import pytest
from certkit.interval import Iv

from robot.planner import TwoSidedClaim, plan_pipette_trajectory
from robot.workspace import (
    Verdict,
    WorkspaceEnvelope,
    classify_point,
    classify_trajectory,
)

ENVELOPE = WorkspaceEnvelope(x=Iv(-1.0, 1.0), y=Iv(-1.0, 1.0), z=Iv(-1.0, 1.0))
ORIGIN = (0.0, 0.0, 0.0)
MARGIN = 0.05


def _trajectory_to(target, margin=MARGIN):
    claim = TwoSidedClaim(lower=0.0, upper=margin, abstained=False)
    return plan_pipette_trajectory(claim, ORIGIN, target, n_waypoints=11)


# --- G1: fully-enclosed trajectory is INSIDE ----------------------------------------------


def test_fully_enclosed_trajectory_is_inside():
    traj = _trajectory_to((0.5, 0.5, 0.5))
    assert classify_trajectory(traj, ENVELOPE) == Verdict.INSIDE


# --- G2 + G4: zero false-safe across a sweep, all three verdicts appear -------------------


def _sweep():
    # target x sweeps from well inside to well outside the [-1, 1] envelope; y, z fixed at 0.
    xs = [-0.5 + 0.1 * i for i in range(26)]  # -0.5 .. 2.0
    return [(x, _trajectory_to((x, 0.0, 0.0))) for x in xs]


def test_zero_false_safe_across_sweep():
    for x, traj in _sweep():
        verdict = classify_trajectory(traj, ENVELOPE)
        # ground truth: does the true (margin-zero) nominal path leave the envelope?
        true_exits = any(abs(x * i / 10) > 1.0 for i in range(11))
        if true_exits:
            assert verdict != Verdict.INSIDE, f"false-safe at x={x}: {verdict}"


def test_sweep_exercises_all_three_verdicts():
    verdicts = {classify_trajectory(traj, ENVELOPE) for _, traj in _sweep()}
    assert verdicts == {Verdict.INSIDE, Verdict.OUTSIDE, Verdict.ABSTAIN}


# --- G3: straddle is ABSTAIN, not INSIDE ---------------------------------------------------


def test_straddling_waypoint_is_abstain_not_inside():
    # nominal point sits exactly on the boundary; margin makes the interval straddle it.
    traj = _trajectory_to((1.0, 0.0, 0.0))
    assert classify_trajectory(traj, ENVELOPE) == Verdict.ABSTAIN


def test_classify_point_straddle_directly():
    verdict = classify_point((1.0, 0.0, 0.0), MARGIN, ENVELOPE)
    assert verdict == Verdict.ABSTAIN


def test_classify_point_definitely_outside():
    verdict = classify_point((1.5, 0.0, 0.0), MARGIN, ENVELOPE)
    assert verdict == Verdict.OUTSIDE


def test_classify_point_definitely_inside():
    verdict = classify_point((0.0, 0.0, 0.0), MARGIN, ENVELOPE)
    assert verdict == Verdict.INSIDE


# --- G5: negative margin rejected -----------------------------------------------------------


def test_negative_margin_rejected():
    with pytest.raises(ValueError):
        classify_point((0.0, 0.0, 0.0), -0.1, ENVELOPE)
