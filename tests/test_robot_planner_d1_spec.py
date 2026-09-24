"""Spec gate for CADE Milestone D1 (specs/SPEC_robot_planner_d1.md)."""

import pytest

from robot.planner import PlanRefusal, Trajectory, TwoSidedClaim, plan_pipette_trajectory

ORIGIN = (0.0, 0.0, 0.0)
TARGET = (1.0, 2.0, 3.0)


# --- G1: abstain yields no trajectory -----------------------------------------------------


def test_abstaining_claim_yields_refusal_not_trajectory():
    claim = TwoSidedClaim(lower=None, upper=None, abstained=True, reason="drift exceeds model bound")
    result = plan_pipette_trajectory(claim, ORIGIN, TARGET)
    assert isinstance(result, PlanRefusal)
    assert result.reason == "drift exceeds model bound"
    assert not isinstance(result, Trajectory)


# --- G2: valid bracket yields a trajectory ------------------------------------------------


def test_valid_bracket_yields_trajectory_tagged_with_upper_bound():
    claim = TwoSidedClaim(lower=0.001, upper=0.02, abstained=False)
    result = plan_pipette_trajectory(claim, ORIGIN, TARGET)
    assert isinstance(result, Trajectory)
    assert result.waypoints[0] == ORIGIN
    assert result.waypoints[-1] == TARGET
    assert result.max_model_error == claim.upper


def test_trajectory_waypoint_count_is_respected():
    claim = TwoSidedClaim(lower=0.0, upper=0.05, abstained=False)
    result = plan_pipette_trajectory(claim, ORIGIN, TARGET, n_waypoints=9)
    assert isinstance(result, Trajectory)
    assert len(result.waypoints) == 9


# --- G3: claim invariants enforced at construction ----------------------------------------


def test_abstained_claim_with_bracket_is_rejected():
    with pytest.raises(ValueError):
        TwoSidedClaim(lower=0.0, upper=1.0, abstained=True, reason="conflicting")


def test_abstained_claim_without_reason_is_rejected():
    with pytest.raises(ValueError):
        TwoSidedClaim(lower=None, upper=None, abstained=True)


def test_non_abstaining_claim_without_bracket_is_rejected():
    with pytest.raises(ValueError):
        TwoSidedClaim(lower=None, upper=None, abstained=False)


def test_non_abstaining_claim_with_inverted_bracket_is_rejected():
    with pytest.raises(ValueError):
        TwoSidedClaim(lower=1.0, upper=0.0, abstained=False)


def test_n_waypoints_below_two_is_rejected():
    claim = TwoSidedClaim(lower=0.0, upper=1.0, abstained=False)
    with pytest.raises(ValueError):
        plan_pipette_trajectory(claim, ORIGIN, TARGET, n_waypoints=1)
