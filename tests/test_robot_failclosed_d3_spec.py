"""Spec gate for CADE Milestone D3 (specs/SPEC_robot_failclosed_d3.md)."""

import importlib.util
import os

import pytest

if importlib.util.find_spec("certkit") is None:
    if os.environ.get("CI"):
        raise RuntimeError(
            "certkit is not installed: the failclosed D3 gate cannot run. "
            'Install the pinned checker with `uv pip install -e ".[certkit]"`.'
        )
    pytest.skip('certkit extra not installed (uv pip install -e ".[certkit]")',
                allow_module_level=True)

from certkit.interval import Iv  # noqa: E402

from robot.failclosed import (  # noqa: E402
    SAFETY_FLOOR,
    DriftEvent,
    MonitorAction,
    ObstacleEvent,
    monitor_trajectory,
    run_failclosed_sweep,
)
from robot.planner import TwoSidedClaim, plan_pipette_trajectory  # noqa: E402
from robot.workspace import WorkspaceEnvelope  # noqa: E402

ENVELOPE = WorkspaceEnvelope(x=Iv(-1.0, 1.0), y=Iv(-1.0, 1.0), z=Iv(-1.0, 1.0))
ORIGIN = (0.0, 0.0, 0.0)
TARGET = (0.5, 0.5, 0.5)


def _baseline(margin=0.05):
    claim = TwoSidedClaim(lower=0.0, upper=margin, abstained=False)
    return plan_pipette_trajectory(claim, ORIGIN, TARGET, n_waypoints=11)


# --- G1: every injected obstacle event halts ------------------------------------------------


def test_obstacle_on_path_halts():
    traj = _baseline()
    midpoint = traj.waypoints[5]
    for radius in (0.01, 0.05, 0.1, 0.2, 0.3, 0.5):
        verdict = monitor_trajectory(traj, ENVELOPE, obstacles=(ObstacleEvent(midpoint, radius),))
        assert verdict.action == MonitorAction.HALT, f"radius={radius} missed"
        assert verdict.reason


# --- G2: every injected drift event past the floor halts ------------------------------------


def test_drift_past_floor_halts():
    traj = _baseline()
    for extra in (0.15, 0.25, 0.4, 0.6, 1.0, 2.0):
        assert traj.max_model_error + extra > SAFETY_FLOOR
        verdict = monitor_trajectory(traj, ENVELOPE, drift=DriftEvent(extra_model_error=extra))
        assert verdict.action == MonitorAction.HALT, f"extra={extra} missed"
        assert verdict.reason


# --- G3: nominal trajectory (no event) continues ---------------------------------------------


def test_nominal_trajectory_continues():
    traj = _baseline()
    assert traj.max_model_error < SAFETY_FLOOR
    verdict = monitor_trajectory(traj, ENVELOPE)
    assert verdict.action == MonitorAction.CONTINUE


# --- G4: sweep report states its own bound honestly -------------------------------------------


def test_failclosed_sweep_reports_zero_misses_with_stated_sample_size():
    report = run_failclosed_sweep(ENVELOPE)
    assert report.events > 0
    assert report.sample_size == report.events
    assert report.misses == 0
    assert report.halts == report.events
