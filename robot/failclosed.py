"""Fail-closed safety monitor for injected obstacle/drift events.

CADE Milestone D3 (see specs/SPEC_robot_failclosed_d3.md). Halts rather than degrades: any
event that pushes model error past SAFETY_FLOOR, or puts an obstacle within reach of the
margin-inflated path, or fails D2's workspace check, stops the monitor from continuing.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from enum import Enum

from robot.planner import Point, Trajectory, TwoSidedClaim, plan_pipette_trajectory
from robot.workspace import Verdict, WorkspaceEnvelope, classify_trajectory

# Calibrated against this spec's own fixtures, not derived from real pipette tolerances.
SAFETY_FLOOR = 0.15


class MonitorAction(Enum):
    CONTINUE = "continue"
    HALT = "halt"


@dataclass(frozen=True)
class MonitorVerdict:
    action: MonitorAction
    reason: str | None = None


@dataclass(frozen=True)
class ObstacleEvent:
    position: Point
    radius: float


@dataclass(frozen=True)
class DriftEvent:
    extra_model_error: float


def _distance(a: Point, b: Point) -> float:
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def monitor_trajectory(
    trajectory: Trajectory,
    envelope: WorkspaceEnvelope,
    obstacles: tuple[ObstacleEvent, ...] = (),
    drift: DriftEvent | None = None,
) -> MonitorVerdict:
    """Halt on floor breach, obstacle intrusion, or a non-INSIDE workspace verdict."""
    effective_margin = trajectory.max_model_error + (drift.extra_model_error if drift else 0.0)
    if effective_margin > SAFETY_FLOOR:
        return MonitorVerdict(
            MonitorAction.HALT,
            reason=(
                f"effective model error {effective_margin:.4f} exceeds "
                f"safety floor {SAFETY_FLOOR}"
            ),
        )

    for obstacle in obstacles:
        for p in trajectory.waypoints:
            if _distance(p, obstacle.position) <= obstacle.radius + effective_margin:
                return MonitorVerdict(
                    MonitorAction.HALT,
                    reason=f"obstacle at {obstacle.position} within reach of waypoint {p}",
                )

    effective_trajectory = replace(trajectory, max_model_error=effective_margin)
    verdict = classify_trajectory(effective_trajectory, envelope)
    if verdict != Verdict.INSIDE:
        return MonitorVerdict(MonitorAction.HALT, reason=f"workspace verdict {verdict.value}")

    return MonitorVerdict(MonitorAction.CONTINUE)


@dataclass(frozen=True)
class FailClosedReport:
    events: int
    halts: int
    misses: int
    sample_size: int


def run_failclosed_sweep(envelope: WorkspaceEnvelope) -> FailClosedReport:
    """Inject a fixed battery of obstacle and drift events; every one must halt.

    A ``misses == 0`` report bounds only this sweep (see spec caveats), not every possible
    obstacle/drift parameterization -- it is restated in ``sample_size`` on purpose.
    """
    origin = (0.0, 0.0, 0.0)
    target = (0.5, 0.5, 0.5)
    baseline_claim = TwoSidedClaim(lower=0.0, upper=0.05, abstained=False)
    baseline = plan_pipette_trajectory(baseline_claim, origin, target, n_waypoints=11)
    midpoint = baseline.waypoints[5]

    obstacle_events = [
        ObstacleEvent(position=midpoint, radius=r) for r in (0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
    ]
    drift_events = [
        DriftEvent(extra_model_error=e) for e in (0.15, 0.25, 0.4, 0.6, 1.0, 2.0)
    ]

    halts = 0
    misses = 0
    for obstacle in obstacle_events:
        verdict = monitor_trajectory(baseline, envelope, obstacles=(obstacle,))
        halts += verdict.action == MonitorAction.HALT
        misses += verdict.action != MonitorAction.HALT
    for drift in drift_events:
        verdict = monitor_trajectory(baseline, envelope, drift=drift)
        halts += verdict.action == MonitorAction.HALT
        misses += verdict.action != MonitorAction.HALT

    events = len(obstacle_events) + len(drift_events)
    return FailClosedReport(events=events, halts=halts, misses=misses, sample_size=events)
