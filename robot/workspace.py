"""Interval-verified boundary checks for planned trajectories against a workspace envelope.

CADE Milestone D2 (see specs/SPEC_robot_workspace_d2.md). Reuses certkit.interval.Iv -- the
same outward-rounded interval primitive certkit_bridge.py uses for solver certificates --
rather than inventing a second, unaudited interval type.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from certkit.interval import Iv

from robot.planner import Point, Trajectory


class Verdict(Enum):
    INSIDE = "inside"
    OUTSIDE = "outside"
    ABSTAIN = "abstain"


_PRIORITY = {Verdict.INSIDE: 0, Verdict.ABSTAIN: 1, Verdict.OUTSIDE: 2}


def _worse(a: Verdict, b: Verdict) -> Verdict:
    return a if _PRIORITY[a] >= _PRIORITY[b] else b


@dataclass(frozen=True)
class WorkspaceEnvelope:
    """An axis-aligned box workspace envelope, one outward-rounded interval per axis."""

    x: Iv
    y: Iv
    z: Iv


def _axis_verdict(axis: Iv, pos: Iv) -> Verdict:
    if axis.encloses(pos):
        return Verdict.INSIDE
    if pos.hi < axis.lo or pos.lo > axis.hi:
        return Verdict.OUTSIDE
    return Verdict.ABSTAIN


def classify_point(p: Point, margin: float, envelope: WorkspaceEnvelope) -> Verdict:
    """Classify a nominal point, inflated by +-margin per axis, against the envelope.

    INSIDE is only returned when the entire margin-inflated interval fits inside the
    envelope on every axis; a definite exit on any axis makes the whole point OUTSIDE,
    and a boundary straddle (with no definite exit) makes it ABSTAIN.
    """
    if margin < 0:
        raise ValueError(f"margin must be >= 0, got {margin}")
    pad = Iv(-margin, margin) if margin > 0.0 else Iv.exact(0.0)
    verdict = Verdict.INSIDE
    for axis, coord in zip((envelope.x, envelope.y, envelope.z), p):
        pos = Iv.exact(coord) + pad
        verdict = _worse(verdict, _axis_verdict(axis, pos))
    return verdict


def classify_trajectory(trajectory: Trajectory, envelope: WorkspaceEnvelope) -> Verdict:
    """Worst-case verdict over every waypoint: OUTSIDE beats ABSTAIN beats INSIDE."""
    verdict = Verdict.INSIDE
    for p in trajectory.waypoints:
        verdict = _worse(verdict, classify_point(p, trajectory.max_model_error, envelope))
    return verdict
