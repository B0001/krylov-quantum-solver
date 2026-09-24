"""Mock wet-lab pipette trajectory planner gated by certabstain's TwoSidedClaim contract.

CADE Milestone D1 (see specs/SPEC_robot_planner_d1.md and ROADMAP_AUTONOMOUS_SCIENCE.md
Phase 3). `certabstain` is a sibling repo, not vendored here -- `TwoSidedClaim` below is a
local, minimal reproduction of its documented contract (a certified [lower, upper] bracket
on model error, or an abstention with a reason), good enough to test refusal semantics.
This is not a real robot driver: no kinematics, no collision checking, no hardware I/O.
"""
from __future__ import annotations

from dataclasses import dataclass

Point = tuple[float, float, float]


@dataclass(frozen=True)
class TwoSidedClaim:
    """A certified two-sided interval claim: a bracket, or an abstention with a reason."""

    lower: float | None
    upper: float | None
    abstained: bool
    reason: str | None = None

    def __post_init__(self) -> None:
        if self.abstained:
            if not self.reason:
                raise ValueError("an abstained claim must carry a reason")
            if self.lower is not None or self.upper is not None:
                raise ValueError("an abstained claim must not carry a bracket")
        else:
            if self.lower is None or self.upper is None:
                raise ValueError("a non-abstaining claim must carry a bracket")
            if self.lower > self.upper:
                raise ValueError(f"invalid bracket: lower {self.lower} > upper {self.upper}")


@dataclass(frozen=True)
class Trajectory:
    """A straight-line mock pipette path, tagged with the claim's certified model error."""

    waypoints: tuple[Point, ...]
    max_model_error: float


@dataclass(frozen=True)
class PlanRefusal:
    """What the planner returns instead of a trajectory when the claim abstains."""

    reason: str


def plan_pipette_trajectory(
    claim: TwoSidedClaim,
    origin: Point,
    target: Point,
    n_waypoints: int = 5,
) -> Trajectory | PlanRefusal:
    """Consume a TwoSidedClaim; refuse on abstention, else emit a straight-line trajectory.

    The caller (a real actuator) is responsible for keeping physical excursions within
    the returned trajectory's `max_model_error`.
    """
    if claim.abstained:
        return PlanRefusal(reason=claim.reason)

    if n_waypoints < 2:
        raise ValueError("n_waypoints must be >= 2 to include both endpoints")

    waypoints = tuple(
        tuple(o + (t - o) * i / (n_waypoints - 1) for o, t in zip(origin, target))
        for i in range(n_waypoints)
    )
    return Trajectory(waypoints=waypoints, max_model_error=claim.upper)
