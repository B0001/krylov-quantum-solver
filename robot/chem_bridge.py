"""Bridge a real, independently-checked certkit certificate into D1's actuation contract.

CADE Milestone D4 (see specs/SPEC_robot_chem_bridge_d4.md). Closes the gap D3 left open:
D1-D3 gate the planner/workspace/fail-closed chain against a local TwoSidedClaim stand-in,
never against a certificate this repo actually produces. `certkit_bridge.Verdict` is a
producer-side record -- its lo/hi fields are populated before the independent checker runs,
so they are present and finite even on an ABSTAIN. This module gates on `ok` alone, never on
lo/hi's mere presence, and treats the certified enclosure's width (not its raw bounds) as the
model-error bound D1's Trajectory expects.

Not re-exported from robot/__init__.py: this module imports certkit_bridge, which requires
the certkit extra and the full solver stack, and D1 must stay importable without either
(same reasoning robot/workspace.py's docstring already gives for the certkit extra alone).
"""
from __future__ import annotations

from certkit_bridge import Verdict
from robot.planner import (
    PlanRefusal,
    Point,
    Trajectory,
    TwoSidedClaim,
    plan_pipette_trajectory,
)


def verdict_to_claim(verdict: Verdict) -> TwoSidedClaim:
    """Map a certkit_bridge.Verdict onto D1's TwoSidedClaim contract.

    `ok=False` (ABSTAIN or a crash) always abstains, regardless of what lo/hi carry -- a
    producer-side Verdict can have a finite bracket attached to an unverified certificate.
    `ok=True` yields a non-abstaining claim whose upper bound is the certified enclosure's
    width (hi - lo): what a verified enclosure actually promises is that the true value is
    confined to an interval this wide, which is the natural stand-in for the model-error
    bound D1's Trajectory carries (see spec section 3 for why, and its caveats).
    """
    if not verdict.ok:
        return TwoSidedClaim(lower=None, upper=None, abstained=True, reason=verdict.line)
    return TwoSidedClaim(lower=0.0, upper=verdict.hi - verdict.lo, abstained=False)


def plan_from_verdict(
    verdict: Verdict,
    origin: Point,
    target: Point,
    n_waypoints: int = 5,
) -> Trajectory | PlanRefusal:
    """Plan a trajectory straight from a real certificate, via verdict_to_claim."""
    claim = verdict_to_claim(verdict)
    return plan_pipette_trajectory(claim, origin, target, n_waypoints=n_waypoints)
