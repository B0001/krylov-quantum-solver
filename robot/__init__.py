"""CADE robot-actuation leg (Milestone D). See ROADMAP_AUTONOMOUS_SCIENCE.md Phase 3.

`robot.workspace` (D2) is not re-exported here: it imports `certkit`, which is an optional
dependency (`uv sync --extra certkit`), and `robot.planner` (D1) must stay importable without
it. Import `robot.workspace` directly once certkit is installed.
"""

from robot.planner import PlanRefusal, Trajectory, TwoSidedClaim, plan_pipette_trajectory

__all__ = ["PlanRefusal", "Trajectory", "TwoSidedClaim", "plan_pipette_trajectory"]
