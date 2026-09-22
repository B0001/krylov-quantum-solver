# SPEC: Robot planner D1 — wire certabstain's TwoSidedClaim into a mock pipette planner

**Status:** IMPLEMENTED

---

## 1. Goal

A mock wet-lab pipette trajectory planner must consume a certified two-sided claim (the
`TwoSidedClaim` contract from `certabstain`, the sibling actuation-safety repo named in
`ROADMAP_AUTONOMOUS_SCIENCE.md` Phase 3) and refuse to emit a trajectory whenever that claim
abstains, carrying the abstention reason forward instead of a plan. This is CADE Milestone D1.

## 2. Background and honest framing

- `certabstain` is not vendored into this repository (unlike `certkit`, which is pulled in as a
  pinned git dependency under `[project.optional-dependencies]`). It is a sibling repo, not
  installable here, so there is nothing to import.
- **What this spec claims:** the *interface and failure semantics* — a planner function that
  takes a two-sided claim (bracket-or-abstain) and either refuses (with reason) or emits a
  trajectory — behaves correctly on both branches.
- **What this spec does not claim:** this is not a real robot driver, not a real
  `certabstain.TwoSidedClaim` (no such import exists here), and the "trajectory" is a straight-line
  mock path with no kinematics, collision checking, or hardware binding. `TwoSidedClaim` here is a
  local, minimal reproduction of the documented contract (a certified `[lower, upper]` bracket on
  model error, or an abstention with a reason) — good enough to test the refusal semantics that
  D2/D3 build on, not a substitute for the real certabstain package.
- This is a systems/control-logic capability, not a numeric physics claim, so it has no FCI/DMRG
  analogue. The "ground truth" `specs/README.md` asks for is the contract itself: refusal must be
  total and reasoned on the abstain branch, and a trajectory must exist and stay tagged with the
  claim's certified bound on the accept branch. Both are directly checkable, which is what makes
  this falsifiable rather than ad hoc.

## 3. Approach

`robot/planner.py` defines `TwoSidedClaim` (validates its own invariants in `__post_init__`:
abstained claims carry a reason and no bracket, non-abstaining claims carry a valid `lower <=
upper` bracket) and `plan_pipette_trajectory(claim, origin, target)`, which returns either a
`Trajectory` (waypoints + `max_model_error` taken from `claim.upper`) or a `PlanRefusal` (the
claim's reason), with no third outcome.

## 4. Public interface

```
robot.planner.TwoSidedClaim(lower, upper, abstained, reason=None)
robot.planner.Trajectory(waypoints: tuple[Point, ...], max_model_error: float)
robot.planner.PlanRefusal(reason: str)
robot.planner.plan_pipette_trajectory(claim, origin, target, n_waypoints=5) -> Trajectory | PlanRefusal
```

## 5. Acceptance criteria (validation gates)

- **G1 — Abstain yields no trajectory.** A claim with `abstained=True` and a reason makes
  `plan_pipette_trajectory` return a `PlanRefusal` carrying that exact reason, never a
  `Trajectory`.
- **G2 — Valid bracket yields a trajectory.** A claim with `abstained=False` and a valid
  `lower <= upper` bracket makes `plan_pipette_trajectory` return a `Trajectory` whose first/last
  waypoints equal `origin`/`target` and whose `max_model_error == claim.upper`.
- **G3 — Claim invariants are enforced at construction.** Constructing a `TwoSidedClaim` that
  mixes abstention with a bracket, or omits the bracket while not abstaining, or omits the reason
  while abstaining, or has `lower > upper`, raises `ValueError` before it ever reaches the planner.

> Definition of done: G1-G3 green in `tests/test_robot_planner_d1_spec.py`.

## 6. Implementation plan (test-first)

1. `tests/test_robot_planner_d1_spec.py` encoding G1-G3 (written first, failing).
2. `robot/planner.py` — minimum code to pass.
3. `make gates` (this spec's test file needs no `pyscf`/`block2`, so process isolation is
   incidental, not load-bearing, here).

## 7. Out of scope

- Real `certabstain` integration (no such package is available to import).
- D2 (simulated wet-lab workspace boundary checks) and D3 (fail-closed on obstacle/drift) —
  separate beads issues, layered on this same `TwoSidedClaim`/`plan_pipette_trajectory` contract.
- Kinematics, collision avoidance, hardware I/O, timing.

## 8. Caveats and risks

- **R1:** if/when the real `certabstain` package becomes available, its actual `TwoSidedClaim`
  may not match this local reproduction field-for-field. Re-verify against the real contract
  before D2/D3 build further on it, and before any of this touches a physical actuator.
- This mock never abstains or brackets anything itself — it only reacts to a claim handed to it.
  A future step must show *something* (a real certabstain instance) actually producing claims.

## 9. Deliverables

- `robot/planner.py` — `TwoSidedClaim`, `Trajectory`, `PlanRefusal`, `plan_pipette_trajectory`.
- `tests/test_robot_planner_d1_spec.py` — gates G1-G3.
