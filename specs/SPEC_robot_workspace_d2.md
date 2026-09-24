# SPEC: Robot workspace D2 — interval-verified boundary checks on planned trajectories

**Status:** IMPLEMENTED

---

## 1. Goal

Classify every planned pipette motion against a simulated workspace envelope as INSIDE,
OUTSIDE, or ABSTAIN, using outward-rounded interval arithmetic, such that **no trajectory that
actually exits the envelope is ever classified INSIDE** across a sweep of motions that stay
inside, straddle the boundary, and exit it. This is CADE Milestone D2, built directly on D1's
`Trajectory`/`TwoSidedClaim` contract (`specs/SPEC_robot_planner_d1.md`).

## 2. Background and honest framing

- **Reused, not reinvented:** the outward-rounded interval type is `certkit.interval.Iv` — the
  same producer-independent primitive `certkit_bridge.py` uses to state solver certificates. Its
  soundness contract (every op returns an interval that *provably* contains the exact result,
  via `nextafter` widening) is exactly what "outward-rounded" in the acceptance criteria means; a
  bespoke interval type here would just be a second, unaudited implementation of the same idea.
- **What this spec claims:** given a trajectory's waypoints and its certified `max_model_error`
  (from D1's `TwoSidedClaim.upper`), inflating each nominal waypoint into a `[x-margin, x+margin]`
  interval per axis and testing envelope containment on that interval — rather than on the bare
  nominal point — never calls a truly-exiting motion safe.
- **What this spec does not claim:** the envelope is an axis-aligned box (no arbitrary geometry),
  "workspace" and "obstacle" are simulated coordinates only (no hardware, no real sensor drift —
  that is D3's `specs/SPEC_robot_failclosed_d3.md`), and ABSTAIN here is a boundary-straddle
  signal, not yet wired to any halting actuator.

## 3. Approach

For each waypoint `p = (x, y, z)` and trajectory margin `r = max_model_error`, build a per-axis
position interval `Iv.exact(coord) + Iv(-r, r)` (outward-rounded by `Iv.__add__`) and compare it
against the corresponding envelope axis interval:

- `axis.encloses(pos)` → **INSIDE** on that axis.
- `pos` disjoint from `axis` (`pos.hi < axis.lo or pos.lo > axis.hi`) → **OUTSIDE** on that axis.
- Otherwise (overlaps but neither encloses nor is disjoint) → **ABSTAIN** on that axis.

A trajectory's verdict is the worst case over all waypoints and axes, ranked
`OUTSIDE > ABSTAIN > INSIDE`: one provably-exiting waypoint makes the whole trajectory OUTSIDE;
absent that, one straddling waypoint makes it ABSTAIN; only a trajectory whose every
margin-inflated waypoint is fully enclosed on every axis is INSIDE.

## 4. Public interface

```
robot.workspace.Verdict            # Enum: INSIDE, OUTSIDE, ABSTAIN
robot.workspace.WorkspaceEnvelope(x: Iv, y: Iv, z: Iv)
robot.workspace.classify_point(p: Point, margin: float, envelope: WorkspaceEnvelope) -> Verdict
robot.workspace.classify_trajectory(trajectory: Trajectory, envelope: WorkspaceEnvelope) -> Verdict
```

## 5. Acceptance criteria (validation gates)

- **G1 — Fully-enclosed trajectory is INSIDE.** A trajectory whose every margin-inflated waypoint
  is enclosed by the envelope on every axis classifies INSIDE.
- **G2 — Zero false-safe across a parametrized sweep.** Sweep a trajectory's target coordinate
  from well inside the envelope to well outside it. For every step where the *true* (margin-zero)
  nominal path leaves the envelope on any axis, `classify_trajectory` never returns INSIDE.
- **G3 — Straddle is ABSTAIN, not INSIDE.** A waypoint placed so its margin-inflated interval
  overlaps the envelope boundary without the nominal point leaving it classifies ABSTAIN, not
  INSIDE — conservatism under uncertainty, not just under a definite exit.
- **G4 — Sweep exercises all three verdicts.** The G2 sweep must produce at least one INSIDE, one
  OUTSIDE, and one ABSTAIN classification (else the sweep doesn't actually cover the boundary).
- **G5 — Negative margin is rejected.** `classify_point` raises `ValueError` for `margin < 0`.

> Definition of done: G1-G5 green in `tests/test_robot_workspace_d2_spec.py`.

## 6. Implementation plan (test-first)

1. `tests/test_robot_workspace_d2_spec.py` encoding G1-G5 (written first, failing).
2. `robot/workspace.py` — minimum code to pass, built on `certkit.interval.Iv` and D1's
   `robot.planner.Trajectory`.
3. `make gates`.

## 7. Out of scope

- Non-box envelope geometry.
- D3 (fail-closed halting on injected obstacle/drift events) — a separate beads issue that
  consumes this module's ABSTAIN/OUTSIDE verdicts as its trigger.
- Any real sensor, actuator, or hardware binding.

## 8. Caveats and risks

- **R1:** the margin here is D1's `max_model_error`, a single scalar radius applied uniformly to
  every waypoint. A real system's uncertainty likely grows along the path (drift accumulates);
  this spec does not model that — see D3's caveats for where drift is actually injected.
- G2/G4's sweep is a fixed, hand-chosen axis and step size; it demonstrates the property on one
  cut through the envelope, not an exhaustive search of all possible trajectories.

## 9. Deliverables

- `robot/workspace.py` — `Verdict`, `WorkspaceEnvelope`, `classify_point`, `classify_trajectory`.
- `tests/test_robot_workspace_d2_spec.py` — gates G1-G5.
