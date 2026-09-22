# SPEC: Robot fail-closed D3 — halt on injected obstacle/drift events

**Status:** IMPLEMENTED

---

## 1. Goal

A safety monitor sitting downstream of D1's planner and D2's workspace check must **halt**
(never silently continue on a degraded plan) whenever an injected obstacle intrusion or drift
event pushes a trajectory's effective model error past a safety floor, or the resulting
margin-inflated path is no longer certifiably INSIDE the workspace envelope. This is CADE
Milestone D3, the last child of the robot-actuation epic (`chem-cjl`).

## 2. Background and honest framing

- **What this spec claims:** over a finite, stated-size sweep of injected obstacle and drift
  events, every single one produces a HALT action from the monitor, and a nominal trajectory with
  no injected event does not halt (the monitor isn't just unconditionally halting, which would
  pass G1 vacuously). The report states `events`, `halts`, `misses`, and `sample_size` explicitly.
- **What this spec does not claim:** zero misses on a sweep of size N is not "the monitor never
  misses" — it is "0/N missed on this sweep," and the sweep's obstacle/drift parameter ranges are
  hand-chosen, not exhaustive. Per the acceptance criteria, this spec states the bound the sweep
  actually supports, not a rounded-up 100%. No real sensor is involved; "drift" and "obstacle" are
  simulated scalar/geometric events layered onto D1/D2's existing contracts.
- **Reused, not reinvented:** the monitor calls straight into D2's `classify_trajectory` /
  `WorkspaceEnvelope` for the boundary check; it does not reimplement interval arithmetic.

## 3. Approach

`ObstacleEvent(position, radius)` and `DriftEvent(extra_model_error)` are injected against a
baseline trajectory (from D1). `monitor_trajectory`:

1. Computes the effective model error = `trajectory.max_model_error + drift.extra_model_error`
   (0 if no drift event). If this exceeds `SAFETY_FLOOR`, **HALT** immediately — a monitor that
   let a known-oversized error pass through to the workspace check would be degrading, not
   fail-closed.
2. Checks every waypoint against every obstacle: if a waypoint lies within
   `obstacle.radius + effective_margin` of `obstacle.position`, **HALT** (the margin-inflated
   pipette tip could physically reach the obstacle).
3. Otherwise, rebuilds the trajectory with the effective margin and defers to D2's
   `classify_trajectory`; anything other than `Verdict.INSIDE` **HALTs**.
4. Only if none of the above fire does the monitor return **CONTINUE**.

`SAFETY_FLOOR` is a calibrated constant (`0.15`, in the same length units as D1/D2's margins,
i.e. roughly 3x the `0.05` margin used in the D1/D2 gates) — picked to sit above nominal
operation and below the drift magnitudes this spec's sweep injects, not derived from any
physical model.

## 4. Public interface

```
robot.failclosed.SAFETY_FLOOR: float
robot.failclosed.MonitorAction            # Enum: CONTINUE, HALT
robot.failclosed.MonitorVerdict(action: MonitorAction, reason: str | None)
robot.failclosed.ObstacleEvent(position: Point, radius: float)
robot.failclosed.DriftEvent(extra_model_error: float)
robot.failclosed.monitor_trajectory(trajectory, envelope, obstacles=(), drift=None) -> MonitorVerdict
robot.failclosed.FailClosedReport(events: int, halts: int, misses: int, sample_size: int)
robot.failclosed.run_failclosed_sweep(envelope) -> FailClosedReport
```

## 5. Acceptance criteria (validation gates)

- **G1 — Every injected obstacle event halts.** For each obstacle placed directly on the nominal
  path, `monitor_trajectory` returns `HALT`.
- **G2 — Every injected drift event past the floor halts.** For each drift event whose
  `extra_model_error` pushes the effective margin above `SAFETY_FLOOR`, `monitor_trajectory`
  returns `HALT`.
- **G3 — Nominal trajectory (no event) continues.** A baseline trajectory with no obstacle, no
  drift, and margin below `SAFETY_FLOOR`, fully inside the envelope, returns `CONTINUE` — the
  monitor is not trivially always-halt.
- **G4 — Sweep report states its own bound honestly.** `run_failclosed_sweep` returns a
  `FailClosedReport` with `misses == 0`, `sample_size == events` (every injected event is
  accounted for), and `events > 0` (a non-vacuous sweep); any nonzero `misses` fails the gate.

> Definition of done: G1-G4 green in `tests/test_robot_failclosed_d3_spec.py`.

## 6. Implementation plan (test-first)

1. `tests/test_robot_failclosed_d3_spec.py` encoding G1-G4 (written first, failing).
2. `robot/failclosed.py` — minimum code to pass, built on D1's `Trajectory` and D2's
   `classify_trajectory`/`WorkspaceEnvelope` (so it also requires the `certkit` extra).
3. `make gates`.

## 7. Out of scope

- Any real obstacle sensor, drift estimator, or actuator halt signal — `ObstacleEvent`/
  `DriftEvent` are directly-injected test fixtures, not derived from simulated physics.
- Recovery/replanning after a halt (this spec only checks that the monitor halts, not what
  happens next).
- This closes the D1→D2→D3 chain under `chem-cjl`; the epic's own acceptance criterion ("a
  verified molecular target flows... through an actuation plan, with every hop carrying a
  certificate that could have failed") is satisfied at the *interface* level by this chain, not by
  wiring in a real chem/certkit-verified molecular target end-to-end, which remains future work.

## 8. Caveats and risks

- **R1:** `SAFETY_FLOOR = 0.15` is calibrated against this spec's own test fixtures, not derived
  from any real pipette tolerance. Do not carry this constant to a real actuator without
  re-deriving it from actual hardware error bounds.
- G4's `misses == 0` is a sweep-size-N result, restated in §2: it bounds this sweep, not all
  possible obstacle/drift parameterizations.

## 9. Deliverables

- `robot/failclosed.py` — `MonitorAction`, `MonitorVerdict`, `ObstacleEvent`, `DriftEvent`,
  `monitor_trajectory`, `FailClosedReport`, `run_failclosed_sweep`.
- `tests/test_robot_failclosed_d3_spec.py` — gates G1-G4.
