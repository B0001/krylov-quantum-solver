# SPEC: Robot chem bridge D4 — wire a real certkit certificate into the actuation chain

**Status:** IMPLEMENTED

---

## 1. Goal

CADE epic `chem-cjl`'s acceptance criterion is a single flow: a verified molecular target moves
from a chem certificate through certkit re-proof into a certabstain-guarded actuation plan, with
every hop carrying a certificate that could have failed. D1-D3 built and gated that chain
(`robot.planner` -> `robot.workspace` -> `robot.failclosed`) against a **local** `TwoSidedClaim`
stand-in, never against a certificate this repo actually produces (see
`SPEC_robot_failclosed_d3.md` section 7). This spec closes that gap: a real, independently-checked
`certkit_bridge.Verdict` — the same object `certkit_bridge.py` emits for an actual molecule — is
translated into a `TwoSidedClaim` and driven, unmodified, through D1's planner, D2's workspace
check, and D3's fail-closed monitor.

## 2. Background and honest framing

- **What this spec claims:** a real `Verdict` from `certkit_bridge.run_case` — `ok=True` with a
  finite `[lo, hi]` enclosure, or `ok=False` (ABSTAIN or a crash) — maps deterministically onto
  D1's `TwoSidedClaim` contract, and the mapped claim drives the existing, unmodified D1/D2/D3
  chain to a real outcome (a trajectory, a workspace verdict, a monitor action) end-to-end, with no
  changes to `robot/planner.py`, `robot/workspace.py`, or `robot/failclosed.py`.
- **What this spec does not claim:** this does not vendor or reproduce `certabstain` — `TwoSidedClaim`
  remains D1's local stand-in for its documented contract. It does not claim the bracket-width ->
  model-error mapping below is physically derived; it is a stated modeling choice (see §3, §8 R1).
  It does not wire a real robot, sensor, or wet-lab synthesis step — "actuation" here still ends at
  a `Trajectory`/`MonitorVerdict` object, exactly as in D1-D3.
- **Reused, not reinvented:** `certkit_bridge.run_case` (unmodified) produces the certificate;
  `robot.planner.plan_pipette_trajectory`, `robot.workspace.classify_trajectory`, and
  `robot.failclosed.monitor_trajectory` (all unmodified) consume the bridged claim. Only the
  bridge itself (`Verdict -> TwoSidedClaim`) is new code.

## 3. Approach

`certkit_bridge.Verdict(name, rule, ok, line, lo, hi)` is a **producer-side** record: `lo`/`hi` are
populated from the numbers passed into `emit()` *before* the independent checker runs, so they are
present and finite even when `ok=False` (confirmed by probing: H2's `certificate_sector` Verdict is
`ok=False` with `lo=-1.8523881735695835, hi=-1.8523881735695833`, a deceptively tight-looking pair
attached to an ABSTAIN). Trusting `lo`/`hi` on their own would silently forward an unverified
number into an actuation plan. The bridge therefore gates on `ok` alone:

- `ok=False` -> an **abstaining** `TwoSidedClaim`, carrying the checker's own `line` as the reason
  (whatever it says — ABSTAIN or CRASH — never re-derived or paraphrased).
- `ok=True` -> a **non-abstaining** `TwoSidedClaim`. D1's `Trajectory.max_model_error` is a single
  non-negative margin, not an energy; the certified enclosure's **width** (`hi - lo`) is what a
  verified certificate actually promises — "the true value is confined to an interval this wide" —
  so it is the natural stand-in for a model-error bound. The bridge emits
  `TwoSidedClaim(lower=0.0, upper=hi-lo, abstained=False)`, matching the shape D3's own
  `run_failclosed_sweep` baseline already uses (`TwoSidedClaim(lower=0.0, upper=0.05, ...)`).

Reference: `certkit_bridge.run_case("H2", ...)`'s own two `ok=True` certificates for the same
molecule give two different, both-correct outcomes downstream — this is the reference the gates
check against, not a synthetic one:

- `certificate_temple` (tight inertia-discharged bound): width **3.7e-9**, far under D3's
  `SAFETY_FLOOR = 0.15` -> the monitor **CONTINUEs**.
- `certificate_gershgorin` (loose, gap-free bound): width **0.161**, *above* `SAFETY_FLOOR` -> the
  monitor **HALTs**, even though the certificate is genuinely `VERIFIED`.

That second case is the finding this spec makes visible rather than smoothing over: **a verified
certificate is not the same thing as a safe-to-actuate one.** D3's floor is a second, independent
gate on top of certkit's own soundness check, and a real chemistry certificate exercises both
outcomes without any synthetic construction.

## 4. Public interface

```
robot.chem_bridge.verdict_to_claim(verdict: certkit_bridge.Verdict) -> robot.planner.TwoSidedClaim
robot.chem_bridge.plan_from_verdict(
    verdict: certkit_bridge.Verdict, origin: Point, target: Point, n_waypoints: int = 5,
) -> Trajectory | PlanRefusal
```

`robot/__init__.py` does not re-export `chem_bridge` (it imports `certkit_bridge`, which requires
the full solver stack plus the `certkit` extra) — same reasoning D2 already applied to `workspace`.

## 5. Acceptance criteria (validation gates)

- **G1 — A `VERIFIED` certificate plans and continues.** H2's `certificate_temple` Verdict
  (`ok=True`, width 3.7e-9) maps to a non-abstaining claim; `plan_from_verdict` returns a
  `Trajectory`; feeding that trajectory through D2's `classify_trajectory` and D3's
  `monitor_trajectory` (nominal envelope, no injected events) yields `Verdict.INSIDE` and
  `MonitorAction.CONTINUE`.
- **G2 — An `ABSTAIN` certificate refuses, and the reason is the checker's own line.** H2's
  `certificate_sector` Verdict (`ok=False`) maps to an abstaining claim whose `reason` is exactly
  `verdict.line`; `plan_from_verdict` returns a `PlanRefusal`, never a `Trajectory`.
- **G3 — A `VERIFIED` but loose certificate still trips the fail-closed floor.** H2's
  `certificate_gershgorin` Verdict (`ok=True`, width 0.161 > `SAFETY_FLOOR`) maps to a
  non-abstaining claim and *does* produce a `Trajectory` (D1 does not know about D3's floor) —
  but `monitor_trajectory` on that trajectory returns `MonitorAction.HALT`. Definition of done:
  this is the epic's "every hop carries a certificate that could have failed" property,
  demonstrated with a real certificate rather than asserted.
- **G4 — The bridge trusts `ok`, never `lo`/`hi`'s mere presence.** A synthetic
  `Verdict(ok=False, lo=-1.0, hi=-1.0, ...)` — a deceptively tight, perfectly finite bracket
  attached to a non-verified verdict, reproducing the exact shape probing found on
  `certificate_sector` — still maps to an abstaining claim. This is the regression gate for the
  footgun in `certkit_bridge.Verdict`'s own field population (§3): a bridge that checked
  "`lo`/`hi` are not `None`/finite" instead of `ok` would pass G1-G3 and still be wrong.

> Definition of done: G1-G4 green in `tests/test_robot_chem_bridge_d4_spec.py`.

## 6. Implementation plan (test-first)

1. `tests/test_robot_chem_bridge_d4_spec.py` encoding G1-G4 (written first, failing — module
   does not exist).
2. `robot/chem_bridge.py` — `verdict_to_claim`, `plan_from_verdict`, built on D1's
   `plan_pipette_trajectory`/`TwoSidedClaim` and consuming `certkit_bridge.Verdict` unmodified.
3. `make gates` (this file needs the `certkit` extra and the full solver stack — same isolation
   requirement as `test_certkit_regression_gate_spec.py`; skips cleanly if the extra is absent,
   same as that file, and raises in CI if so).

## 7. Out of scope

- Vendoring or reproducing the real `certabstain` package — `TwoSidedClaim` is still D1's stand-in.
- Choosing *which* certificate route (temple vs. gershgorin) to actuate on when a case emits more
  than one `ok=True` Verdict — that policy decision is a caller concern, not this bridge's.
- Recalibrating `SAFETY_FLOOR` against real chemistry-certificate widths; §3's G3 outcome is
  reported honestly, not treated as a bug to fix by raising the floor.
- A real robot, sensor, or wet-lab synthesis step.

## 8. Caveats and risks

- **R1:** "certified enclosure width = model error" is a stated modeling choice (§3), not derived
  from a physical map between quantum-chemistry energy uncertainty and pipette positioning error.
  Do not carry this mapping to a real actuator without re-deriving it.
- G3's HALT is specific to `SAFETY_FLOOR = 0.15` (itself already flagged in
  `SPEC_robot_failclosed_d3.md` R1 as calibrated against fixtures, not hardware) colliding with
  this particular Gershgorin width for this particular molecule/active space. It demonstrates that
  verified-but-loose certificates can trip the floor; it is not a general claim that the
  Gershgorin route is unsafe.
- G1-G4 exercise `run_case("H2", ...)` once per test session (module-level cache, ~a few seconds);
  this is a real, small quantum-chemistry computation, not a mock.

## 9. Deliverables

- `robot/chem_bridge.py` — `verdict_to_claim`, `plan_from_verdict`.
- `tests/test_robot_chem_bridge_d4_spec.py` — gates G1-G4.
