# chem-0cb — Two robot D2/D3 spec tests import certkit unconditionally, breaking full-suite collection

Status: **closed, fixed**. Acceptance criteria met exactly as stated.

## What changed

- `tests/test_robot_workspace_d2_spec.py` — moved `from certkit.interval import Iv` behind the
  same guard used by `test_certkit_regression_gate_spec.py` /
  `test_robot_chem_bridge_d4_spec.py`: `importlib.util.find_spec("certkit") is None` →
  `pytest.skip(..., allow_module_level=True)` (with the pre-existing repo convention of a hard
  `RuntimeError` under `CI=true`, so a broken CI environment still fails loud instead of silently
  skipping). `Iv`, `robot.planner`, and `robot.workspace` imports moved below the guard with
  `# noqa: E402` (matching the sibling files' style).
- `tests/test_robot_failclosed_d3_spec.py` — identical guard applied; `Iv`, `robot.failclosed`,
  `robot.planner`, `robot.workspace` imports moved below it.
- No other files touched. `specs/` untouched — this is a test-collection bug, not a spec change.

Diff is two files, 34 insertions / 7 deletions total (`git diff --stat` on the commit below).

## Pre-registered acceptance criteria (from `bd show chem-0cb`) and pass/fail

1. **Both files guard their certkit import the same way the two sibling gate files already do.**
   ✅ Done — same `find_spec` + `CI` hard-fail + `pytest.skip(allow_module_level=True)` pattern,
   byte-for-byte structurally identical to `test_robot_chem_bridge_d4_spec.py`'s guard.
2. **`uv run pytest --collect-only -q tests/` reports zero collection errors without the certkit
   extra installed.** ✅ Confirmed. Environment for this run: `uv sync --extra dmrg --extra test`
   (no `--extra certkit`); `uv run python -c "import certkit"` fails with
   `ModuleNotFoundError: No module named 'certkit'`, confirming the repro precondition.

   Command and result:
   ```
   $ uv run pytest --collect-only -q tests/
   ...
   810 tests collected in 2.27s
   ```
   Exit code `0`. Before the fix (reproduced first, matching the bead's own reproduction): 2
   collection errors (`tests/test_robot_workspace_d2_spec.py:4` and
   `tests/test_robot_failclosed_d3_spec.py:3`, both `ModuleNotFoundError: No module named
   'certkit'`).

## Additional verification (not required by the acceptance criteria, done for confidence)

- `uv run pytest -q tests/test_robot_workspace_d2_spec.py
  tests/test_robot_failclosed_d3_spec.py` → `2 skipped in 0.03s` (clean skip, not an error, and
  the modules no longer error at import time even when run directly).
- `uv run pytest -q tests/test_certkit_regression_gate_spec.py
  tests/test_robot_chem_bridge_d4_spec.py` → `2 skipped in 0.22s` — the two sibling files that
  already had the guard are unaffected, so all four certkit-gated spec files now behave
  consistently.
- `uv run ruff check tests/test_robot_workspace_d2_spec.py
  tests/test_robot_failclosed_d3_spec.py` (via `uv sync --extra dev` to get ruff) → `All checks
  passed!`
- `uv run pytest -q tests/test_reference_energies.py` → `9 passed` — sanity check that the
  environment itself (PySCF/qiskit path) is healthy and this change didn't disturb an unrelated
  gate.
- Did not install the `certkit` extra to check the *positive* path (guard falls through, tests
  actually run against real certkit). No network egress to
  `git+https://github.com/B0001/certkit@v0.2.0` was attempted in this container, and the bead's
  acceptance criteria only concern the no-certkit collection path. The positive path is
  structurally identical code to the two already-guarded sibling files, which do exercise it in
  environments where certkit is installed — I did not re-verify that path here since it was not
  touched by this bead and is out of scope.

## What I decided not to do, and why

- Did not run the full `make gates` / `scripts/run_gates.sh` (all ~90 spec files, including DMRG
  isolation). The acceptance criteria is specifically about collection, not full-suite execution,
  and a full run would burn a large chunk of the shared 8-CPU/16GB compute envelope on gates
  wholly unrelated to this bug. Ran targeted gates instead (see above).
- Did not touch `pyproject.toml`'s certkit extra pin or add certkit to the `test` extra — the bug
  is that the guard was missing, not that certkit should be a hard dependency; the existing design
  (optional extra, skip when absent) is correct and this fix just makes two files conform to it.

## What I could not verify

- The certkit-installed / non-skip code path in the two fixed files (i.e., that `Iv`, D2, D3
  classify/monitor logic still work correctly when certkit *is* present) — unverified in this
  session, no certkit extra installed. Low risk: the only change below the guard is dedent-free
  relocation of pre-existing import statements plus `# noqa: E402` comments; no logic in the test
  bodies was touched.

## Commands to reproduce

```bash
cd /workspace
uv sync --extra dmrg --extra test          # matches documented setup, no certkit
uv run pytest --collect-only -q tests/     # -> 810 tests collected, 0 errors, exit 0
uv run pytest -q tests/test_robot_workspace_d2_spec.py tests/test_robot_failclosed_d3_spec.py
                                            # -> 2 skipped
```

## Git

Committed on `sandbox/chem-0cb` (this branch was cut from `main` for this task):

```
e88a66f Guard certkit imports in robot D2/D3 spec gates
```

Tree is otherwise clean (`git status --short` shows nothing outstanding). Per this run's
instructions, I did not push — the host pushes and a human merges.

## Bead

Closed `chem-0cb` with the evidence above (`bd close chem-0cb`).
