# chem-qqu handoff

## Bug

`tests/test_hubbard_lieb_wu_spec.py::test_G8_dmrg_free_fermion_and_dimer_anchors` failed in a
fresh checkout (no `./.dmrg_tmp/` yet), passing only on re-run. Reported as reproducing on
`origin/main` too.

## Root cause

`G8` calls `dmrg_energy_extrapolated(..., scratch="./.dmrg_tmp/g8_free")` /
`scratch="./.dmrg_tmp/g8_dimer")` — explicit scratch paths. In
`hybrid_quantum_solver/dmrg_reference.py`'s `_scratch_dir`, the `scratch is not None` branch
passed the path straight to `DMRGDriver` **without creating it**. block2's C++ layer opens
scratch files by path and raises `RuntimeError: SparseMatrix:save_data on
'./.dmrg_tmp/g8_free/F.MPS.KET.0' failed` rather than creating missing parent directories itself
(confirmed by reading `pyblock2/driver/core.py` — `os.makedirs` never appears there; only
`restart_dir` gets one). First run in a fresh tree: `./.dmrg_tmp/g8_free` doesn't exist → fails.
Any later run (this test's own leftover dirs, or another test that happened to create
`./.dmrg_tmp` first) → passes, matching the "fails fresh, passes on re-run" report exactly.

## Fix

`hybrid_quantum_solver/dmrg_reference.py`, `_scratch_dir`: when `scratch` is explicit, `os.makedirs(scratch, exist_ok=True)` before yielding it. Still never torn down (unchanged contract —
callers namespace and reuse their own scratch; see `test_dmrg_scratch_isolation_spec.py` G2,
which already asserts this and still passes). One-line change plus a docstring note explaining
why block2 needs it. Diff:

```diff
     if scratch is not None:
+        os.makedirs(scratch, exist_ok=True)
         yield scratch
         return
```

## Verification

Reproduced the exact failure first, on a truly fresh tree (moved `.dmrg_tmp/` aside):

```
$ mv .dmrg_tmp /tmp/dmrg_tmp_backup
$ uv run pytest tests/test_hubbard_lieb_wu_spec.py -k G8 -v
...
RuntimeError: SparseMatrix:save_data on './.dmrg_tmp/g8_free/F.MPS.KET.0' failed.
1 failed, 4 deselected in 4.82s
```

After the fix, with `.dmrg_tmp/` removed again each time (`rm -rf .dmrg_tmp`) to simulate a fresh
checkout:

```
$ rm -rf .dmrg_tmp && uv run pytest tests/test_hubbard_lieb_wu_spec.py -v
tests/test_hubbard_lieb_wu_spec.py::test_G5_every_vendored_point_passes_the_preregistered_check PASSED
tests/test_hubbard_lieb_wu_spec.py::test_G5_table_covers_the_preregistered_grid PASSED
tests/test_hubbard_lieb_wu_spec.py::test_G6_recorded_tdl_numbers_and_acceptance PASSED
tests/test_hubbard_lieb_wu_spec.py::test_G7_stalled_ladder_is_now_rejected_by_truncation_regime_alone PASSED
tests/test_hubbard_lieb_wu_spec.py::test_G8_dmrg_free_fermion_and_dimer_anchors PASSED
5 passed, 1 warning in 59.06s
```

No re-run needed — single pass on first attempt, which is the acceptance criterion (G8 green on a
fresh checkout, first try).

Regression checks, all also run with `.dmrg_tmp/` removed first:

```
$ rm -rf .dmrg_tmp && uv run pytest tests/test_dmrg_scratch_isolation_spec.py -v
test_G1_default_scratch_is_unique_per_call_and_cleaned_up PASSED
test_G2_explicit_scratch_is_passed_through_and_never_cleaned_up PASSED
test_G3_two_concurrent_dmrg_processes_use_distinct_scratch_and_leave_none_behind PASSED
3 passed in 4.49s
```
(G2 is the gate that specifically pins "explicit scratch untouched, never cleaned up" — still
holds; `os.makedirs(..., exist_ok=True)` doesn't touch an already-existing dir's contents.)

```
$ rm -rf .dmrg_tmp && uv run pytest tests/test_reference_energies.py tests/test_krylov_convergence.py tests/test_noise_resilience.py -q
20 passed, 14 warnings in 14.91s
```

```
$ uv run pytest tests/test_dmrg_reference.py -v
test_fci_energy_matches_casci PASSED
test_reference_energy_auto_returns_exact PASSED
test_dmrg_energy_requires_block2 SKIPPED   (block2 is installed in this env, so the "raises
                                             without block2" gate correctly skips itself)
2 passed, 1 skipped
```

```
$ uv run ruff check hybrid_quantum_solver/dmrg_reference.py
All checks passed!
```

Did not run the full `make gates` / `scripts/run_gates.sh` sweep — that includes hour-scale
DMRG/benchmark gates (hchain TDL, N-B-N ladders) unrelated to this bug and explicitly flagged in
CLAUDE.md as something to leave to the user's own terminal given the compute envelope here (1
worker, ~8 CPU / 16 GB shared). The gates actually exercising the changed code path
(`test_hubbard_lieb_wu_spec.py` G8, `test_dmrg_scratch_isolation_spec.py` G1-G3) were run directly
and are the ones that matter for this fix.

## Hardware / wall time

`uname -a`: `Linux cc48cb1e5c80 7.0.14-linuxkit x86_64` (Docker VM), `nproc`: 8. All test runs
above were under a minute each (G8's DMRG run ~59s, everything else single-digit seconds) — none
crossed the 1-minute reporting threshold except the G8 run itself, which is unrelated to the fix
(same DMRG cost before and after; the fix only changes whether the directory gets created, not
the solve).

## Scope / what I did not do

- Did not touch `DEFAULT_SCRATCH_ROOT`'s (the `scratch=None`) branch — it already calls
  `os.makedirs(DEFAULT_SCRATCH_ROOT, exist_ok=True)` before `tempfile.mkdtemp`, so it was never
  affected by this bug.
- Did not change any caller (`benchmark_hubbard_lieb_wu.py`, `nbn_dmrg_reference.py`) — the bead's
  own suggested fix locations were "the test/fixture, or dmrg_reference's per-call scratch
  handling"; fixing it once in `_scratch_dir` covers every explicit-scratch caller instead of
  patching each test/script site individually, and is where `test_dmrg_scratch_isolation_spec.py`
  already specifies the contract for this helper.
- Did not open a new bead — this was fully in scope and is resolved.

## Files changed

- `hybrid_quantum_solver/dmrg_reference.py` (the fix, 4 lines incl. docstring)
- `.beads/issues.jsonl` (bd metadata from claiming/closing chem-qqu — not hand-edited)

## Git status (tree left dirty per instructions — do not push)

```
 M .beads/issues.jsonl
 M hybrid_quantum_solver/dmrg_reference.py
```

Proposed commit (not run — git policy for this session is commit-but-not-push, see below):

```
git add hybrid_quantum_solver/dmrg_reference.py .beads/issues.jsonl
git commit -m "dmrg_reference: create explicit scratch dir if missing (chem-qqu)

block2 opens scratch files by path and fails with a RuntimeError rather
than creating missing parent directories. test_hubbard_lieb_wu_spec.py's
G8 passes explicit scratch=./.dmrg_tmp/g8_free and g8_dimer, so a fresh
checkout (no .dmrg_tmp/ yet) failed the first run and only passed on
re-run once some earlier test had created the directory tree.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

## Bead status

Closed `chem-qqu` with evidence above (fresh-checkout G8 green on first try, no regressions in
the scratch-isolation or core validation gates).
