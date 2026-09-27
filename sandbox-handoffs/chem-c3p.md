# chem-c3p handoff — concurrent DMRG gates sharing `./.dmrg_tmp` scratch

## Status: fix was already implemented and largely tested by a prior (unfinished) worker
session on this checkout. This session verified it end-to-end against the bead's exact
acceptance criteria, fixed the environment needed to actually run block2, and closes the bead.

## What was already in the working tree at claim time

`hybrid_quantum_solver/dmrg_reference.py` (uncommitted, `git diff` shows it) had already been
changed to:

- Add `DEFAULT_SCRATCH_ROOT = "./.dmrg_tmp"` and a `_scratch_dir(scratch)` contextmanager:
  - `scratch=None` (new default for both `dmrg_energy` and `dmrg_energy_extrapolated`, was
    previously the literal string `"./.dmrg_tmp"`): creates `tempfile.mkdtemp(prefix=f"pid{os.getpid()}_",
    dir=DEFAULT_SCRATCH_ROOT)`, yields it, `shutil.rmtree(..., ignore_errors=True)` on exit.
  - `scratch=<explicit path>`: passed straight through, untouched, never cleaned up (callers that
    already namespace their own scratch — `nbn_dmrg_reference.py`, `benchmark_hubbard_lieb_wu.py`
    — rely on it surviving the call).
- `dmrg_energy` / `dmrg_energy_extrapolated` signatures changed `scratch: str = "./.dmrg_tmp"` ->
  `scratch: Optional[str] = None`, driver construction moved inside `with _scratch_dir(scratch) as
  scratch_dir:`.

`tests/test_dmrg_scratch_isolation_spec.py` (untracked, new) was also already present:
- G1: two `_scratch_dir(None)` calls in the same process get distinct dirs under
  `DEFAULT_SCRATCH_ROOT`, both cleaned up on exit.
- G2: an explicit `scratch=` path is passed through unchanged and NOT cleaned up.
- G3 (requires block2): spawns two real OS *processes*, both calling `dmrg_energy` with the
  default scratch on the 2-site Hubbard dimer, blocks until both have printed their scratch dir
  (so there's a genuine concurrency window), asserts the dirs differ, both children exit 0 with
  the correct energy (`E0 = (U - sqrt(U^2+16t^2))/2`, t=1, U=4 -> -2.472135955...), and asserts no
  stale scratch remains after.

I did not write any of the above — it was in the tree when I claimed the bead. I read it, did not
assume it worked, and verified it below.

## What I did this session

1. **Fixed the environment to actually exercise block2.** The committed `.venv` in this checkout
   is a macOS (`cpython-3.13-macos-aarch64-none`) venv baked by the repo owner's laptop — unusable
   on this Linux x86_64 sandbox (`.venv/bin/python` is a dangling symlink to a `/Users/...` path).
   `uv run` was pointed at `UV_PROJECT_ENVIRONMENT=/tmp/venv` (already set in this container) and I
   ran `uv sync --extra dmrg --extra test --extra dev` there, which pulled a real
   `block2==0.5.3` linux-gnu build (`block2.cpython-313-x86_64-linux-gnu.so`), pyscf, qiskit, pytest,
   ruff. This is local to `/tmp/venv`, nothing under version control changed for this.

2. **Ran the dedicated gate.**
   ```
   $ export UV_PROJECT_ENVIRONMENT=/tmp/venv
   $ uv run pytest -q tests/test_dmrg_scratch_isolation_spec.py -v
   tests/test_dmrg_scratch_isolation_spec.py::test_G1_default_scratch_is_unique_per_call_and_cleaned_up PASSED
   tests/test_dmrg_scratch_isolation_spec.py::test_G2_explicit_scratch_is_passed_through_and_never_cleaned_up PASSED
   tests/test_dmrg_scratch_isolation_spec.py::test_G3_two_concurrent_dmrg_processes_use_distinct_scratch_and_leave_none_behind PASSED
   3 passed in 4.82s
   ```

3. **Ran the bead's literal acceptance criteria**: two DMRG-backed spec gates concurrently, via
   the real `scripts/run_gates.sh` machinery (not a hand-rolled repro), against the two files the
   bead names:
   ```
   $ rm -rf .dmrg_tmp
   $ GATE_RUN="uv run" GATE_JOBS=2 GATE_NO_CACHE=1 \
     GATE_GLOB='tests/test_hchain_largen2_spec.py tests/test_hubbard_lieb_wu_spec.py' \
     bash scripts/run_gates.sh
   gates: 2 files, 2 parallel processes, cache=off
   PASS  tests/test_hubbard_lieb_wu_spec.py  (86s)
   PASS  tests/test_hchain_largen2_spec.py  (304s)
   all spec gates passed
   ```
   While both were running I inspected `.dmrg_tmp/` mid-flight and saw distinct entries: `g8_free`
   (the explicit, caller-owned scratch `test_hubbard_lieb_wu_spec.py` passes) and
   `pid443_254kxk5b` (the auto-generated per-call dir for `test_hchain_largen2_spec.py`, which uses
   the new default). After both finished, `pid443_...` was gone (cleaned up); only the
   pre-existing, intentionally-persistent `g8_free`/`g8_dimer` explicit dirs remained — i.e. no
   *stale* scratch, only caller-owned scratch that was never supposed to be cleaned up by this fix
   (see G2). This is the exact failure mode chem-mjz's session hit (`test_hchain_largen2_spec.py`
   hanging next to another default-scratch DMRG gate) and it no longer collides.

4. **Regression-checked the direct DMRG reference tests** (not new, but touch the same code path):
   ```
   $ uv run pytest -q tests/test_dmrg_reference.py
   2 passed, 1 skipped in 3.63s
   ```

5. **Lint.**
   ```
   $ uv sync --extra dmrg --extra test --extra dev   # pulled ruff==0.15.22 into /tmp/venv
   $ uv run ruff check hybrid_quantum_solver/dmrg_reference.py tests/test_dmrg_scratch_isolation_spec.py
   All checks passed!
   ```

6. **Checked every caller of `dmrg_energy`/`dmrg_energy_extrapolated` in the repo** for a scratch
   arg that could be affected by the `str -> Optional[str] = None` signature change: all pass
   `scratch=` by keyword or omit it entirely; none pass it positionally, so nothing else in the
   repo needed touching (`nbn_dmrg_reference.py`, `benchmark_hubbard_lieb_wu.py`,
   `benchmark_dmrg.py`, `benchmark_dmrg_large.py`, `benchmark_nbn.py`, `study_be2.py`,
   `benchmark_hchain_tdl.py`, and every `tests/test_*_spec.py` that calls either function).

7. Cleaned up the scratch this session's own runs left in `.dmrg_tmp/` (`rm -rf .dmrg_tmp`,
   git-ignored, no tracked-file effect).

## Pre-registered acceptance criteria (from the bead) — result

> Two DMRG-backed spec gates run concurrently (e.g. GATE_JOBS=2 over
> test_hchain_largen2_spec.py + test_hubbard_lieb_wu_spec.py) both pass; scratch dirs are distinct
> per process; no stale scratch left behind.

**PASSED**, all three clauses, per step 3 above.

## What I did not do / could not verify

- I did not run the full `make gates` / `make test` suite (many gates are unrelated to this bead,
  several are long-running per-file processes, and this container is a shared 8-CPU/16GB box —
  running everything wasn't necessary to verify this bead and would have eaten a lot of wall time
  for no additional evidence). I scoped verification to the two DMRG gates the bead names plus the
  dedicated new gate plus `test_dmrg_reference.py`.
- I did not touch or investigate the other uncommitted changes already in the tree
  (`benchmark_nbn.py`, `nbn_dmrg_reference.py`, `nbn_low_spin.py`, `shift_both_sides.py`, the NBN
  specs, the untracked `scripts/gen_nbn_lower_minimum_chk.py`, `specs/SPEC_nbn_scf_determinism.md`,
  `specs/nbn_mp-2634.cif`, `specs/nbn_scf_reference.chk`, `tests/test_nbn_czw_lower_minimum_spec.py`,
  `tests/test_nbn_scf_determinism_spec.py`) — those belong to other in-flight beads (NBN work per
  the recent commit log: chem-y6w/chem-1yr/chem-czw/chem-c3p mentioned together in
  `b5da300`) and are out of this bead's scope. `git status` before and after my session is
  identical for those paths.
- No `specs/SPEC_dmrg_scratch_isolation.md` was written. This repo's own tree has multiple
  `test_*_spec.py` gates with no companion `SPEC_*.md` (e.g. `test_certchem_core_spec.py`,
  `test_hubbard_lieb_wu_spec.py`, `test_excited_bounds_spec.py`) — a bug-fix bead with acceptance
  criteria stated directly in the bead and encoded directly in the gate test's docstring/asserts is
  consistent with that existing precedent, so I didn't add one.
- Hardware/timing for the record: `uname -a` -> `Linux 360ef9904553 7.0.12-linuxkit #1 SMP PREEMPT
  Thu Aug 27 14:02:21 UTC 2026 x86_64 GNU/Linux`; `nproc` -> 8. The concurrent 2-gate run took
  ~304s wall (the slower of the two), matching the bead's own note that `test_hchain_largen2_spec.py`
  alone took 254s — i.e. running next to another DMRG gate no longer costs a 30+ minute hang, just
  the expected extra wall time from sharing 8 CPUs.

## Files changed for this bead (all uncommitted, tree left dirty per git policy)

- `hybrid_quantum_solver/dmrg_reference.py` — the fix (pre-existing in tree, verified).
- `tests/test_dmrg_scratch_isolation_spec.py` — the acceptance gate (pre-existing in tree,
  verified, new/untracked file).

## Suggested commands for a human to commit (NOT run by me — git policy)

```
git add hybrid_quantum_solver/dmrg_reference.py tests/test_dmrg_scratch_isolation_spec.py
git commit -m "dmrg_reference: per-call scratch dir by default, fixing concurrent-gate collisions (chem-c3p)"
```
(Everything else currently dirty/untracked in the tree belongs to other beads and is intentionally
left untouched/uncommitted here.)

## Bead disposition

Closing chem-c3p as done: the fix in the tree satisfies the bead's acceptance criteria exactly as
stated, verified against real block2 concurrency (not a mock), with the exact two gate files named
in the bead.
