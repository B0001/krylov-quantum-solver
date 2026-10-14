# Handoff: chem-s6c

## Disposition

**Fixed and closing.** `reachability.py`'s own `scf_symmetry_status` and `hf_orbital_irreps` are
now deterministic across ambient thread count, by the same pin/restore pattern chem-7cx applied to
`build_molecular_hamiltonian`. Gate file green: `tests/test_reachability_scf_determinism_spec.py`,
5/5 passed. `ruff check .` clean. Full spec: `specs/SPEC_reachability_scf_determinism.md`.

One out-of-scope finding filed as a new bead (**chem-50w**), not fixed here — see "Out-of-scope
findings" below.

## Root cause

Identical mechanism to chem-7cx (`specs/SPEC_scf_branch_determinism.md`): PySCF's own
multi-threaded OpenMP integral evaluation is not bit-associative. On a near-degenerate geometry,
this can tip a free (`symmetry=False`) RHF solve onto a genuinely different locally-stable branch
between otherwise-identical calls, depending only on how many threads happen to be configured —
not just low-digit numerical noise in the energy. `reachability.py` builds its own `gto.M`/
`scf.RHF` objects directly in two functions `build_molecular_hamiltonian` never touches:
`scf_symmetry_status` (free + symmetry-enforced RHF) and `hf_orbital_irreps` (its own
symmetry-enforced RHF for orbital labelling, plus an internal call to `scf_symmetry_status`).

Directly reproduced before fixing anything — square H4, sto-3g, `a = 1.10`, `conv_tol = 1e-6`:

```
$ for n in 1 2 4 8; do OMP_NUM_THREADS=$n uv run python -c "
import reachability
print($n, reachability.scf_symmetry_status('H 0 0 0; H 1.10 0 0; H 1.10 1.10 0; H 0 1.10 0', conv_tol=1e-6))"
done
1 (True, -0.07159833033797103)
2 (False, -4.026112776500668e-12)
4 (True, -0.07159833033797103)
8 (True, -0.07159850445902016)
```

`OMP_NUM_THREADS=2` flips the **classification itself** (`is_broken`), not just the tail digits of
`dE`. This is the exact mechanism chem-7cx's own investigation notes (filed while looking at
chem-3fd) predicted would be present in `scf_symmetry_status` — this bead confirms it directly.

## Fix

`reachability.py` — pin `pyscf.lib.num_threads(1)` around every RHF `.kernel()` call in both
functions, restore the ambient count in a `finally`:

```python
def scf_symmetry_status(atom, basis="sto-3g", conv_tol=TIGHT_SCF_CONV_TOL):
    from pyscf import gto, lib, scf

    ambient_threads = lib.num_threads()
    lib.num_threads(1)
    try:
        free = scf.RHF(gto.M(atom=atom, basis=basis, symmetry=False, verbose=0))
        free.conv_tol = conv_tol
        free.kernel()
        enforced = scf.RHF(gto.M(atom=atom, basis=basis, symmetry=True, verbose=0))
        enforced.conv_tol = conv_tol
        enforced.kernel()
    finally:
        lib.num_threads(ambient_threads)
    d_e = float(free.e_tot - enforced.e_tot)
    return (d_e < -SYMMETRY_BREAK_TOL), d_e
```

`hf_orbital_irreps` calls `scf_symmetry_status` internally (already pinned by the above) and
additionally pins its own extra `mf_sym.kernel()` call used for the labelling step, with the same
pin/restore pattern.

`reachability._driver_orbitals`'s own, separate unpinned `drv.run_pyscf()` call is **not** touched
— see "Out-of-scope findings" below.

## SDD artifacts

- `specs/SPEC_reachability_scf_determinism.md` — full spec: goal, background/honest framing
  (explicitly separates "now deterministic" from "now physically correct" — the latter is
  chem-g14's territory, untouched here), approach, G1-G4 acceptance gates, out-of-scope, caveats
  (R1-R3), results, deliverables.
- `tests/test_reachability_scf_determinism_spec.py` — G1-G4, pyscf-only (imports `reachability`, a
  top-level script, not `hybrid_quantum_solver`; does not import block2/qiskit-aer), runs in its
  own process per the `test_*_spec.py` isolation convention.

Verbatim gate run (this session, final confirmation):

```
$ uv run pytest tests/test_reachability_scf_determinism_spec.py -v
tests/test_reachability_scf_determinism_spec.py::test_G1_deterministic_across_ambient_thread_counts PASSED [ 20%]
tests/test_reachability_scf_determinism_spec.py::test_G2_deterministic_across_repeats_separate_processes PASSED [ 40%]
tests/test_reachability_scf_determinism_spec.py::test_G3_scf_symmetry_status_deterministic_same_process_and_threads_restored PASSED [ 60%]
tests/test_reachability_scf_determinism_spec.py::test_G3_hf_orbital_irreps_deterministic_same_process_and_threads_restored PASSED [ 80%]
tests/test_reachability_scf_determinism_spec.py::test_G4_affected_witnesses_stabilize_pass_or_fail PASSED [100%]
========================= 5 passed in 68.96s (0:01:08)
```

```
$ uv run ruff check .
All checks passed!
```

## Pre-registered acceptance criteria and outcome

| Gate | Claim | Outcome |
|---|---|---|
| G1 | `scf_symmetry_status` on the square-H4 witness (`a=1.10`, `conv_tol=1e-6`) is bit-identical across fresh subprocesses at `OMP_NUM_THREADS` in `{1,2,4,8}` | **PASS** — before the fix, thread count 2 gave `is_broken=False` vs. `True` for 1/4/8; after, all four agree |
| G2 | Ten fresh subprocess calls at the default ambient thread count are pairwise bit-identical | **PASS** |
| G3 | Ten same-process, uncached calls to each of `scf_symmetry_status` and `hf_orbital_irreps` are pairwise bit-identical, and `pyscf.lib.num_threads()` is restored to its pre-call value after every call | **PASS** (both functions) |
| G4 | The 5 `tests/test_symmetry_reachability_spec.py` gates that call these two functions directly on the square-H4 family stop flipping pass/fail across 5 repeated full-file subprocess runs | **PASS** — same 5/14 fail every time, before and after the fix (see "What this does NOT claim") |

## What this does NOT claim (read before reusing these numbers)

- **Not a correctness fix for which branch is "right."** On this host, the now-deterministic
  answer for `scf_symmetry_status` at square H4 `a = 1.05/1.10/1.35`, `conv_tol = TIGHT_SCF_CONV_TOL`
  (the default), is `is_broken = True` — which disagrees with `specs/SPEC_symmetry_reachability.md`'s
  macOS-recorded witness table (`available` at those geometries). That is the **same**
  cross-platform/pyscf-version divergence `specs/SPEC_scf_branch_determinism.md` already found and
  assigned to **chem-g14** — not a new defect, and not something this bead's fix changes.
  `tests/test_symmetry_reachability_spec.py` still fails the same 5/14 tests after this fix as
  before; G4 only asserts that failure is now **stable**, not that it's correct.
- **`test_eigenstate_reachability_spec.py` is unaffected by this fix, by design — not a residual
  bug.** Verified by grep: no gate in that file imports `scf_symmetry_status`,
  `hf_orbital_irreps`, or `symmetry_filter_available`. This contradicts the bead's own parenthetical
  ("affects test_G4 ... uses scf_symmetry_status via hf_orbital_irreps through the diagnostic
  path") — that claim does not match current code; `test_G4_*` there goes through
  `orbital_parity_bounds` → `_driver_orbitals`, a different, unpinned function (see below).
  Measured directly: one full run of that file with the fix reverted and one with it applied gave
  16 vs. 12 failures, but the specific failing IDs differed in ways inconsistent with this fix
  being the cause (e.g. `test_G2_decision_equals_symmetry_adapted_fci[1.05-1e-09]` failed only in
  the *with-fix* run) — i.e. that file's pass/fail pattern is independently flaky from
  `build_molecular_hamiltonian`'s own still-unpinned nondeterminism on this branch (chem-7cx's fix
  lives on sibling branch `sandbox/chem-7cx`, not merged here), not from this change.
  `tests/test_reachability_tolerance_spec.py` shows the identical pattern (also doesn't import
  these two functions; three repeated unpinned-baseline runs gave 6, 4, and 7 failures with
  different IDs each time) — logged as the same pre-existing noise, not evidence against this fix.

## Out-of-scope findings filed separately

- **chem-50w** (new, filed this session, open, P2): `reachability._driver_orbitals`'s own
  independent, unpinned `drv.run_pyscf()` call (used by `orbital_parity_bounds`/
  `hf_symmetry_sector`, not by the two functions this bead fixes) shares the same thread-
  nondeterminism mechanism. Evidence it's live: chem-0yh's investigation (recorded in chem-s6c's
  notes) found `test_eigenstate_reachability_spec.py::test_G2_decision_equals_symmetry_adapted_fci[1.19]`'s
  `in_sector` jumping from length 10 to 36 across fresh-process repeats, even with
  `build_molecular_hamiltonian`'s output and `_fci`'s output both bit-identical across the same
  repeats. Not fixed here because `_driver_orbitals` already carries chem-g14's own bounded
  retry + discrete sign-gauge search machinery for a related mechanism, and a drive-by three-line
  pin there would interact with that existing logic in ways that deserve their own
  before/after-retry-count measurement and spec, not a change riding on this bead's title.

## What could not be verified

- Only tested on this container (Linux x86_64, 8 cores, OpenBLAS,
  `Linux 42e21aeeac3f 7.0.14-linuxkit #1 SMP PREEMPT Fri Sep 18 10:19:32 UTC 2026 x86_64 GNU/Linux`).
  Not cross-checked on a second machine or BLAS backend (same caveat
  `SPEC_scf_branch_determinism.md` R1 records for `build_molecular_hamiltonian`).

## Files changed

- `reachability.py` — the fix (both functions).
- `specs/SPEC_reachability_scf_determinism.md` — new spec.
- `tests/test_reachability_scf_determinism_spec.py` — new gate file, G1-G4.
- `sandbox-handoffs/chem-s6c.md` — this file.
- `.beads/issues.jsonl` — bead state (chem-s6c claimed/closed, chem-50w filed) via `bd`'s own export.

Not included: `.claude/settings.local.json` (pre-existing, unrelated, untracked).

## Suggested next commands (not run by me — git policy for this session is commit-only, no push)

```bash
git status
git log --oneline -5
```
