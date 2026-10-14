# SPEC: `reachability.py`'s own RHF calls (`scf_symmetry_status`, `hf_orbital_irreps`) are deterministic regardless of ambient thread count

**Status:** IMPLEMENTED — gate green (`tests/test_reachability_scf_determinism_spec.py`).

---

## 1. Goal

**Claim:** `reachability.scf_symmetry_status(atom, basis, conv_tol)` and
`reachability.hf_orbital_irreps(atom, basis, conv_tol)` return **bit-identical** results for a
fixed `(atom, basis, conv_tol)`, regardless of the ambient `pyscf.lib.num_threads()`/
`OMP_NUM_THREADS` setting the caller leaves in place, or of repeated calls at a fixed setting.

Before the fix, this was false. Direct repro, square H4 (sto-3g, `a = 1.10`, `conv_tol = 1e-6`,
captured 2026-10-14 on this container, `scf_symmetry_status` unpinned):

| `OMP_NUM_THREADS` | `(is_broken, dE)` |
|---|---|
| 1 | `(True, -0.07159833033797103)` |
| 2 | `(False, -4.026112776500668e-12)` |
| 4 | `(True, -0.07159833033797103)` |
| 8 | `(True, -0.07159850445902016)` |

Thread count 2 flips the **classification itself** (`is_broken` False instead of True) on the exact
same `(atom, basis, conv_tol)` input — not just noise in the last few digits of `dE`. This is
chem-s6c's own repro of the mechanism chem-7cx's notes (filed while investigating chem-3fd)
predicted would be present here: "square H4 a=1.20's 'broken RHF' classification
(`scf_symmetry_status`) flipped True/False across runs."

## 2. Background and honest framing

- Filed as **chem-s6c** (P2), a direct follow-up to **chem-7cx**
  (`specs/SPEC_scf_branch_determinism.md`), which fixed the identical mechanism in
  `hybrid_quantum_solver.molecular_hamiltonian.build_molecular_hamiltonian`: PySCF's own
  multi-threaded OpenMP integral evaluation is not bit-associative, so a near-degenerate SCF can
  land on a different locally-stable solution between otherwise-identical calls, depending only on
  how many threads happen to be configured when the call runs. chem-7cx isolated the mechanism by
  toggling `OMP_NUM_THREADS` vs `OPENBLAS_NUM_THREADS` independently (only the former reproduces
  the flip) and deliberately left `reachability.py`'s own raw `scf.RHF(gto.M(...))` calls out of
  that bead's scope, filing this one instead (see that spec's §7, "Out of scope").
- **`reachability.py` is explicitly NOT part of the `hybrid_quantum_solver` package** (CLAUDE.md,
  Architecture) and builds its own `gto.M`/`scf.RHF` objects directly in two functions
  `build_molecular_hamiltonian` never touches: `scf_symmetry_status` (two RHF solves, free and
  symmetry-enforced) and `hf_orbital_irreps` (one more RHF solve, symmetry-enforced, used for
  orbital labelling; it also calls `scf_symmetry_status` internally).
- **What you can claim if the gate passes:** any future call to `scf_symmetry_status` or
  `hf_orbital_irreps` for a fixed `(atom, basis, conv_tol)` returns the same answer every time, on
  this host, independent of how many threads PySCF happens to be configured for when it is called.
  `specs/SPEC_symmetry_reachability.md`'s own `scf_symmetry_status`/`hf_orbital_irreps` call sites
  (`symmetry_filter_available`, and every `test_symmetry_reachability_spec.py` gate) stop depending
  on a coin flip for their answer.
- **What you cannot claim:** that the now-deterministic branch is the "physically correct" one, or
  that it matches `specs/SPEC_symmetry_reachability.md`'s documented 9-geometry table (which records
  `a = 1.05, 1.10, 1.35` as symmetric/available). On **this** host, the deterministic answer for
  `scf_symmetry_status(_square_h4(1.05 | 1.10 | 1.35))` at the default `conv_tol = TIGHT_SCF_CONV_TOL
  = 1e-13` is consistently `is_broken = True` (measured §9) — the SAME cross-platform/pyscf-version
  divergence `specs/SPEC_scf_branch_determinism.md` already found and left to **chem-g14**, not a
  new defect. `tests/test_symmetry_reachability_spec.py::test_G1...` and
  `::test_G3_available_means_the_symmetric_solution_is_the_RHF_solution[1.05|1.10|1.35]` fail
  identically before and after this fix (§9) — this spec makes that failure **stable**, it does not
  turn it green. That remains chem-g14's territory.
- This fix does **not** touch `reachability._driver_orbitals`'s own unpinned `drv.run_pyscf()` call
  (used by `orbital_parity_bounds`/`hf_symmetry_sector`, not by `scf_symmetry_status`/
  `hf_orbital_irreps`) — chem-0yh's notes flagged that as a second, independent nondeterminism
  source reaching `test_eigenstate_reachability_spec.py::test_G2[1.19]`, already carrying chem-g14's
  own bounded retry/gauge-search logic for a related (but distinct) mechanism. Investigated here
  (§8 R2) and left as a separate, newly-filed follow-up — not the same functions this bead's title
  names, and `_driver_orbitals` already has its own blast-radius-scoped fix machinery that a
  three-line thread pin would interact with in ways worth its own spec.

## 3. Approach

Pin `pyscf.lib.num_threads(1)` around each RHF `.kernel()` call inside `scf_symmetry_status` and
`hf_orbital_irreps`, restoring the ambient count in a `finally` — the identical pattern chem-7cx
used for `build_molecular_hamiltonian` and chem-y6w for `benchmark_nbn._tight_scf`.
`hf_orbital_irreps` calls `scf_symmetry_status` internally (already pinned) and additionally pins
its own extra `mf_sym.kernel()` call for the labelling step. **Reference for this claim:** repeated,
varied-ambient-thread reruns of the exact functions under test, compared bit-for-bit against each
other — a determinism claim about the code, not an energy claim about the molecule, so no FCI/DMRG
ground truth is needed (same reasoning `SPEC_scf_branch_determinism.md` §3 and
`SPEC_nbn_scf_determinism.md` §3 used).

## 4. Public interface

No new public interface and no signature change. `scf_symmetry_status` and `hf_orbital_irreps`
pin/restore `pyscf.lib.num_threads` internally around their RHF `.kernel()` calls.

## 5. Acceptance criteria (validation gates)

`tests/test_reachability_scf_determinism_spec.py` — pyscf only (imports `reachability`, which is a
top-level script, not `hybrid_quantum_solver`); **must not** share a process with block2/DMRG tests
— own file, per the `test_*_spec.py` isolation convention (CLAUDE.md).

- **G1 — DEFINITION OF DONE: cross-ambient-thread determinism.** Fresh `scf_symmetry_status` calls
  on square H4 (`a = 1.10`, `conv_tol = 1e-6`, the exact witness in §1) at ambient
  `OMP_NUM_THREADS` in `{1, 2, 4, 8}`, each in its own subprocess, give a **bit-identical**
  `(is_broken, dE)` (`==`, not `pytest.approx`) across all four. Before the fix this witness flipped
  `is_broken` itself between thread counts 1/4/8 (`True`) and 2 (`False`) — not just digit noise.
- **G2 — repeat-call determinism at one ambient setting (regression floor).** Ten fresh
  `scf_symmetry_status` calls, each its own subprocess at the default (unset) ambient thread count,
  on the same witness, are pairwise bit-identical.
- **G3 — same-process repeat-call determinism, both functions, thread count restored.** Ten fresh,
  uncached calls to `scf_symmetry_status` AND ten to `hf_orbital_irreps`, each inside ONE process
  (no subprocess spawn), are pairwise bit-identical, and leave the process's ambient
  `pyscf.lib.num_threads()` unchanged after every call — the `finally`-restore is load-bearing: a
  caller elsewhere in the same process must not find its own thread count silently pinned to 1.
  Run on the witness (`a = 1.10`, broken) and on LiH (`a` n/a, not broken, exercises the
  `hf_orbital_irreps` labelling branch with a real irrep list returned).
- **G4 — the affected spec's gates stop flipping pass/fail.** The five
  `tests/test_symmetry_reachability_spec.py` gates that call `scf_symmetry_status`/
  `hf_orbital_irreps` directly on the square-H4 family
  (`test_G1_filter_is_available_exactly_when_the_artifact_is_present`,
  `test_G3_unavailable_means_RHF_found_a_lower_solution[1.2]`,
  `test_G3_available_means_the_symmetric_solution_is_the_RHF_solution[1.05|1.10|1.35]`) are run as
  five separate `pytest` subprocess invocations of the whole file; this gate asserts the pass/fail
  OUTCOME of each is identical across all five (it does not assert which outcome — chem-g14 governs
  whether it should be GREEN on this host; this spec only requires it stop being a coin flip, same
  framing as `SPEC_scf_branch_determinism.md`'s G4).

## 6. Implementation plan (test-first)

1. `tests/test_reachability_scf_determinism_spec.py` encoding G1–G4 — G1 RED before the fix
   (verified manually, §1's captured table, by reverting the `lib.num_threads` pin and rerunning
   with `OMP_NUM_THREADS` varied).
2. `reachability.py`: pin `lib.num_threads(1)` around each RHF `.kernel()` call in
   `scf_symmetry_status` and `hf_orbital_irreps`, restore the ambient count in a `finally`.
3. Green (§9).

## 7. Out of scope

- **`reachability._driver_orbitals`'s own unpinned `drv.run_pyscf()` call** (used by
  `orbital_parity_bounds`/`hf_symmetry_sector`, not by the two functions this spec covers) —
  chem-0yh's notes flagged this as a second, independent nondeterminism source reaching
  `test_eigenstate_reachability_spec.py::test_G2[1.19]`. Investigated (§8 R2): extending the same
  pin there is plausible but interacts with `_driver_orbitals`'s existing bounded retry + discrete
  sign-gauge search (chem-g14, `specs/SPEC_rebuild_gauge_search.md`) in a way that deserves its own
  spec rather than a drive-by three-line change riding on this bead's title. Filed as
  **chem-50w**.
- Re-deriving which branch is "correct" relative to `specs/SPEC_symmetry_reachability.md`'s
  documented 9-geometry table, or to `specs/SPEC_eigenstate_reachability.md`'s macOS-measured
  witness values — not this spec's claim (see §2; chem-g14's territory).
- `test_eigenstate_reachability_spec.py` as a whole: none of its gates import or call
  `scf_symmetry_status`/`hf_orbital_irreps` (verified by grep, §8 R1) — its own flakiness is driven
  entirely by `build_molecular_hamiltonian` (chem-7cx, fixed on a sibling branch not yet merged
  here) and `_driver_orbitals`/`_fci` (chem-g14/chem-0yh), none of which this bead touches.
- Cross-host / cross-BLAS-backend verification (`SPEC_scf_branch_determinism.md` R1's caveat
  applies identically here).

## 8. Caveats and risks

- **R1:** verified by grep that no gate in `tests/test_eigenstate_reachability_spec.py` imports
  `scf_symmetry_status`, `hf_orbital_irreps`, or `symmetry_filter_available` — the bead's own
  parenthetical ("affects test_G4 ... uses scf_symmetry_status via hf_orbital_irreps through the
  diagnostic path") does not match the current code (`test_G4_*` there calls
  `orbital_parity_bounds`, which goes through `_driver_orbitals`, a different, unpinned function —
  see §7). Measured directly: running `tests/test_eigenstate_reachability_spec.py` once without
  this fix and once with it gave 16 vs 12 failures respectively, but the specific failing set
  differed in ways inconsistent with this fix being the cause (e.g. a NEW failure,
  `test_G2_decision_equals_symmetry_adapted_fci[1.05-1e-09]`, appeared only in the "with fix" run)
  — i.e. that file's pass/fail pattern is independently flaky from `build_molecular_hamiltonian`'s
  own still-unpinned nondeterminism on this branch, not from this change. This spec's G4 therefore
  targets `test_symmetry_reachability_spec.py` instead, whose gates do call the two functions this
  spec fixes directly. `tests/test_reachability_tolerance_spec.py` shows the identical pattern
  (also doesn't import these two functions; 3 repeated unpinned-baseline runs gave 6, 4, and 7
  failures with different test IDs each time) — same unrelated, pre-existing
  `build_molecular_hamiltonian` flakiness, not evidence against this fix.
- **R2:** this fix does NOT make `test_symmetry_reachability_spec.py` fully green on this host (see
  §2) — G4 deliberately only checks for stability, not correctness. A rerun that shows any of G1-G4
  flip is itself evidence of a regression, not a reason to retry and move on (same discipline as
  `SPEC_scf_branch_determinism.md` R2 / `SPEC_nbn_scf_determinism.md` R2).
- **R3:** verified only on this container (x86_64, 8 cores, OpenBLAS). Not tested on a second
  machine (`SPEC_scf_branch_determinism.md` R1 records the same caveat for
  `build_molecular_hamiltonian`).

## 9. Results

**Linux x86_64 container, 8 cores, pyscf, 2026-10-14.** `uname -a`:
`Linux 42e21aeeac3f 7.0.14-linuxkit #1 SMP PREEMPT Fri Sep 18 10:19:32 UTC 2026 x86_64 GNU/Linux`.
Wall time for the full gate file: see handoff (`sandbox-handoffs/chem-s6c.md`); G4's five full-file
subprocess spawns of `test_symmetry_reachability_spec.py` dominate the gate file's own wall time.

- `tests/test_reachability_scf_determinism_spec.py -q`: **G1-G4 pass** — see handoff for the
  verbatim line.
- Before the fix (manual check, pin reverted): `OMP_NUM_THREADS` in `{1, 2, 4, 8}` on the witness
  gave `(True, -0.07159833033797103)`, `(False, -4.026112776500668e-12)`,
  `(True, -0.07159833033797103)`, `(True, -0.07159850445902016)` — the exact table in §1, captured
  directly, not inferred.
- After the fix: the same four thread counts give `(True, -0.07159833033797103)` every time (G1).
- `tests/test_symmetry_reachability_spec.py -q`, run three times with the pin reverted and once
  more with it applied: the same 5 of 14 tests fail every time, before and after
  (`test_G1_filter_is_available_exactly_when_the_artifact_is_present`,
  `test_G3_unavailable_means_RHF_found_a_lower_solution[1.2]`,
  `test_G3_available_means_the_symmetric_solution_is_the_RHF_solution[1.05|1.10|1.35]`) — the
  cross-platform divergence chem-g14 already tracks, unaffected (made provably stable, not fixed)
  by this change.

## 10. Deliverables

- `reachability.py` — `scf_symmetry_status` and `hf_orbital_irreps` pin/restore
  `pyscf.lib.num_threads` internally around their RHF `.kernel()` calls.
- `tests/test_reachability_scf_determinism_spec.py` — G1-G4.
- **chem-50w** (new, filed, open) — whether the same pin should extend to `_driver_orbitals`'s
  `drv.run_pyscf()` call, deferred here per §7.
