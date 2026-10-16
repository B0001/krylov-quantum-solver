# SPEC: `reachability._driver_orbitals`'s own unpinned RHF solve shares chem-7cx/chem-s6c's thread-nondeterminism mechanism

**Status:** IMPLEMENTED — gate green (`tests/test_driver_orbitals_scf_determinism_spec.py`).

---

## 1. Goal

**Claim:** `reachability._driver_orbitals(mh)`'s `drv.run_pyscf()` call, for a **fixed** `mh`,
returns a **bit-identical** rebuild outcome (match/no-match, active MO coefficients, and the number
of retries/sign-gauge searches the retry loop actually runs) regardless of the ambient
`pyscf.lib.num_threads()`/`OMP_NUM_THREADS` setting in effect when `_driver_orbitals` is called, or
of repeated calls at one fixed setting.

Before the fix, this was false. Direct repro, square H4 (sto-3g, `a = 1.35`, `conv_tol = 1e-9`, `mh`
held fixed, captured 2026-10-14 on this container, ambient thread count held fixed at 2, ten
successive calls to the unpinned `_driver_orbitals(mh)`):

| call # | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| found? | yes | yes | yes | yes | **NO** | **NO** | yes | yes | yes | yes |
| `run_pyscf` calls | 3 | 1 | 3 | 1 | 5 | 5 | 2 | 2 | 5 | 5 |

Same fixed `(mh, ambient threads)` input, ten fresh calls, and both the yes/no classification and
the retry count vary run to run — not just low-digit noise. (Separately, holding `mh` fixed and
varying only the ambient thread count at call time also flips the outcome: `a = 1.19`,
`conv_tol = 1e-6`, ambient `{1, 2, 4, 8}` gave found = `{yes, NO, yes, NO}`.)

## 2. Background and honest framing

- Filed as **chem-50w** (P2), explicitly spun out of **chem-s6c**
  (`specs/SPEC_reachability_scf_determinism.md`), which pinned `pyscf.lib.num_threads(1)` around
  `reachability.scf_symmetry_status`'s and `reachability.hf_orbital_irreps`'s own RHF `.kernel()`
  calls but deliberately left `_driver_orbitals`'s independent `PySCFDriver(...).run_pyscf()` call
  untouched, because that function already carries chem-g14's bounded retry + discrete per-orbital
  sign-gauge search (`specs/SPEC_rebuild_gauge_search.md`) for a **related but distinct**
  mechanism (a rebuild landing on a different, or gauge-equivalent, MO basis than the operator it is
  verifying against) and chem-s6c's own notes asked for this interaction to get its own spec and its
  own before/after retry-count measurement rather than a drive-by three-line change.
- `_driver_orbitals` is reached from `orbital_parity_bounds` → `hf_symmetry_sector` →
  `symmetry_allowed` → `reachable_mask`'s diagnostic path, and from
  `tests/test_eigenstate_reachability_spec.py::test_G2_decision_equals_symmetry_adapted_fci` and
  `::test_G4_*`/`test_G5_*`/`test_G6_*` directly. chem-0yh's investigation (recorded in chem-s6c's
  notes) found `test_G2_decision_equals_symmetry_adapted_fci[1.19]` still flaky across fresh-process
  repeats even with `build_molecular_hamiltonian` and `_fci` both already fixed on their own
  branches — `in_sector` jumped from length 10 to length 36 (a real classification flip) while
  `build_molecular_hamiltonian`'s own output matrix (md5) and `_fci`'s own output were both
  bit-identical across the same repeats. `_driver_orbitals`'s unpinned `drv.run_pyscf()` is the only
  remaining unpinned RHF call on that path, confirmed by direct repro above.
- **What you can claim if the gate passes:** for a fixed `mh`, repeated or cross-thread-count calls
  to `_driver_orbitals(mh)` give the identical answer, and the retry loop's own behavior (how many
  `run_pyscf` calls it makes) stops depending on ambient thread count. This removes `_driver_orbitals`
  as a *source* of nondeterminism distinct from whatever nondeterminism the function it rebuilds
  against (`build_molecular_hamiltonian`'s own free RHF call, chem-7cx's territory, not yet merged
  on this branch) still has.
- **What you cannot claim:** that `_driver_orbitals` now always finds a match, or that
  `test_eigenstate_reachability_spec.py`'s remaining failures clear. §9's isolated measurement
  (`mh` built once under pinned threads, so the OUTER build's own separate nondeterminism is held
  out of the picture — see R1) shows the retry/sign-search loop settles to **1 `run_pyscf` call,
  every time**, for every witness this module's gate file exercises. But in the *unisolated* case —
  `mh` built fresh under whatever ambient thread count happens to be live, the ordinary way every
  call site in this repo actually uses it — `_driver_orbitals`'s answer can still depend on ambient
  threads, because `build_molecular_hamiltonian`'s own RHF call (which sets `mh.hf_energy` and
  `mh.qubit_hamiltonian`, the targets `_driver_orbitals` rebuilds against) is **itself** unpinned on
  this branch (chem-7cx's fix lives on a sibling sandbox branch not merged here). `_driver_orbitals`
  now pins its OWN rebuild to a single, deterministic branch (whatever `num_threads=1` produces) —
  but if the thing it is checking against was built on a *different* branch, a deterministic rebuild
  of the wrong branch still correctly, deterministically, refuses. That is chem-7cx's gap, not a new
  one; this spec does not claim to close it and §9 R2 measures it directly so the boundary is not
  hand-waved.

## 3. Approach

Pin `pyscf.lib.num_threads(1)` around the one `drv.run_pyscf()` call inside `_driver_orbitals`'s
retry loop (`reachability.py`, inside `rebuild_from`'s caller), restoring the ambient count in a
`finally` — the identical pattern chem-7cx/chem-s6c/chem-0yh used for the same mechanism elsewhere.
The sign-gauge search itself makes **no** further SCF call (`rebuild_from` only calls
`drv.to_problem()` on already-computed `mo_coeff`, a fresh AO2MO transform, not a new RHF solve), so
nothing else in the function needs touching. **Reference for this claim:** repeated, varied-ambient-
thread reruns of `_driver_orbitals` itself (with `mh` held fixed, built once under a controlled
thread count to isolate this function from `build_molecular_hamiltonian`'s own separate,
out-of-scope nondeterminism), compared bit-for-bit — a determinism claim about this code, not an
energy claim about the molecule, so no FCI/DMRG ground truth is needed (same reasoning
`SPEC_reachability_scf_determinism.md` §3 used).

## 4. Public interface

No new public interface and no signature change. `_driver_orbitals` pins/restores
`pyscf.lib.num_threads` internally around its one `drv.run_pyscf()` call.

## 5. Acceptance criteria (validation gates)

`tests/test_driver_orbitals_scf_determinism_spec.py` — pyscf + qiskit-nature only (imports
`reachability` and `hybrid_quantum_solver.molecular_hamiltonian`), **must not** share a process
with block2/DMRG tests — own file, per the `test_*_spec.py` isolation convention.

- **G1 — DEFINITION OF DONE: cross-ambient-thread determinism, fixed `mh`.** `mh` is built once
  (under `pyscf.lib.num_threads(1)`, to hold the comparison target itself fixed and isolate this
  gate from `build_molecular_hamiltonian`'s own separate nondeterminism — see §2/R1) for each of the
  three witnesses (`a = 1.19, conv_tol = 1e-6`; `a = 1.10, conv_tol = 1e-6`; `a = 1.35,
  conv_tol = 1e-9`). `_driver_orbitals(mh)` is then called once per ambient thread count in
  `{1, 2, 4, 8}` (each its own subprocess, `lib.num_threads` set before the call): the match/no-match
  outcome AND the active MO coefficients returned (`np.array_equal`, bit-identical) are the same
  across all four, for every witness.
- **G2 — repeat-call determinism at one fixed ambient setting (regression floor).** Ten fresh calls
  to `_driver_orbitals(mh)` (same fixed `mh` as G1, `a = 1.35, conv_tol = 1e-9` witness, ambient
  threads fixed at 2 to match §1's captured table), each its own subprocess, are pairwise
  bit-identical.
- **G3 — same-process repeat-call determinism, thread count restored.** Ten fresh, uncached calls to
  `_driver_orbitals(mh)` inside ONE process (no subprocess spawn) are pairwise bit-identical, and
  leave the process's ambient `pyscf.lib.num_threads()` unchanged after every call — the
  `finally`-restore is load-bearing, same as chem-s6c's G3.
- **G4 — the retry loop stops firing when the only reason it fired was thread noise.** For the fixed
  `mh` of G1, across all four ambient thread counts, `_driver_orbitals` makes exactly **one**
  `PySCFDriver.run_pyscf` call (monkeypatched call-counter) — i.e. the retry loop never actually
  retries for any of this module's own gated witnesses once its own RHF solve is pinned. This is the
  explicit before/after retry-count measurement chem-s6c's notes asked for (§1's table is the
  "before"; this gate is the "after").
- **G5 — no regression on `SPEC_rebuild_gauge_search.md`'s own mechanics.** The existing
  `tests/test_rebuild_gauge_search_spec.py` (G3-G5: cheap-when-unneeded call counting, the energy
  pre-filter's retry-not-false-accept behavior, and the oversized-active-space refusal) all still
  pass unmodified — the energy-mismatch and oversized-active-space paths are exercised by
  monkeypatching `e_tot`/the module's own bound constants, neither of which the thread pin touches.

> G1 and G4 are the definition of "done": G1 closes the determinism gap chem-s6c's notes named, G4
> is the concrete retry-count answer the bead asked for.

## 6. Implementation plan (test-first)

1. `tests/test_driver_orbitals_scf_determinism_spec.py` encoding G1-G5 — G1/G4 RED before the fix
   (verified manually, §1's captured table and §9, by reverting the `lib.num_threads` pin via
   `git stash` and rerunning).
2. `reachability.py`: pin `lib.num_threads(1)` around `_driver_orbitals`'s one `drv.run_pyscf()`
   call, restore the ambient count in a `finally`.
3. Green (§9). Confirm `tests/test_rebuild_gauge_search_spec.py` (G5) and
   `tests/test_eigenstate_reachability_spec.py` (regression floor, not re-pinned here) still behave
   the same or better.

## 7. Out of scope

- **`build_molecular_hamiltonian`'s own (first, non-rebuilt) free RHF call** landing on a different
  branch run-to-run — confirmed still present and unmerged on this branch (chem-7cx, a sibling
  sandbox branch). §2/R2 measures its effect on `_driver_orbitals`'s *unisolated* behavior directly
  rather than hand-waving past it, but fixing it is not this spec's job.
- Re-deriving whether the now-deterministic rebuild branch is "physically correct" relative to
  `specs/SPEC_eigenstate_reachability.md`'s documented witness values, or clearing
  `test_eigenstate_reachability_spec.py::test_G2[1.19]`'s flakiness end-to-end — that gate's
  remaining flakiness is a *combination* of this mechanism and chem-7cx's, and clears only once both
  land on the same branch (chem-7cx's territory per §2).
- Raising `_MAX_SIGN_SEARCH_ORBITALS`/`_MAX_REBUILD_RETRIES`, or re-opening
  `SPEC_rebuild_gauge_search.md`'s own retry/search mechanics — untouched here (G5 pins that they
  still pass as-is).
- Cross-host / cross-BLAS-backend verification (same caveat `SPEC_reachability_scf_determinism.md`
  R3 records).

## 8. Caveats and risks

- **R1:** G1-G4 deliberately build `mh` under `pyscf.lib.num_threads(1)` *in the test itself*, not
  via a production-code fix to `build_molecular_hamiltonian` (that fix is chem-7cx's, on a sibling
  branch not merged here). This isolates "does `_driver_orbitals`'s own RHF call add
  nondeterminism on top of a fixed comparison target" (chem-50w's actual question) from "is the
  comparison target itself reproducible" (chem-7cx's question) — conflating the two would make a
  pass here look like more than it is. §9 R2 shows the honest, unisolated picture: with `mh` built
  the ordinary way (ambient threads, no pin) and only `_driver_orbitals`'s own call pinned,
  `a = 1.19, conv_tol = 1e-6` still flips found across ambient thread counts `{1, 2, 4, 8}` because
  the SCF branch `mh` itself landed on also depends on those same ambient threads. That is not a
  failure of this fix — a deterministically-pinned rebuild of the *wrong* branch still, correctly,
  deterministically fails to match — it is the visible edge of chem-7cx's separate, still-open gap.
- **R2:** verified only on this container (x86_64, 8 cores, OpenBLAS, pyscf 2.13.1). Not tested on a
  second machine.
- **R3:** this fix does not and cannot make `_driver_orbitals` succeed on every input — an active
  space larger than `_MAX_SIGN_SEARCH_ORBITALS`, or a genuinely different (non-degenerate) SCF
  solution that persists across all `_MAX_REBUILD_RETRIES` attempts, still correctly returns `None`
  (G5 confirms this machinery is untouched).

## 9. Results

**Linux x86_64 container, 8 cores, pyscf 2.13.1, qiskit-nature, 2026-10-14.** `uname -a`:
`Linux 5404a8aa3661 7.0.14-linuxkit #1 SMP PREEMPT Fri Sep 18 10:19:32 UTC 2026 x86_64 GNU/Linux`.
Wall time: see handoff (`sandbox-handoffs/chem-50w.md`) for the verbatim gate-run lines and timing.

- **Before the fix** (`git stash` on just `reachability.py`'s pin), isolated `mh` (built once under
  `num_threads=1`), ambient threads varied `{1, 2, 4, 8}` at call time, `a = 1.19, conv_tol = 1e-6`:
  found = `{yes, NO, yes, NO}` — §1's second claim.
- **Before the fix**, same isolated `mh`, ambient threads held FIXED at 2, ten successive calls,
  `a = 1.35, conv_tol = 1e-9`: found = `{yes,yes,yes,yes,NO,NO,yes,yes,yes,yes}`,
  `run_pyscf` call counts = `{3,1,3,1,5,5,2,2,5,5}` — §1's table, captured directly.
- **After the fix**, identical isolated-`mh` setup, ambient threads `{1, 2, 4, 8}`,
  `a = 1.19, conv_tol = 1e-6`: found = `yes` every time, `run_pyscf` calls = 1 every time (G1/G4).
- **After the fix**, ambient fixed at 2, ten successive calls, `a = 1.35, conv_tol = 1e-9`:
  found = `yes` every time, `run_pyscf` calls = 1 every time (G2).
- **R2 (unisolated, honest picture — not a gate, recorded here as the explicit scope boundary)**:
  `mh` built the ORDINARY way (`OMP_NUM_THREADS` set before the whole process starts, no thread pin
  around the build itself, matching every real call site on this branch today), only
  `_driver_orbitals`'s own call pinned, `a = 1.19, conv_tol = 1e-6`: `OMP_NUM_THREADS` in
  `{1, 2, 4, 8}` gave found = `{yes, yes, NO, NO}` with `run_pyscf` call counts
  `{1, 1, 5, 5}` — i.e. still thread-count-dependent, because `build_molecular_hamiltonian`'s own
  unpinned free RHF call (chem-7cx, unmerged here) selects a different branch per ambient thread
  count before `_driver_orbitals` ever runs, and a deterministically-pinned rebuild of a branch that
  does not match that choice still, correctly, deterministically fails. This is the gap named in
  §2/§7/R1, not a regression.
- `tests/test_driver_orbitals_scf_determinism_spec.py -q`: **G1-G5 pass** — see handoff for the
  verbatim line.
- `tests/test_rebuild_gauge_search_spec.py -q`: unchanged, still green (G5 regression check) — see
  handoff.

## 10. Deliverables

- `reachability.py` — `_driver_orbitals` pins/restores `pyscf.lib.num_threads` internally around its
  one `drv.run_pyscf()` call.
- `tests/test_driver_orbitals_scf_determinism_spec.py` — G1-G5.
- This spec documents, rather than closes, the residual unisolated gap that is chem-7cx's to fix.
