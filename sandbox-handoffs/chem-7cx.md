# chem-7cx handoff

## Disposition: fixed and gated, with two narrower follow-ups filed

`build_molecular_hamiltonian`'s free RHF is now deterministic: repeated builds of the same
`(atom, basis, charge, spin, conv_tol)` give a bit-identical `qubit_hamiltonian` and `hf_energy`,
on this host, regardless of ambient thread count. The bead's two named witnesses
(`test_G1_residue_level_identified_warned_and_kept[1.1-1e-06-1e-08]`,
`test_G9_control_tight_scf_stays_on_the_ag_level[1.1]`) no longer flip pass/fail across repeated
full-suite runs — they now fail the *same way* every time. Fixing reproducibility is not the same
claim as fixing which branch is "platform-correct" relative to the spec's macOS-recorded witness
values; that remains chem-g14's separate, pre-existing, still-open scope (see below).

## Root cause

PySCF's own OpenMP-threaded integral evaluation is not bit-associative. On a near-degenerate
geometry (square H4, sto-3g), that noise is occasionally enough to tip the free
(`symmetry=False`) RHF solve onto a different locally-stable solution between otherwise-identical
calls. Isolated by toggling one env var at a time:

```
$ OMP_NUM_THREADS=1 <10x repro>        # deterministic (fixes it)
$ OPENBLAS_NUM_THREADS=1 <10x repro>   # still flips (does not fix it)
```

This places the noise in PySCF's own compiled integral-evaluation extensions, not in
numpy/BLAS linear algebra — confirming the hypothesis in the bead's own text.

## Fix

`hybrid_quantum_solver/molecular_hamiltonian.py`: pin `pyscf.lib.num_threads(1)` around the
`PySCFDriver(...).run()` call, restore the ambient count in a `finally`. Identical pattern to the
already-validated `benchmark_nbn._tight_scf` fix (chem-y6w, `specs/SPEC_nbn_scf_determinism.md`).
No signature change, no new parameters, zero blast radius for callers:

```diff
-    problem = PySCFDriver(atom=atom, basis=basis, charge=charge, spin=spin,
-                          conv_tol=effective_conv_tol).run()
+    ambient_threads = lib.num_threads()
+    lib.num_threads(1)
+    try:
+        problem = PySCFDriver(atom=atom, basis=basis, charge=charge, spin=spin,
+                              conv_tol=effective_conv_tol).run()
+    finally:
+        lib.num_threads(ambient_threads)
```

Also added a `Note:` to the function's docstring (see the diff in `git show` for the full text)
documenting the mechanism, the measurement, and the explicit non-claim about cross-platform
branch-matching.

## SDD artifacts

- `specs/SPEC_scf_branch_determinism.md` — goal, background, approach, 4 acceptance gates (G1-G4),
  out-of-scope, caveats, results. Modeled on `SPEC_nbn_scf_determinism.md`'s structure.
- `tests/test_scf_branch_determinism_spec.py` — the gate file. pyscf + qiskit only (no block2),
  own file per the `test_*_spec.py`-runs-in-its-own-process convention; G1/G2/G4 use subprocesses
  regardless.

```
$ uv run pytest tests/test_scf_branch_determinism_spec.py -v
tests/test_scf_branch_determinism_spec.py::test_G1_deterministic_across_ambient_thread_counts PASSED
tests/test_scf_branch_determinism_spec.py::test_G2_deterministic_across_repeats_separate_processes PASSED
tests/test_scf_branch_determinism_spec.py::test_G3_deterministic_across_repeats_same_process_and_threads_restored PASSED
tests/test_scf_branch_determinism_spec.py::test_G4_known_flipping_witnesses_stabilize_pass_or_fail PASSED
4 passed, 10 warnings in 59.87s
```

- G1: fresh builds at ambient `OMP_NUM_THREADS` 2 and 8 → bit-identical digest.
- G2: 10 fresh builds, separate processes, default ambient threads → pairwise bit-identical
  (directly reproduces the bead's own original 10-trial experiment, which found 9/10 vs 1/10；
  now 10/10).
- G3: 10 builds in ONE process (the `_system()`-`lru_cache` masked case) → bit-identical, and
  `lib.num_threads()` is restored to exactly what this process found it as after every call.
- G4: the bead's two named witness tests, run 5x in separate pytest subprocess invocations → same
  outcome every time (does not assert which outcome — that's chem-g14's territory).

## Regression check: the full `test_eigenstate_reachability_spec.py` file

Before the fix (`git stash`, 2 full-file runs): different test IDs failed run to run — e.g.
`test_G9_control_tight_scf_stays_on_the_ag_level[1.19]` failed in run 1 but passed in run 2, while
`test_G2_decision_equals_symmetry_adapted_fci[...1.19...-1e-09]` did the opposite. 14 failures
each run, not the same 14.

After the fix, 3 full-file runs:

```
run 1: 17 failed, 32 passed, 145.42s
run 2: 17 failed, 32 passed, 128.97s
run 3: 16 failed, 33 passed, 133.09s
```

16 of these failing tests are now **perfectly stable** (identical outcome all 3 runs), including
both of the bead's own named witnesses. Only one pair still flips:
`test_G2_decision_equals_symmetry_adapted_fci` at a=1.19 (both `conv_tol` parametrizations) —
run 1 failed the `1e-09` variant, run 2 failed the `1e-06` variant, run 3 failed neither.

Investigated this residual flake directly rather than assume my fix was incomplete:

```
$ uv run python /tmp/repro_a119.py   # build_molecular_hamiltonian itself, a=1.19, both conv_tol, 4x
# bit-identical digest every trial — build_molecular_hamiltonian is NOT the source here.

$ uv run python /tmp/repro_fci.py    # tests/test_eigenstate_reachability_spec.py's own _fci()
                                      # helper's raw scf.RHF(mol, symmetry="D2h") call, unpinned, 6x
-1.7071494724561918 / ...926 / ...926 / ...926 / ...918 / ...926   # ~1e-13 bit noise

$ uv run python /tmp/repro_fci_pinned.py   # same call with lib.num_threads(1) pinned, 6x
-1.7071494724561926 (identical all 6 trials)
```

`_fci()` (`tests/test_eigenstate_reachability_spec.py:66-82`) makes its own independent,
unpinned `scf.RHF` call that `build_molecular_hamiltonian` never touches. It shares chem-7cx's
exact root-cause mechanism (confirmed above: unpinned → ~1e-13 noise, pinned → bit-identical), and
at a=1.19's near-degenerate geometry that tiny noise is apparently amplified through the FCI
Davidson solve enough to occasionally flip a classification right at a threshold boundary — even
though `build_molecular_hamiltonian`'s own output for the same geometry is rock solid. **Filed as
chem-0yh** rather than fixed here: it's a third, independent unpinned call site, outside
`build_molecular_hamiltonian`, which is what chem-7cx's title and fix are about.

## Bonus: `test_scf_conv_tol_spec.py` also stabilizes (not named in the bead)

```
before fix (3 runs): 3, then 6, then 6 failures — different test IDs each time
after fix (3 runs):  4 failed, 10 passed — identical IDs, all 3 runs:
  test_G1_clean_geometry_does_not_move[H4 square 1.2]
  test_G2_witness_moves_but_the_eigenvalue_does_not[H4 square 1.35]
  test_G2_witness_moves_but_the_eigenvalue_does_not[H4 square 1.19]
  test_G3_residue_collapses_at_the_tight_tolerance[H4 square 1.35]
```

Same mechanism, same fix, affected this file too — not claimed in the bead title but confirmed
the same way (repeated before/after trials, not a single-run comparison, since the before-fix
state is itself too flaky to trust a single sample).

## Out-of-scope findings filed separately (not fixed here)

- **chem-g14** (pre-existing, open): which branch is "platform-correct" relative to the spec's
  macOS-measured witness values. The bead's own notes already record that on macOS 27, main
  passes `test_eigenstate_reachability_spec` 49/49 and `test_scf_conv_tol_spec` 14/14, while on
  this Linux container they fail 14-18 and 5-6. Determinism (this bead) and platform-correctness
  (chem-g14) are separate, honestly-separated claims — fixing the former does not fix the latter,
  and this session did not attempt to.
- **chem-s6c** (filed this session, before this handoff): `reachability.py`'s own raw `scf.RHF`
  calls (`scf_symmetry_status`, `hf_orbital_irreps`) share the identical mechanism but are a
  different function/file outside `build_molecular_hamiltonian`.
- **chem-0yh** (filed this session, see above): `tests/test_eigenstate_reachability_spec.py`'s own
  `_fci()` helper has the same unpinned-call mechanism, causing the one residual `test_G2[1.19]`
  flake documented above.

## Pre-registered acceptance criteria and outcome

| Criterion | Outcome |
|---|---|
| Root cause of the branch nondeterminism identified | **Done** — PySCF's own OpenMP-threaded integral evaluation, not BLAS/linalg (OMP vs OPENBLAS isolation test). |
| A fix implemented in `build_molecular_hamiltonian` | **Done** — `pyscf.lib.num_threads(1)` pinned/restored, matching the chem-y6w precedent. |
| New SDD spec + gate test, gates pass | **Done** — `specs/SPEC_scf_branch_determinism.md` + `tests/test_scf_branch_determinism_spec.py`, 4/4 passing (verified 3x). |
| Bead's two named witnesses stop flipping across repeated runs | **Done** — confirmed via G4 (5 subprocess trials) and via 3 full before/after file runs. |
| Broader blast radius checked, not just the two named tests | **Done** — full-file regression on `test_eigenstate_reachability_spec.py` (16/17 previously-flaky tests now stable) and `test_scf_conv_tol_spec.py` (now fully stable). |
| Out-of-scope discoveries filed, not fixed here | **Done** — chem-s6c, chem-0yh. |
| `make lint` / `uv run ruff check .` clean | **Done** — one issue found and fixed (unused `pytest` import in the new spec file), now clean. |
| Core validation gates still pass | **Done** — `test_reference_energies.py`, `test_krylov_convergence.py`, `test_noise_resilience.py`: 20 passed. |

## What was not run, and why

The repo has 115 `tests/test_*_spec.py` gate files (`make gates` runs every one, each in its own
process — `make` itself is not installed in this container, so I ran the equivalents directly via
`uv run pytest`). Running the entire 115-file suite is disproportionate for this change: the fix
is a pure internal implementation detail (thread pinning around one already-existing call, no
signature or behavior change to the branch-selection logic itself) whose only observable effect is
determinism, and its blast radius is bounded to callers of `build_molecular_hamiltonian` on
near-degenerate geometries. I instead verified: the 3 core validation gates, the new spec's own
gate (3x), and full repeated-run regression checks on the two test files most directly exercising
this code path at the flagged geometries. I did not attempt a full `make gates` sweep.

## What could not be verified

- Cross-platform behavior (macOS, different PySCF/BLAS builds) — only tested on this Linux
  container (Debian 12, `nproc`=8). The spec and this handoff are explicit that the fix is a
  same-host reproducibility guarantee, not a cross-platform result-matching guarantee.
- Whether pinning `lib.num_threads(1)` has a measurable wall-clock cost for larger active spaces
  than square H4 (sto-3g, 8 qubits) — not benchmarked; out of scope for a reproducibility bug fix.

## Files changed

```
$ git diff --stat
 hybrid_quantum_solver/molecular_hamiltonian.py | 26 +++++++++++++++++++++++---
 1 file changed, 23 insertions(+), 3 deletions(-)
```

Plus new files `specs/SPEC_scf_branch_determinism.md`, `tests/test_scf_branch_determinism_spec.py`,
this handoff, `.beads/issues.jsonl` (bead notes/close + chem-s6c/chem-0yh creation), and
`bd update chem-7cx --notes=...`. `.claude/settings.local.json` (pre-existing untracked, unrelated)
is deliberately not committed.

## Suggested commands for the human

```bash
git status
git log --oneline -5
```

Work is committed to `sandbox/chem-7cx` in this session per the session's git policy — not
pushed, merged, or rebased. The host pushes; a human merges.
