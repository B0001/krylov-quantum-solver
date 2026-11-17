# chem-a1t handoff

## What changed

Bare `python` in runtime/example invocation strings (outside chem-5fl's "Run:" docstring-header
scope) changed to `uv run python`, per CLAUDE.md's uv-run rule.

The bead named four occurrences; I fixed those plus a re-grep turned up three more genuine
example/Usage-style mentions in the same spirit. All are doc/comment/print-string text only — no
executable logic changed.

1. `sweep_hybrid_solver.py:240` — `print("Usage: python sweep_hybrid_solver.py [1|2]")` ->
   `print("Usage: uv run python sweep_hybrid_solver.py [1|2]")`
2. `sweep_hybrid_solver.py:7,9` — module-docstring inline examples (`python sweep_hybrid_solver.py
   1` / `2`, not under a `Run:` label) -> `uv run python sweep_hybrid_solver.py 1` / `2`
3. `hybrid_quantum_solver/chemistry_gateway.py:135` —
   `print("Usage: python -m hybrid_quantum_solver.chemistry_gateway ...")` -> `uv run python -m ...`
4. `run_senseforge_sweep.py:5` — docstring inline example `python run_senseforge_sweep.py ...` ->
   `uv run python run_senseforge_sweep.py ...`
5. `run_real_material_simulation.py:5` — docstring inline example -> `uv run python
   run_real_material_simulation.py ...`
6. `nb3x8_gaps.py:54` — prose example `` ``python nb3x8_gaps.py`` `` -> `` ``uv run python
   nb3x8_gaps.py`` ``
7. `compile_from_gcs.py:22` — `Usage:\n    python compile_from_gcs.py [--bucket NAME] [--verify]`
   -> `uv run python compile_from_gcs.py ...`
8. `tests/test_senseforge_spec.py:301` — reproduction comment `# Produced by \`python
   run_senseforge_sweep.py ...\`.` -> `uv run python run_senseforge_sweep.py ...`

## Re-grep (command run, as the bead specifies)

```
grep -rn 'python ' --include='*.py' . | grep -v '.claude/worktrees' | grep -v 'uv run python\|__pycache__' \
  | grep -iv 'python3\.\|pythonpath\|sys\.executable\|python_requires\|python_version\|language'
```

Remaining hits after my fixes, all correctly out of scope:

- `benchmark_krylov.py:9`, `benchmark_nbn.py:25`, `benchmark_dmrg_large.py:19`,
  `study_be2.py:20`, `benchmark_resources.py:11`, `benchmark_n2.py:14`, `benchmark_dmrg.py:26`,
  `benchmark_hchain_tdl.py:18`, `certchem/examples/n2_quickstart.py:3` — all `^Run:` docstring
  headers. That's chem-5fl's territory (bead text: "outside chem-5fl's scope (which covered 'Run:'
  docstring headers only)"); several of these still read bare because chem-5fl's fix commits
  (`11987ef`, `3873ed4`) live on a different sandbox branch not yet merged into this one — not a
  regression I introduced, and not this bead's job to fix.
- `tests/test_scf_conv_tol_spec.py:9` — already `uv run --no-sync python -c ...`, fine as-is.
- `tests/test_gate_runner.py:43` — `# \`$RUN python -m pytest\` with RUN="" needs \`python\` to be
  this interpreter.` This is explaining PATH/interpreter-substitution semantics for a test
  fixture's subprocess env, not an instruction to copy-paste and run `python ...` — left alone.
- `.claude/worktrees/**` — stale leftover worktree directories from prior sessions, not tracked by
  git (`git ls-files .claude/worktrees` is empty); not part of the repo.

I did not touch the `.claude/worktrees/` copies or the `Run:` headers — both are other beads'
territory or stale untracked scratch dirs, not this bead's scope.

## Validation

- `uv run ruff check .` -> `All checks passed!`
- `uv run pytest tests/test_senseforge_spec.py -q` -> `30 passed, 72 warnings in 4.42s` (the one
  spec file whose test source I edited — a comment-only change, confirmed still green)
- `uv run python -m py_compile` on all seven edited `.py` files -> exits clean

No gate run across the full `tests/test_*_spec.py` set was needed: every edit is a string literal
inside a docstring, a comment, or a `print()` argument — zero control-flow or assertion changes, so
nothing beyond the touched file's own tests could regress. Did not run `make gates` / `make test`
(expensive DMRG/multi-process suite) since nothing executable changed; flagging as unverified if a
reviewer wants the full-suite confirmation.

## Pre-registered criteria (from the bead)

All four named instances fixed, plus the explicitly requested re-grep sweep performed and its three
extra genuine hits (`run_senseforge_sweep.py`, `run_real_material_simulation.py`, `nb3x8_gaps.py`,
`compile_from_gcs.py` — four, not three, counting both docstring-example files) also fixed.
`tests/test_senseforge_spec.py`'s comment fixed for consistency since it's a copy-paste-style
reproduction instruction, same category as the others.

## What I decided not to do

- Did not touch the `Run:`-header lines (`benchmark_*.py`, `study_be2.py`,
  `certchem/examples/n2_quickstart.py`) — explicitly chem-5fl's scope per the bead text, even
  though some still read bare on this branch (see re-grep section above for why).
- Did not touch `tests/test_gate_runner.py:43` — it's an explanatory comment about subprocess PATH
  substitution, not a runnable example a user would copy.
- Did not clean up `.claude/worktrees/` — untracked scratch directories, out of scope for a
  docstring-wording bead.

## What I could not verify

- Whether chem-5fl's branch will merge cleanly with this one given the overlapping-but-disjoint
  `Run:`-header edits — not this bead's concern, but worth the merging human knowing both branches
  touch nearby doc lines in the same files.
- Full `make gates` / `make test` run not executed (see Validation) — confident it's unaffected
  given the change is string-literal-only, but "confident" is not "measured".

## Git status (uncommitted, per the sandbox git policy for this task)

```
$ git status --porcelain
 M .beads/issues.jsonl
 M compile_from_gcs.py
 M hybrid_quantum_solver/chemistry_gateway.py
 M nb3x8_gaps.py
 M run_real_material_simulation.py
 M run_senseforge_sweep.py
 M sweep_hybrid_solver.py
 M tests/test_senseforge_spec.py
?? .claude/settings.local.json
?? sandbox-handoffs/chem-a1t.md
```

Suggested commit (human runs, per this task's git-on-this-run instructions):

```
git add compile_from_gcs.py hybrid_quantum_solver/chemistry_gateway.py nb3x8_gaps.py \
  run_real_material_simulation.py run_senseforge_sweep.py sweep_hybrid_solver.py \
  tests/test_senseforge_spec.py sandbox-handoffs/chem-a1t.md .beads/issues.jsonl
git commit -m "chem-a1t: uv-run-ify bare 'python' in Usage:/example strings outside chem-5fl's Run: scope"
```
