# Handoff: chem-5fl

Docstring `Run: python <script>.py` headers across top-level scripts contradicted
CLAUDE.md's "never invoke a bare `python`/`python3`, use `uv run`" rule.

## What I changed

9 files' module docstrings: `Run:  python ...` -> `Run:  uv run python ...` (only the
`Run:` line touched in each; spacing/formatting otherwise untouched):

- `benchmark_krylov.py`
- `benchmark_n2.py`
- `benchmark_resources.py`
- `study_be2.py`
- `benchmark_nbn.py`
- `benchmark_hchain_tdl.py`
- `benchmark_dmrg.py`
- `benchmark_dmrg_large.py`
- `certchem/examples/n2_quickstart.py`

`git diff --stat`:
```
 benchmark_dmrg.py                  | 2 +-
 benchmark_dmrg_large.py            | 2 +-
 benchmark_hchain_tdl.py            | 2 +-
 benchmark_krylov.py                | 2 +-
 benchmark_n2.py                    | 2 +-
 benchmark_nbn.py                   | 2 +-
 benchmark_resources.py             | 2 +-
 certchem/examples/n2_quickstart.py | 2 +-
 study_be2.py                       | 2 +-
 9 files changed, 9 insertions(+), 9 deletions(-)
```

7 of the 16 files named in the bead description were **already compliant** when I
started (a prior session must have fixed them, presumably during chem-5iq's adjacent
pip/conda audit): `hchain_tdl_analysis.py`, `probe_hchain_localize.py`,
`hubbard_tdl_analysis.py`, `certkit_bridge.py`, `benchmark_hubbard_lieb_wu.py`,
`scripts/spec_pm3_subspace_eta_bound.py`, `scripts/spec_qpe_precision_bound.py` — I
left these untouched (already correct, no diff).

## Acceptance criteria check

Criterion: `grep -rn '^Run:.*python ' --include='*.py' . ` returns nothing outside
quarantined/historical files, where every match is `uv run python` not bare `python`.

Command run after the edits:
```
$ grep -rn '^Run:.*python ' --include='*.py' . | grep -v '.claude/worktrees' | grep -v 'uv run python'
(no output)
```
PASS — every remaining `Run:` header line in the tree (excluding the stray
`.claude/worktrees/agent-a854e794faebb6afe/` copies, which are a separate worktree's
files, not this repo's tracked tree) now reads `uv run python`.

Full list of all `Run:` lines post-fix, for review (all contain `uv run python`):
```
$ grep -rn '^Run:.*python ' --include='*.py' . | grep -v '.claude/worktrees'
benchmark_resources.py:11:Run:  uv run python benchmark_resources.py   ->  prints a resource table and writes data/krylov_resources.csv
benchmark_n2.py:14:Run:  uv run python benchmark_n2.py   ->  prints a table and writes data/n2_dissociation.csv
certkit_bridge.py:12:Run:  uv pip install -e ".[certkit]"  then  uv run python certkit_bridge.py [H2|H4|H4-stretched|N2]
study_be2.py:20:Run:  uv run python study_be2.py   ->  prints the curve + binding energies, writes data/be2_curve.csv
hchain_tdl_analysis.py:17:Run:  uv run python hchain_tdl_analysis.py [table.csv]   (default: the vendored spec table)
benchmark_nbn.py:25:Run:  uv run python benchmark_nbn.py   ->  prints a table and writes data/nbn_active_space.csv
benchmark_krylov.py:9:Run:  uv run python benchmark_krylov.py   ->  prints a table and writes data/krylov_benchmark.csv
benchmark_dmrg.py:26:Run:  uv run python benchmark_dmrg.py   ->  prints a table and writes data/dmrg_ladder.csv
probe_hchain_localize.py:11:Run:  uv run python probe_hchain_localize.py [--n 20] [--D 400] [--R 1.8] [--cap-min 90]
benchmark_hubbard_lieb_wu.py:16:Run:  nohup uv run python benchmark_hubbard_lieb_wu.py --route open --U 4 --L 20,40,60,80,100 \
hubbard_tdl_analysis.py:26:Run:  uv run python hubbard_tdl_analysis.py [table.csv]   (default: the vendored spec table)
benchmark_hchain_tdl.py:18:Run:  uv run python benchmark_hchain_tdl.py [--threads N] [--localize] [--R bohr]
benchmark_dmrg_large.py:19:Run:  uv run python benchmark_dmrg_large.py [--threads N]   ->  table + data/dmrg_large_hchain.csv
certchem/examples/n2_quickstart.py:3:Run:  uv run python -m certchem.examples.n2_quickstart
scripts/spec_pm3_subspace_eta_bound.py:30:Run: uv run python scripts/spec_pm3_subspace_eta_bound.py
scripts/spec_qpe_precision_bound.py:10:Run: uv run python scripts/spec_qpe_precision_bound.py
```

## Validation

- `uv run python -c "import ast; ast.parse(open(f).read())"` over all 9 edited files —
  all parse cleanly (no syntax breakage from the edit).
- `uv sync --extra dmrg --extra test --extra dev` (ruff wasn't installed at session
  start; synced per CLAUDE.md).
- `uv run ruff check .` over the whole repo: **All checks passed!** (docstring-only
  change, as expected, but ran the full repo lint rather than just the touched files).
- Did not run `make gates` / `make test` — this change touches only comment/docstring
  text, no importable code, control flow, or test-visible behavior. Grepped
  `tests/` and `specs/` for any reference to the old `Run: python ...` strings
  (`grep -rn "python benchmark\|python study_be2\|python certchem" tests/ specs/`) —
  no hits, so no test depends on the old wording. I did not execute any of the 9
  benchmark/study scripts themselves (several require DMRG/pyscf active-space runs
  that are multi-minute+); the change is textual only and carries no risk to their
  runtime behavior.

## What I decided not to do, and why

- Did **not** touch the 7 already-compliant files in the bead's list — editing them
  would have been a no-op diff; confirmed via the grep above that they already say
  `uv run python`.
- Did **not** fix two unrelated bare-`python` strings I found while grepping:
  `sweep_hybrid_solver.py:240` (`print("Usage: python sweep_hybrid_solver.py
  [1|2]")`) and `hybrid_quantum_solver/chemistry_gateway.py:135` (`print("Usage:
  python -m hybrid_quantum_solver.chemistry_gateway ...")`), plus two inline
  example lines in `sweep_hybrid_solver.py`'s docstring body (not under a `Run:`
  label). These match the same underlying problem (bare `python` contradicting the
  uv-run rule) but are a different pattern (`Usage:` runtime print statements, not
  `Run:` docstring headers) and weren't in chem-5fl's named file list or its
  acceptance-criteria grep. Filed as **chem-a1t** rather than doing it under this
  bead's scope.
- Did **not** touch the `.claude/worktrees/agent-a854e794faebb6afe/` copies —
  that's a separate git worktree, not this repo's tracked tree; files there are
  outside scope and will be reviewed/cleaned up independently.

## What I could not verify

- I did not execute the 9 changed scripts end-to-end to confirm the (pre-existing,
  unchanged) functional behavior still matches their docstrings — out of scope for
  a comment-only fix, and several require significant compute (DMRG ladders,
  active-space CASCI sweeps) per this repo's compute-envelope rules. The `ast.parse`
  check above confirms the edit didn't break syntax; I'm relying on that plus the
  minimal, mechanical nature of the diff (single string substitution inside a
  docstring) rather than a full run.

## Git

Tree is dirty and ready to commit (not committed — per this session's git policy,
the host commits/pushes). Suggested commands for the human:

```bash
git add benchmark_krylov.py benchmark_n2.py benchmark_resources.py study_be2.py \
        benchmark_nbn.py benchmark_hchain_tdl.py benchmark_dmrg.py \
        benchmark_dmrg_large.py certchem/examples/n2_quickstart.py \
        sandbox-handoffs/chem-5fl.md
git commit -m "chem-5fl: fix bare-python 'Run:' docstring headers to 'uv run python'

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

## Beads

- `chem-5fl`: closing as done, criteria met (see grep output above).
- `chem-a1t`: new bead filed for the two `Usage:` print-string occurrences +
  docstring-body example lines found during this audit, deliberately left out of
  this bead's scope.
