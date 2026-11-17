# chem-66q handoff

## Disposition: closed, implemented as specified

chem-7ko's n=10/12 `PQG+T2sym` runs (SPEC_hchain_t2sym) solved for hours and never saved the SCS
dual vector, so re-checking a certified bound later (e.g. after any change to
`rigorous_lower_bound`'s rounding) meant a 2-4h re-solve just to get `y` back. Separately, a hard
solver failure surfaced three steps downstream as `TypeError: unsupported operand type(s) for +:
'NoneType' and 'float'` at `val + e_nuc` in `main()`, instead of naming what actually failed.

## What changed

`vsdp_hchain.py`:
- **`require_solved(val, y, status, label)`** (new) — rejects `val is None`, `y is None`, a
  non-finite `val` (inf/nan), or a status outside `{optimal, optimal_inaccurate, solved,
  solved_inaccurate}`, raising a labeled `RuntimeError` naming the status and value. Wired into
  both `solve()` (cvxpy/CLARABEL path) and `solve_scs()` (the direct-SCS path used for `--sym`),
  so a hard failure in either backend is caught at the point of failure, not at `val + e_nuc`
  several lines later in `main()`.
- **`recertify(n, t2, dual_y, e_nuc)`** (new) — rebuilds the SDP problem data via `build()`
  (seconds) and calls `certify()` on a previously-solved, already-serialized dual vector, without
  calling `solve()`/`solve_scs()` again. Checks the supplied `e_nuc` against the rebuilt
  integrals' `e_nuc` and raises `ValueError` on mismatch (catches a dual paired with the wrong
  `(n, t2)` rather than silently recertifying garbage).
- **`main()`**: every freshly-solved row now also stores `e_nuc` and `dual_y` (`y.tolist()`, a
  plain JSON-round-trippable list). New `--recertify` flag: for each `n`, loads the stored row's
  `dual_y`/`e_nuc` from `--output` and calls `recertify()` instead of re-solving; raises a clear
  `RuntimeError` naming the file and `n` if a row predates this change and has no `dual_y`.

`specs/SPEC_hchain_vsdp_bracket.md`: documented the new public interface (`recertify`,
`--recertify`), added gates G6/G7 to §5, and a dated note in §9 making explicit that this does
**not** change any number in the §6 results table — only how future rows are saved/re-checked.

`tests/test_hchain_vsdp_bracket_spec.py`: added
- **G6** `test_g6_recertify_from_stored_dual_without_resolving` — solves n=4 once, round-trips `y`
  through `.tolist()` (simulating the JSON write/read), calls `recertify(4, False, stored_y,
  e_nuc)` (which internally rebuilds the problem but never calls `solve()` again), and asserts the
  recertified bound exactly equals what `certify()` gave at solve time, and is `<=` FCI.
- **G6b** `test_g6b_recertify_rejects_mismatched_e_nuc` — a dual stored against the wrong `e_nuc`
  must raise `ValueError`, not silently recertify.
- **G7** `test_g7_hard_solver_failure_raises_clearly` — directly exercises `require_solved` with
  `None`, non-finite, and missing-dual inputs (all must raise `RuntimeError` with a readable
  message) and confirms a genuinely-solved input passes through unchanged.

## Pre-registered criteria and outcome

| Criterion (from the bead) | Outcome |
|---|---|
| Save the dual next to each row | **Done** — `dual_y` (+ `e_nuc`) on every freshly-solved row in `main()`. |
| Re-certifying later must not need a 2-4h re-solve | **Done** — `recertify()` / `--recertify` redo only `build()`; demonstrated live below (n=2, CLARABEL) and gated at n=4 (G6). |
| A hard cvxpy failure gets a clear error, not an unhelpful TypeError at `val + e_nuc` | **Done** — `require_solved` raises a labeled `RuntimeError` at the point of failure in both `solve()` and `solve_scs()`; gated by G7. |
| Add a gate that re-certifies a stored dual against FCI | **Done** — G6 (n=4, below the FCI ceiling it's checked against). |

## Commands and their output

Gate file, isolated process (per CLAUDE.md's block2/pyscf isolation rule — this file imports
pyscf but not block2, still run standalone via the harness):

```
$ GATE_RUN="uv run" GATE_JOBS=1 GATE_GLOB='tests/test_hchain_vsdp_bracket*_spec.py' bash scripts/run_gates.sh
gates: 1 files, 1 parallel processes, cache=on
PASS  tests/test_hchain_vsdp_bracket_spec.py  (86s)
all spec gates passed
```

Same file via plain pytest (10 cases: G1 x2 params, G2, G3 x2 params, G4, G5, G6, G6b, G7):

```
$ uv run pytest -q tests/test_hchain_vsdp_bracket_spec.py
..........                                                               [100%]
10 passed in 85.26s (0:01:25)
```

Lint:

```
$ uv run ruff check .
All checks passed!
```

Live smoke test of the actual CLI path (n=2, CLARABEL, ~3s) confirming the dual is saved and
`--recertify` reproduces the identical `certified_lb` from the JSON file alone:

```
$ uv run python vsdp_hchain.py --ns 2 --solver CLARABEL --output /tmp/vsdp_smoke.json
{"n": 2, ..., "certified_lb": -1.118733937010524, ..., "dual_y": [...88 floats...]}

$ uv run python vsdp_hchain.py --ns 2 --solver CLARABEL --output /tmp/vsdp_smoke.json --recertify
{"n": 2, ..., "certified_lb": -1.118733937010524, ...}   # identical, no re-solve
```

And the clear-error path on a row stripped of `dual_y` (simulating a pre-chem-66q row):

```
$ uv run python vsdp_hchain.py --ns 2 --solver CLARABEL --output /tmp/vsdp_smoke_nodal.json --recertify
RuntimeError: n=2: no stored dual in /tmp/vsdp_smoke_nodal.json to recertify (row predates chem-66q;
re-run a full solve once to get one)
```

(Both temp files were scratch and are not part of the diff; `data/` is gitignored anyway.)

## Hardware / wall time

`uname -a`: `Linux 1ce41d8d7c38 7.0.14-linuxkit #1 SMP PREEMPT ... x86_64 GNU/Linux`, `nproc`: 8.
The 85s pytest run and 86s gate run above are the only commands that took over a minute; both are
pre-existing SDP solves at n=2/4/6 that G2/G3/G5/G6 already did before this bead (this bead adds no
new expensive solve — n stays small in all new/changed test code, by design, so the point of the
feature — skipping an hours-long re-solve — doesn't itself require running one in the gate).

## What I decided not to do, and why

- **Did not touch `--table`'s existing semantics.** It already re-pairs a stored row with a fresh
  upper bound without touching the lower bound; conflating it with recertification would change
  existing behavior silently. `--recertify` is a separate, additive flag instead.
- **Did not re-run the actual n=10/12 `--sym` solves** to populate `dual_y` retroactively for
  existing `data/hchain_vsdp_t2sym.json` rows (if a worker's local copy has one) — `data/` is
  gitignored and not present in this checkout, and re-solving n=12 is exactly the 2-4h cost this
  bead exists to avoid paying for a bookkeeping change. The spec's §9 note says existing rows have
  no `dual_y` and `--recertify` reports that honestly rather than guessing or re-solving silently.
- **Did not change `rigorous_lower_bound`/`certify`'s numerics.** Out of scope — this bead is
  about persistence and error clarity, not the bound's math. No number in SPEC_hchain_vsdp_bracket
  §6 changed.
- **Did not try to force a real external solver failure** (e.g. feeding CLARABEL/SCS adversarial
  data to trigger `SOLVER_ERROR`/`infeasible` for real) for G7. Traced cvxpy's and SCS's actual
  failure-status code paths first (`cvxpy/reductions/solution.py`'s `failure_solution`, which does
  set `opt_val = None` for `SOLVER_ERROR`; SCS's `info['status']` for non-`solved` cases) to
  confirm `require_solved`'s checks match real failure shapes, then unit-tested `require_solved`
  directly with those shapes — deterministic and fast, versus trying to coerce a real solver into
  a flaky failure mode.

## What could not be verified

- Whether a real n=10/12 `PQG+T2sym` solve's dual vector (thousands of floats, SU(2)×reflection
  adapted blocks) round-trips through `json.dump`/`json.load` without precision loss at the scale
  `rigorous_lower_bound`'s rounding bookkeeping needs — not re-run here (hours, see above). CPython
  `json` uses the shortest round-trip `repr` for floats since 3.1, which is lossless in principle;
  the n=2 smoke test above exercises the exact same code path (`y.tolist()` -> `json.dump` ->
  `json.load` -> `np.asarray(..., float)` -> `certify`) and reproduces `certified_lb` bit-for-bit,
  which is strong evidence but not a direct n=12-scale confirmation.

## Files changed

```
$ git diff --stat
 specs/SPEC_hchain_vsdp_bracket.md      | 22 ++++++++++++++--
 tests/test_hchain_vsdp_bracket_spec.py | 46 ++++++++++++++++++++++++++++++++
 vsdp_hchain.py                         | 48 +++++++++++++++++++++++++++++++---
 3 files changed, 112 insertions(+), 4 deletions(-)
```

Plus `.beads/issues.jsonl` (bead claim/close export) and this handoff file.

## Suggested commands for the human

```bash
git add vsdp_hchain.py specs/SPEC_hchain_vsdp_bracket.md tests/test_hchain_vsdp_bracket_spec.py \
        sandbox-handoffs/chem-66q.md .beads/issues.jsonl
git commit  # see the actual commit(s) made this session for message text
git status
```

Not pushed, merged, or rebased by this session — work is committed to `sandbox/chem-66q` here, the
host pushes and a human merges, per this run's git policy.
