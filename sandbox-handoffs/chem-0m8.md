# chem-0m8 handoff

## Disposition: fixed (test-only change, gate-robustness). Closing.

`test_G2_decision_equals_symmetry_adapted_fci`'s square-H4 cases and
`test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged` /
`test_G4_naive_labels_would_flag_a_physical_level` hardcoded a specific RHF-branch outcome
(symmetric vs. symmetry-broken) for specific bond lengths, when chem-7cx shows that branch choice
is platform/BLAS/LAPACK-dependent, not a function of the code.

**This bead went through two passes in this container, and the first pass's conclusion was wrong.**
That first pass (notes below, superseded) decided `test_G2`'s square-H4 over-inclusion was a
*structural* fact independent of RHF branch, and shipped it as an `xfail`. A later review of that
diff on macOS 27 (Apple M3) found all 8 of those cases **`xpass`** there (`41 passed, 8 xpassed`):
the over-inclusion only happens when the C2y bit (D2h bit 0) is *excluded* (RHF breaks symmetry, or
an exactly-degenerate pair rotates arbitrarily under the untracked D4h-vs-D2h gap) — exactly the
branch-dependent mechanism `specs/SPEC_eigenstate_reachability.md` section 8 R2 (chem-3fd) already
named. The `xfail` was a wrong, overly strong claim. This handoff reflects the **final, fixed**
state, already committed on this branch (`49c9bda`, on top of `80be1d7` which merged in
chem-g14's `_driver_orbitals` fix and chem-3fd's R2 text):

1. **`test_G2`** now computes `orbital_parity_bounds(mh)` itself and checks whether the C2y bit
   (bit 0) is trusted (`beta < 1.0`) on *this* run. If trusted: assert exact FCI agreement (as
   before). If not trusted: assert the conservative direction only — every true Ag level is still
   present and never flagged forbidden, but extra non-Ag levels may also be kept (R2's documented
   over-inclusion). Either branch is now a defined pass, not an assumption.
2. **`test_G4`**'s hardcoded `a=1.20`/`a=1.40` cases call a new `_require_broken_rhf(a)` helper that
   checks `scf_symmetry_status` first and `pytest.skip`s with the measured `dE` when the free RHF
   did not actually break symmetry at that `a` on this host, instead of asserting a hardcoded
   outcome or crashing on `orbital_parity_bounds(mh) is None`.
3. R5 was added to the spec documenting the first pass's wrong claim and why — left on the record,
   not deleted, per CLAUDE.md's "never quietly delete a number/claim."
4. `specs/SPEC_eigenstate_reachability.md` section 8 R2's own text already carries the fix
   (`replacing an earlier xfail-based version — see R5 below`); R5 documents the supersession.

Host this session: Linux `4bb2b2fd0fb0` 7.0.14-linuxkit x86_64, 8 cores (`nproc`), pyscf 2.13.1,
`uv run python 3.13.15`. All commands below took well under a minute except the full-file run
(~2 minutes, noted).

## Pre-registered criteria and outcome

| Criterion | Outcome |
|---|---|
| `test_G2`'s square-H4 cases pass or cleanly branch, never assert a hardcoded RHF-branch outcome | **PASS.** Code at `tests/test_eigenstate_reachability_spec.py:132-159` computes `orbital_parity_bounds` and branches on `c2y_trusted` per-run. |
| `test_G4`'s hardcoded bond lengths skip with a clear reason when the branch premise doesn't hold, instead of crashing or asserting wrongly | **PASS.** `_require_broken_rhf` (lines 200-208) calls `scf_symmetry_status` and skips with the measured `dE` when `broken=False`. |
| `test_G2`/`test_G4` gates pass or cleanly skip on this host, repeatably | **PASS** — see gate-run lines below: 3 repeated `-k "G2 or G4"` runs, identical `12 passed, 2 skipped` every time. No xfail/xpass ambiguity remains (the superseded first pass's 8 `xfailed` are gone; these 8 cases now `PASS` outright on this host because the free RHF landed on the C2y-trusted branch for every bond length tried here this session). |
| Whole-file suite does not regress outside this bead's scope | **No new failures introduced.** Whole-file run: `4 failed, 43 passed, 2 skipped` — all 4 failures are `test_G1`/`test_G9_control`, chem-7cx's free-RHF-nondeterminism territory (already filed, open), never `test_G2`/`test_G4`. |
| `ruff check .` passes | **PASS** — `All checks passed!` |

## Gate-run lines (verbatim, this host, this session, after the fix)

```
$ uv run pytest tests/test_eigenstate_reachability_spec.py -q
4 failed, 43 passed, 2 skipped, 37 warnings in 132.69s (0:02:12)
```
The 4 failures: `test_G1_residue_level_identified_warned_and_kept[1.1-1e-06-1e-08]`,
`[1.1-1e-06-1e-10]`, `[1.19-1e-06-1e-10]`, `test_G9_control_tight_scf_stays_on_the_ag_level[1.1]` —
all chem-7cx's free-RHF basin nondeterminism, not this bead's symptom (confirmed: none of
`test_G2`/`test_G4` appear).

Isolated, repeated 3x to confirm `test_G2`/`test_G4` stability specifically:

```
$ uv run pytest tests/test_eigenstate_reachability_spec.py -k "G2 or G4" -q      # x3, identical each time
12 passed, 2 skipped, 35 deselected, 12 warnings in 20.24s
12 passed, 2 skipped, 35 deselected, 12 warnings in 20.67s
12 passed, 2 skipped, 35 deselected, 12 warnings in 20.74s
```

The 2 skips are both `test_G4`'s `a=1.20` cases
(`test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged[1.2]`,
`test_G4_naive_labels_would_flag_a_physical_level`) — this host's free RHF does not break symmetry
at `a=1.20` (consistent with chem-3fd's and this bead's own first-pass finding: 20/20 and 15/15
trials respectively never broke there). `a=1.40` (`test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged[1.4]`)
and `test_G4_lih_c2_bit_is_exact` pass outright. All 10 `test_G2` cases pass outright (the
C2y-trusted branch of the new conditional) on every one of the 3 repeated runs — i.e. this host's
free, unsymmetrized RHF solve for square H4 at 1.05/1.10/1.19/1.35 landed on a C2y-trusted MO basis
every single time this session, so the conservative/over-inclusion branch of the new assertion was
not exercised here. (It was exercised, and is the reason this bead's first pass was wrong, on the
macOS 27 host in the 2026-10-13 review noted in the bead and in spec R5: `41 passed, 8 xpassed`
pre-fix, now presumably all outright `PASS` there too under the branch-aware check, though this
session did not have access to that host to re-confirm post-fix.)

```
$ uv run ruff check .
All checks passed!
```

## What changed (already committed on this branch before this session started; this session verified it)

- `tests/test_eigenstate_reachability_spec.py`:
  - `G2_CASES` dropped the `xfail` marks entirely; it's now a plain list of `(atom, conv_tol)` pairs.
  - `test_G2_decision_equals_symmetry_adapted_fci` computes `orbital_parity_bounds(mh)`, checks
    `c2y_trusted = bounds is not None and bounds[0][1] < 1.0`, and asserts exact FCI agreement when
    trusted, or the conservative containment property (every true Ag level present, none flagged
    forbidden, possible extras) when not.
  - `_require_broken_rhf(a)` helper added; both `test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged`
    and `test_G4_naive_labels_would_flag_a_physical_level` call it before touching
    `orbital_parity_bounds`, with an explicit skip reason + measured `dE` when the premise fails.
- `specs/SPEC_eigenstate_reachability.md` section 8: R2's text updated to describe the branch-aware
  `test_G2` check; new R5 entry records the first pass's wrong `xfail`-based claim, why it was wrong
  (the macOS xpass), and that it's superseded rather than deleted.

No production code (`reachability.py`, `hybrid_quantum_solver/`) was touched in either pass — this
bead is a test/spec-only fix, per its own framing.

## What this session did (vs. what was already there)

The code and spec changes above were already committed on this branch (`49c9bda`) when this
session started — a prior worker's second pass, after the bead's own "2026-10-13 review" notes
(in `bd show chem-0m8`) identified the first pass's `xfail` as wrong and specified the fix
(drop the xfail and R5-as-written, make G2 branch on C2y-bit trust, keep the G4
`_require_broken_rhf` change). This session: re-derived nothing, trusted nothing on the reviewer's
word alone, and instead **independently re-ran the gates from a clean `uv sync`** on this Linux
container to confirm the committed fix actually holds here (it does — see gate-run lines above),
ran `ruff check .`, confirmed the spec text and test code are mutually consistent (R2/R5 describe
exactly what the committed test file does), and confirmed `chem-0m8` is not itself blocked by any
open dependency (it blocks `chem-g14`, not the reverse — `bd dep tree chem-0m8` and `bd blocked`
both checked). Updated this handoff file, which still described the superseded first pass, to match
the final state.

## What I decided not to do, and why

- **Did not fix chem-g14** (`_driver_orbitals`'s rebuild-reproducibility check) or **chem-7cx**
  (free RHF basin nondeterminism) — both already filed, both open, both out of this bead's scope
  (the bead's own title names only `test_G2`/`test_G4`'s hardcoded-branch assertions). They are the
  reason `test_G1`/`test_G9_control` still fail on this host; confirmed those 4 failures are
  attributable to chem-7cx by name (free-RHF basin selection), not a chem-0m8 symptom, by checking
  that none of the failing test IDs are `test_G2`/`test_G4`.
- **Did not run the full `make gates`/90-spec suite.** This bead's diff touches exactly one file
  (`tests/test_eigenstate_reachability_spec.py`) plus the spec prose; re-running the one gate file
  that changed, 4 times total (once whole-file, 3x isolated to `G2 or G4`), is the right-sized
  check for a test-only diff — running DMRG-isolated process gates unrelated to this file would not
  add evidence about this bead's claim and would burn a large fraction of this container's compute
  envelope for no signal.
- **Did not try to reproduce the macOS 27 review's `xpass` finding on this Linux host.** Can't: this
  host's free RHF for square H4 never landed on the C2y-untrusted branch in any of the runs this
  session did (3 full `G2 or G4` runs, all 10 `test_G2` cases outright `PASS` via the trusted
  branch each time) — the conservative/over-inclusion branch of the new assertion is real code
  (it's exercised by whatever host produced the original 20-vs-10 symptom this bead was filed
  over), but unexercised by this session's own runs. Flagged under "what I could not verify."

## What I could not verify

- **The conservative/over-inclusion branch of `test_G2`'s new check** (`c2y_trusted == False`) was
  never exercised on this host this session — every square-H4 bond length tried (1.05, 1.10, 1.19,
  1.35, at both `conv_tol`s) landed on the C2y-trusted branch every time, across 3 repeated runs.
  The branch's assertions (containment, never-forbidden) are exercised by chem-3fd's own
  investigation (cited in spec R2) and were the original motivating symptom, but this session adds
  no new confirmation that the `else` branch's specific assertions are correct beyond reading the
  code and the prior investigation's record.
- **Whether `test_G4[1.20]` would ever observe `broken=True` on a different host/BLAS/LAPACK** —
  not tested beyond this one host, same caveat the first pass noted.
- **Whether `scf_symmetry_status`'s independent fresh RHF solve (used by `_require_broken_rhf`) can
  disagree with the specific `mh` built moments later by `_system(_sq(a))` inside the same test** —
  two independent SCF calls, each with its own basin selection. No disagreement observed in any run
  this session, but this is empirical, not a guarantee.

## Files changed (relative to `main`)

```
$ git diff main...HEAD --stat -- specs/SPEC_eigenstate_reachability.md tests/test_eigenstate_reachability_spec.py
 specs/SPEC_eigenstate_reachability.md      | 15 ++++++-
 tests/test_eigenstate_reachability_spec.py | 62 +++++++++++++++++++++++++---
 2 files changed, 72 insertions(+), 5 deletions(-)
```

Already committed on this branch as `dd185af` (first pass) and `49c9bda` (second pass, on top of
`80be1d7`'s merge of chem-g14/chem-3fd). This session made no further code changes — only verified
and rewrote this handoff file.

## Suggested commands for the human

```bash
git add sandbox-handoffs/chem-0m8.md
git commit   # message: chem-0m8: update handoff to reflect the final branch-aware fix, not the superseded xfail pass
git status
```

`bd close chem-0m8` — closing now, since the committed fix is verified stable on this host (3x
repeated `12 passed, 2 skipped` on `test_G2`/`test_G4`, no hardcoded-branch assertions remain,
`ruff` clean) and the bead is not blocked by any open dependency.
