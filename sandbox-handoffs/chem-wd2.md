# chem-wd2 handoff

## What this was

`specs/BACKLOG.md` had two `[ ]` (open) entries whose questions were already answered by later
`[x]` (done) entries in the same file, per the bead's own reproduction. Pure documentation hygiene:
no code or spec content needed to change, just flip the two stale checkboxes and point them at the
entries/specs that answer them, per the file's own status-key convention
(`[ ] open · [~] specced · [x] done (link the spec) · [-] killed`, `specs/BACKLOG.md:8`).

## What I changed

`specs/BACKLOG.md`, two edits, both pure prepended cross-references — no existing text removed:

1. Line 96 (`Does SPEC_subspace_floor_resolvability's mechanism survive?`): `[ ]` → `[x]`, with a note
   pointing at "The floor-guard mechanism survives — the ~1e-4 level is PHYSICAL, not the SCF
   artifact" (now at line ~201), which is explicitly filed as the follow-up to this exact question
   (same `conv_tol=1e-13` re-run, same 8 witnesses) and confirms the mechanism as originally written.

2. Line 224 (`The krylov_refine stub is not a marginal tightening — it may moot the block
   certificate`): `[ ]` → `[x]`, with a note pointing at "CLOSED 2026-08-01 — the stub is implemented
   and it works in SELF mode" (line ~163) and at `SPEC_chained_overlap.md`, whose header literally
   says "Closes the backlog hypothesis". Re-verified in code during this session (see below) that the
   stub is gone and the function is implemented.

Verified before editing (per the bead's own evidence, re-checked independently this session):
```
$ find . -maxdepth 1 -name krylov_refine.py
(nothing)
$ grep -rn "def refine_via_lanczos" --include=*.py .
hybrid_quantum_solver/certified_overlap/krylov_refine.py:51:def refine_via_lanczos(
```
`hybrid_quantum_solver/certified_overlap/krylov_refine.py` is a real docstring-first implementation,
not a `NotImplementedError` stub — confirms entry 2's premise is stale.

```
$ head -4 specs/SPEC_chained_overlap.md
# SPEC: the chained overlap bound — the stub, implemented, and it works in self mode

**Status:** IMPLEMENTED once gates green. Closes the backlog hypothesis *"The `krylov_refine` stub is
not a marginal tightening"*.
```
Confirms the spec's own header claims closure of exactly this backlog entry.

Full diff (also `git diff -- specs/BACKLOG.md` on this branch):
```diff
-- [ ] **Does `SPEC_subspace_floor_resolvability`'s mechanism survive?** *(follow-up to the
++ [x] **Does `SPEC_subspace_floor_resolvability`'s mechanism survive?** *(ANSWERED below — see the
++  "The floor-guard mechanism survives" entry, filed as the direct follow-up to this exact question with
++  the same conv_tol=1e-13 re-run and the same witnesses: the mechanism is confirmed as originally
++  written, PHYSICAL not a spurious SCF artifact.)* *(follow-up to the
...
-- [ ] **The `krylov_refine` stub is not a marginal tightening — it may moot the block certificate**
++ [x] **The `krylov_refine` stub is not a marginal tightening — it may moot the block certificate**
++  *(SUPERSEDED — see the "CLOSED 2026-08-01 — the stub is implemented and it works in SELF mode" entry
++  above and [`SPEC_chained_overlap.md`](SPEC_chained_overlap.md), whose header states it "closes the
++  backlog hypothesis" this entry states; `krylov_refine.py` no longer exists at the top level and
++  `refine_via_lanczos` is implemented at
++  `hybrid_quantum_solver/certified_overlap/krylov_refine.py:51`, not a stub.)*
```

## Gates / tests

No code changed — only a markdown file (`specs/BACKLOG.md`) was edited, and the change is additive
prose (checkbox flip + a cross-reference note), so no `tests/test_*_spec.py` gate exercises this
content. There is no automated linter over `BACKLOG.md`'s prose/checkbox format in this repo (checked:
`grep -rln BACKLOG tests/ specs/README.md` only turns up docstring mentions in unrelated spec tests,
not a format validator). I did not run `make gates` / `make test` since nothing they cover changed;
running the full DMRG-isolated suite for a two-line doc edit would not exercise the change and would
burn the compute envelope for no signal. `ruff` was not available in this container's environment
(`uv run ruff check .` → `Failed to spawn: ruff`, no venv had it installed) and is irrelevant to a
markdown-only diff regardless.

## Pre-registered criteria (from the bead)

"The two stale `[ ]` lines in `specs/BACKLOG.md` are updated to reflect that they are already
answered (cross-referenced to the entries/specs that answer them), consistent with the file's own
status-key convention. No other content changes."

- Both lines updated to `[x]` with cross-references: **PASS**.
- No other content changed: **PASS** — `git diff` above touches only the two checkbox markers plus
  one inserted parenthetical note per entry; nothing else in the 400+ line file moved.

## What I decided not to do, and why

- Did not touch the *other* struck-through-but-still-`[ ]` entries in the file (e.g. line 106 "Unify
  the reachability tolerance", line 211 "UPDATE 2026-07-30... BLOCKED") even though they use similar
  "superseded"/"original text follows" language. The bead named exactly two entries by content and
  approximate line number; those two are still genuinely open/blocked per their own text (not fully
  answered elsewhere), so re-flagging them would be scope creep beyond what was reproduced and
  asked for. If they're also considered stale, that's a separate bead.
- Did not renumber or re-verify every other `[ ]`/`[x]` pairing in the 400-line file — out of scope;
  the bead's reproduction was specific to these two.

## What I could not verify

Nothing outstanding — both cross-references were independently re-derived from the current repo
state (file existence checks, grep, and reading the referenced entries/spec header in full) rather
than taken on the bead's word.

## Git

Working tree is dirty with the one-file diff above, on branch `sandbox/chem-wd2`, ready to commit.
Suggested commands for the human reviewer:
```
git add specs/BACKLOG.md sandbox-handoffs/chem-wd2.md
git commit -m "docs(backlog): flip two stale [ ] entries to [x], cross-reference their answers"
```
I did not commit per the git policy for this session type — actually note: this session's *task*
instructions say to commit on `sandbox/chem-wd2` directly (superseding the general "do not commit"
policy). Committing now as instructed.

Bead `chem-wd2` closed with acceptance criteria met.
