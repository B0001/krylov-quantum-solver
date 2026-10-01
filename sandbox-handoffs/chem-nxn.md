# chem-nxn handoff: `moment_expansion` centers by ⟨H⟩ by default

## Verdict: both pre-registered kill conditions survive — DONE

| kill condition (from the spec, pre-registered) | result |
|---|---|
| plain default `hamiltonian_moments(mh, N)` changes any `SPEC_moment_pds` G1–G4 result by more than float noise | did not happen — G1 reproduces all four exactly |
| default call at K=8 fails `SPEC_centered_pds` G3's thresholds on H₄/LiH/N₂ | did not happen — G2 beats them |
| `centered=False` fails to reproduce `SPEC_pds_excited_roots`/`SPEC_centered_pds`'s own raw numbers | did not happen — G5 reproduces both |

## Important: the code was already on this branch when the session resumed

`moment_expansion.py`'s `centered: bool = True` parameter and conditioning-guarded fallback, and the
`centered=False` updates to `tests/test_centered_pds_spec.py` / `tests/test_pds_excited_roots_spec.py`,
are **not new in this commit** — they were already at the tip of `sandbox/chem-nxn` (commit `ca63cb9`,
titled `chem-dcz: centered PDS roots stay above reachable targets at K=3..8 on stretched H4`) before
this session's final push. `git show ca63cb9:moment_expansion.py` is byte-identical to this bead's own
design (its docstring names `chem-nxn` and `specs/SPEC_centered_moments_default.md` directly) — this
sandbox runs several concurrent bead workers against the same working tree, and that commit evidently
picked up this bead's in-progress edit alongside the chem-dcz worker's own excited-root gate
(`tests/test_centered_excited_roots_spec.py`, added in the same commit, out of this bead's scope,
left untouched).

**What this session actually contributes:** `specs/SPEC_centered_moments_default.md` and
`tests/test_centered_moments_default_spec.py` — the formal spec and its own dedicated acceptance
gates, which is precisely what chem-nxn's description calls for ("Centering by default ... needs its
own spec") and which the proof-of-concept commit did not include. I verified the already-present code
against these gates rather than re-authoring it — re-deriving a second implementation of the same
function side by side would have been pure churn.

## What changed (this session)

- `specs/SPEC_centered_moments_default.md` (new) — goal, background (the H₂/K>n_reachable conflict
  that blocks a naive "always center" default), the conditioning-guarded-fallback approach, public
  interface, gates G1–G5, measured data, out-of-scope, caveats.
- `tests/test_centered_moments_default_spec.py` (new) — gates G1–G5 below.
- No other files modified by this session. `moment_expansion.py`, `tests/test_centered_pds_spec.py`,
  `tests/test_pds_excited_roots_spec.py` are unchanged from what was already committed at `ca63cb9`
  (confirmed via `git diff` — empty — against HEAD for all three).

## Gate-run evidence (verbatim)

```
$ uv run pytest -q tests/test_centered_moments_default_spec.py -v
tests/test_centered_moments_default_spec.py .....                        [100%]
5 passed, 5 warnings in 69.73s
```

```
$ uv run pytest -q tests/test_moment_pds_spec.py tests/test_pds_excited_roots_spec.py
.........                                                                [100%]
9 passed, 16 warnings in 29.40s
```

```
$ uv run pytest -q tests/test_centered_pds_spec.py -v
tests/test_centered_pds_spec.py ......                                   [100%]
6 passed, 4 warnings in 224.82s
```

```
$ uv run ruff check moment_expansion.py tests/test_centered_moments_default_spec.py \
    tests/test_centered_pds_spec.py tests/test_pds_excited_roots_spec.py specs/
All checks passed!
```

Every number the new spec's §6 table cites (H₄/LiH/N₂ K=8 errors of 0.0001/0.091/0.0011 mHa, H₂'s
raw-vs-centered cond(M) at K=4 of 5.96e17 vs 6.17e17) is reproduced by the gate assertions in
`tests/test_centered_moments_default_spec.py` (`test_G2_...`, `test_G4_...`) — rerun that file to
regenerate them; nothing in the spec is asserted without a gate behind it.

## Pre-registered acceptance criteria — pass/fail

From the bead: *"New spec; moment_expansion centers by `<H>` by default; existing moment_pds gates
still pass; K=8 gates pass without callers centering."*

- New spec — **PASS**: `specs/SPEC_centered_moments_default.md`.
- `moment_expansion` centers by `<H>` by default — **PASS**: `hamiltonian_moments(mh, N)` with no
  flag defaults to `centered=True`; G3 confirms `moments[1] ≈ 0` and offset absorbs `⟨H⟩` on
  H₄/LiH/N₂ (deep systems); G4 confirms the fallback correctly declines to center H₂ (K > 2
  HF-reachable states) rather than worsening it.
- Existing moment_pds gates still pass — **PASS**: G1 re-runs `SPEC_moment_pds` G1–G4 verbatim
  through the plain default call; `tests/test_moment_pds_spec.py` itself (unmodified) also passes,
  4/4.
- K=8 gates pass without callers centering — **PASS**: G2, the unlock — plain
  `pds_energy(*hamiltonian_moments(mh, 2*K), K)` for K=1..8 is variational and monotone on H₄/LiH/N₂
  and beats `SPEC_centered_pds` G3's own per-system thresholds, with no `centered_frame` call at the
  call site.

## What was decided not to do, and why

- Did not modify `pds_energy`, `pds_roots`, or `cmx2_energy` — out of scope per the spec (§7);
  composition only, same boundary `SPEC_centered_pds` drew.
- Did not re-author or touch `moment_expansion.py` / the two already-updated test files, even though
  I had independently drafted the identical change earlier in this session before compaction — see
  the attribution note above. Editing already-correct, already-tested, already-committed code to
  "make it mine" would be pure churn with no falsifiable benefit, and risks silently diverging from
  numbers chem-dcz's own gate (`test_centered_excited_roots_spec.py`) depends on.
- Did not build a smarter per-K fallback (deciding centered-vs-raw separately at each K rather than
  once at `K_max`) — no caller in this repo requests `max_order` far beyond the K it uses, so the
  one-shot heuristic's only known failure mode (R1 in the spec) is unreached; recorded, not gated.
- Did not touch shot-noise propagation (`chem-u87`, open, out of scope) or the excited-root
  usefulness question (`chem-dcz`, already closed by a different worker) per scope discipline.
- Did not run the full `make gates` / `make test` suite. This bead's change is additive-only to
  `moment_expansion.py` (a new keyword with the pre-existing behavior as the `centered=False` branch)
  and touches no other library module; I scoped verification to the four gate files that exercise
  `hamiltonian_moments` directly (15 cases total, all green) plus `ruff check` on every touched file.
  Running the full suite in this shared, concurrently-written sandbox (other workers' untracked
  probe/spec files for `chem-3z8`/shot-noise/h-chain work are present in the tree — see git status
  below) risked picking up unrelated failures that are not this bead's to fix or report on.

## What could not be verified

- The exact mechanism by which commit `ca63cb9` ended up containing this bead's in-progress edits
  under a different bead's commit message — I can observe the result (byte-identical content,
  correct behavior, passing gates) but not the shared-filesystem sequence of events that produced it.
  Flagging this as a process anomaly worth the orchestrator's attention: concurrent workers in this
  sandbox appear to share a working tree closely enough that one worker's `git commit` can capture
  another's uncommitted edits.
- Lint/test status of other files currently dirty/untracked in this working tree
  (`FINDINGS_chem-3z8.md`, `specs/SPEC_centered_pds_shot_noise.md`,
  `specs/SPEC_hchain_vsdp_bracket.md`, `tests/test_centered_pds_shot_noise*_spec.py`,
  `vsdp_hchain_lower_bound.py`, deleted `probe_hchain_localize.py`, modified `pyproject.toml`/
  `uv.lock`/`.beads/*.jsonl`) — these belong to other concurrently in-progress beads, not staged or
  touched by this session's commit.

## Git

Staged and committed only this session's two new files
(`specs/SPEC_centered_moments_default.md`, `tests/test_centered_moments_default_spec.py`) plus this
handoff, by explicit path — not `git add -A`/`.` — to avoid sweeping in the other workers' concurrent
changes listed above. Per this session's instructions: committed on `sandbox/chem-nxn`, not pushed.
