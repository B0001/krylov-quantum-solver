# chem-nxn handoff — `hamiltonian_moments` centers by ⟨H⟩ by default

## What this session did

The 2026-10-12 state of this branch already carried the implementation (`moment_expansion.py`'s
`centered: bool = True` + conditioning-guarded fallback, `specs/SPEC_centered_moments_default.md`,
`tests/test_centered_moments_default_spec.py`, and the `centered=False` updates to the two callers
main gained since (`test_centered_pds_roots_spec.py`, `test_centered_pds_shots_spec.py`)) but was
**only verified on macOS**. This session's job, per the bead's own notes, was the three remaining
items:

1. Run the six affected gate files on Linux.
2. Check no other `hamiltonian_moments` caller needs a `centered=False` update.
3. Add the `specs/BACKLOG.md` entry.

No code change was needed — verification confirmed the existing implementation. I did not touch
`moment_expansion.py` or any test file this session.

## Environment

```
$ uname -a
Linux cf6bf70e64f5 7.0.14-linuxkit #1 SMP PREEMPT Sun Sep 20 16:14:09 UTC 2026 x86_64 GNU/Linux
$ nproc
8
$ MemTotal: 16353600 kB  (~16 GB)
```
`uv sync --extra dmrg --extra test --extra dev` completed cleanly (scipy 1.15.3 Linux wheel, no
Accelerate/darwin special-casing needed).

## Gate results (Linux, one pytest process per file, as `run_in_chem.sh`/`make gates` require)

```
$ uv run pytest -q tests/test_centered_moments_default_spec.py
.....                                                    [100%]
5 passed, 5 warnings in 67.15s

$ uv run pytest -q tests/test_moment_pds_spec.py
....                                                     [100%]
4 passed, 10 warnings in 11.86s

$ uv run pytest -q tests/test_pds_excited_roots_spec.py
.....                                                    [100%]
5 passed, 6 warnings in 2.69s

$ uv run pytest -q tests/test_centered_pds_spec.py
......                                                   [100%]
6 passed, 4 warnings in 186.11s

$ uv run pytest -q tests/test_centered_pds_roots_spec.py
....                                                     [100%]
4 passed in 36.55s

$ uv run pytest -q tests/test_centered_pds_shots_spec.py
.......                                                  [100%]
7 passed, 2 warnings in 99.65s
```

Total: **31/31 passed, 0 failed**, matching the macOS pass counts recorded in the bead's
2026-10-12 notes (5, 4, 5, 6, 4, 7) exactly. All warnings are the pre-existing qiskit-nature
`BlueprintCircuit` deprecation notice, unrelated to this change.

`uv run ruff check .` (whole repo, including the five touched test files, `moment_expansion.py`,
and `specs/SPEC_centered_moments_default.md`): **All checks passed.**

I did not run the full `make gates`/`scripts/run_gates.sh` (~90 spec files) — no `.py` file changed
this session (only `specs/BACKLOG.md`), and the caller audit below shows the blast radius is exactly
the six files already run. Re-running the full suite on a cold cache (no prior run in this
container) would mean executing all ~90 gates on a shared 8 vCPU/16 GB box, which is disproportionate
to a change with zero code diff this session. **Unverified:** the rest of the spec suite beyond
these six files, on this container, this session. The 2026-10-12 notes already cover the
implementation's own correctness; nothing in that implementation changed here.

## Caller audit (acceptance criterion: "check no other hamiltonian_moments caller needs
centered=False")

```
$ grep -rln "hamiltonian_moments\|moment_expansion" --include="*.py" . \
    | grep -v "\.git\|\.claude/worktrees\|^tests/\|^moment_expansion.py$"
centered_pds_shots.py

$ grep -n "moment_expansion\|hamiltonian_moments" centered_pds_shots.py
33:from moment_expansion import pds_energy
```

`centered_pds_shots.py` imports only `pds_energy` (which is frame-agnostic — it just solves the
moment linear system and finds roots), not `hamiltonian_moments` itself, so it is not affected by
the new default. Every tracked caller of `hamiltonian_moments` lives in `tests/` and was already
either updated to pass `centered=False` explicitly (the raw-frame studies:
`test_centered_pds_spec.py`, `test_centered_pds_roots_spec.py`, `test_centered_pds_shots_spec.py`,
`test_pds_excited_roots_spec.py`) or left on the new default deliberately
(`test_moment_pds_spec.py`, `test_centered_moments_default_spec.py` — the whole point of G1/G2 is
that the plain default call reproduces/improves on those numbers). `grep` also matched a dozen
`.claude/worktrees/agent-*` directories left over from other sandbox sessions on this host; those
are not part of this repo's working tree and were excluded from the audit and from this bead's
scope.

## BACKLOG entry

Added to `specs/BACKLOG.md` under "Method rungs", immediately after the existing `chem-70c`/
`chem-dcz` entries it builds on:

> CLOSED (chem-nxn) — `hamiltonian_moments` centers by ⟨H⟩ by default, conditioning-guarded.
> SURVIVES. [...] Every existing gate reproduces its original numbers [...] H₂'s K≥3 exhaustion
> boundary correctly falls back to raw: centered cond(Hankel) at K=4 is 6.17e17 vs 5.96e17 raw [...]
> No other `hamiltonian_moments` caller in the repo needed a `centered=False` update.

(Full text in the diff; numbers are quoted from `specs/SPEC_centered_moments_default.md` §6/§G4,
not re-derived this session.)

## Pre-registered criteria vs. outcome

| Criterion (from the bead) | Outcome |
|---|---|
| New spec exists | `specs/SPEC_centered_moments_default.md`, status CLOSED, present before this session |
| `moment_expansion` centers by ⟨H⟩ by default | Confirmed in code (`centered: bool = True`) and by G3 (moments[1] ≈ 0, offset absorbs ⟨H⟩) |
| Existing `moment_pds` gates still pass | `test_moment_pds_spec.py` 4/4 green on Linux, unmodified file |
| K=8 gates pass without callers centering | `test_centered_moments_default_spec.py` G2 — H₄/LiH/N₂ variational+monotone, PDS(8) within `SPEC_centered_pds` thresholds — 5/5 green |

All four criteria hold, now verified on Linux (previously only macOS).

## What I decided not to do, and why

- **Did not re-run `make gates` / the full ~90-file spec suite.** No `.py` changed this session;
  the caller audit above shows the change's blast radius is exactly the six files run. Running the
  full suite cold would burn significant time/compute on gates this bead's diff cannot affect.
- **Did not re-implement or modify `moment_expansion.py`, the spec, or any test file.** The prior
  worker's implementation and spec were already correct and complete; this session's job was
  verification, not re-authoring. (The spec's own attribution note already explains the shared-tree
  provenance of that code; I did not re-litigate it.)
- **Did not touch the two `.claude/worktrees/agent-*/...` copies of `moment_expansion.py` that
  `grep -r` surfaced** (stale worktrees from unrelated sandbox sessions on this host, several
  `chem-7ko`/`chem-8tl`/etc. branches) — out of scope for this repo's working tree and this bead.

## Commands to reproduce

```bash
cd /workspace
uv sync --extra dmrg --extra test --extra dev
uv run pytest -q tests/test_centered_moments_default_spec.py
uv run pytest -q tests/test_moment_pds_spec.py
uv run pytest -q tests/test_pds_excited_roots_spec.py
uv run pytest -q tests/test_centered_pds_spec.py
uv run pytest -q tests/test_centered_pds_roots_spec.py
uv run pytest -q tests/test_centered_pds_shots_spec.py
uv run ruff check .
```

## Git state at handoff

Branch `sandbox/chem-nxn`, tip `a0f17d8` plus this session's commit(s): `specs/BACKLOG.md` entry
and `.beads/issues.jsonl` (export side effect of `bd update --claim`/`bd close`). Tree is otherwise
clean. I did not commit, push, merge, or rebase beyond what this task's git-policy override
authorizes (commit on this branch only); the host pushes.

Per the task instructions, I am committing this handoff and the BACKLOG change on this branch
myself, and closing the bead once committed.
