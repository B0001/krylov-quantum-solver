# chem-dy4 handoff

## Outcome: closed with `--reason` (premise stale / superseded) — no code changed

chem-dy4 predicted that `tests/test_reachability_tolerance_spec.py::test_G2_...` and
`::test_G5_constant_governs_the_certified_reference` would fail on Linux because PR #62 made
`reachable_eigenpairs` **veto** symmetry-forbidden levels, which would change what
`exact_reachable_overlap` returns at the square-H4 a=1.1 witness. That mechanism no longer exists:
PR #64 (already merged into `main` at `b0d82fd`, visible in this tree) reverted the veto's effect
on the mask. Current `reachability.py` docstring, verbatim:

> "Vetoing that level is nevertheless WRONG for a reference ... A reference must describe what the
> data targets, so `reachable_mask` IS the population cut."

So G2/G5 cannot be failing for the reason the bead names, because the thing it names (the veto
changing the mask) was already undone before I started. I verified this by reading
`reachability.py` and `hf_overlap_certificate.exact_reachable_overlap` (still a plain population
cut through `reachable_eigenpairs`, no symmetry veto applied to the mask).

## What I actually found (reproduce before you file)

Ran the live gate on this container to check anyway:

```
uv run pytest tests/test_reachability_tolerance_spec.py -v
```

Result: **6 failed, 6 passed** (set below), but re-running the identical command with no code
changes repeatedly gives a **different failure count each time** (6/6, 7/5, 6/6, 6/6 across 5
runs, pinned-threads and not). Representative run:

```
FAILED test_G1_two_thresholds_select_different_ground_states
FAILED test_G2_the_two_specs_d1_references_disagree
FAILED test_G3_a_level_sits_between_the_two_thresholds
FAILED test_G4_thresholds_agree_away_from_the_witness[1.2]
FAILED test_G5_constant_governs_the_certified_reference
FAILED test_G6_the_disputed_amplitude_is_an_SCF_convergence_residue
6 failed, 6 passed, 14 warnings in 9.49s
```

Root cause (confirmed, not the veto): the spec's witness, square H4 a=1.1 Å, is supposed to sit on
the RHF **symmetric** branch with a tiny residue (`p0 = 5.07e-10` at default `conv_tol`, recorded in
`specs/SPEC_reachability_tolerance.md` §2b/§10 from an earlier Linux session). On **this** container
it deterministically lands on the **broken-symmetry** branch instead:

```
$ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    uv run python -c "
from reachability import scf_symmetry_status
print(scf_symmetry_status('H 0 0 0; H 1.1 0 0; H 1.1 1.1 0; H 0 1.1 0'))"
(True, -0.0716...)        # spec's own session recorded (False, ~-5e-15) here
```

`amp2[0]` (lowest level's HF population) comes out to `0.4692280862...` every time I checked it
with threads pinned to 1 — a real, O(1) physical overlap, not a residue — so this is not the
tiny-residue BLAS-thread-nondeterminism chem-ffe already logged (that one swings ~1e-51 to ~1e-8).
This is a full branch flip, deterministic on this machine, different from the number the spec
recorded on whatever Linux machine wrote it.

Separately, *without* thread pinning the exact pass/fail subset of this same file varies run to
run (confirmed 3+ times) — that part **is** the chem-ffe phenomenon (unpinned OpenBLAS thread count
nondeterminism near a degenerate SCF), now also confirmed in this file (chem-ffe's own notes only
name `test_scf_conv_tol_spec.py` / `test_eigenstate_reachability_spec.py`).

I also spot-checked the three sibling spec files chem-ffe's notes flagged as having the same kind
of macOS/Linux divergence, to see how far it spreads:

```
uv run pytest tests/test_eigenstate_reachability_spec.py -q   # 17 failed, 28 passed
uv run pytest tests/test_scf_conv_tol_spec.py -q               # 5 failed, 9 passed
uv run pytest tests/test_symmetry_reachability_spec.py -q      # 5 failed, 9 passed
```

All four files are hitting variants of the same underlying issue: square-H4's RHF
symmetric/broken-symmetry branch selection at these particular geometries is not stable across
machines (OpenBLAS DYNAMIC_ARCH/Haswell here, numpy 2.4.6, scipy 1.15.3, pyscf 2.13.1) and, in some
cases, not stable run-to-run on the same machine either.

## Hardware / versions (this container)

```
Linux 9e3695e5bd59 7.0.14-linuxkit x86_64, nproc=8
numpy 2.4.6, scipy 1.15.3 (OpenBLAS 0.3.31.188.0 DYNAMIC_ARCH Haswell, MAX_THREADS=64)
pyscf 2.13.1, qiskit 2.5.0, qiskit_nature 0.8.0
```

## Why I did not try to fix it

The bead's prescribed fix ("update them to record the veto as the fix") targets a mechanism that
no longer exists. The failure that *does* exist is a pre-existing, already-tracked-elsewhere
problem (chem-ffe) that's much bigger than this bead's scope: picking a Linux-stable witness
geometry (or a branch-selection-robust reformulation) for four spec files is a scientific decision
requiring per-spec re-validation against new reference numbers, not a tolerance or assertion edit.
Doing that now would be scope creep beyond what chem-dy4 asked for and would duplicate/preempt
chem-ffe's remit.

## Filed

`chem-8ls` — "Square-H4 witness specs don't reproduce their recorded Linux numbers; RHF branch
selection is non-deterministic on this container" (P2, bug), with the repro commands above, the
full list of affected gates per file, and the acceptance criterion (deterministic pass across >=5
repeated runs, or re-pinned witnesses with new recorded numbers). Left open for whoever takes on
the Linux-witness-stability work; it supersedes/extends chem-ffe's existing note for this specific
file.

## Beads

- `chem-dy4`: **closed**, `--reason`: premise superseded — the predicted veto-vs-mask mechanism was
  already reverted by PR #64 before this session started; the Linux failures that do exist have a
  different, larger, already-partially-tracked cause, filed as `chem-8ls`.
- `chem-8ls`: created, open, P2.
- `chem-ffe`: left untouched (not claimed, not edited) — `chem-8ls` references it and extends its
  finding to a file chem-ffe's own description doesn't name.

## Gate / test runs (verbatim tails)

```
$ uv run pytest tests/test_reachability_tolerance_spec.py -v  (one of several runs, counts vary)
6 failed, 6 passed, 14 warnings in 9.49s

$ uv run pytest tests/test_eigenstate_reachability_spec.py -q
17 failed, 28 passed, 35 warnings in 57.85s

$ uv run pytest tests/test_scf_conv_tol_spec.py -q
5 failed, 9 passed, 22 warnings in 18.38s

$ uv run pytest tests/test_symmetry_reachability_spec.py -q
5 failed, 9 passed, 20 warnings in 11.43s
```

## What I did not do / could not verify

- Did not touch `specs/SPEC_reachability_tolerance.md`, `reachability.py`,
  `hf_overlap_certificate.py`, or any test file — no code changes this session.
- Did not attempt to find a Linux-stable witness geometry for any of the four files — that's
  `chem-8ls`'s job, not this bead's.
- Did not run `make gates` / `make test` end to end (would re-trigger the same block2/pyscf process
  isolation concerns across the whole suite for no new information beyond what's above; the
  specific gates named by chem-dy4 and its neighbors are covered above).
- Did not determine *why* this container's OpenBLAS build picks the opposite RHF branch from the
  session that recorded the spec's numbers (plausible candidates: OpenBLAS DYNAMIC_ARCH kernel
  dispatch, numpy/scipy minor-version drift, or just float summation order sensitivity at a
  genuinely near-degenerate SCF) — flagged as an open question in `chem-8ls`, not run down further.

## Git status

No files changed. Tree is clean except pre-existing untracked `.claude/settings.local.json` (not
mine, not touched).

```
$ git status
On branch sandbox/chem-dy4
Untracked files:
  (use "git add <file>..." to include in what will be committed)
	.claude/settings.local.json
nothing added to commit but untracked files present (use "git add" to track)
```

Nothing to commit for this bead. The handoff file itself (`sandbox-handoffs/chem-dy4.md`) and the
beads DB update are the only artifacts; the beads update is captured in Dolt per the usual beads
workflow (`bd close` below), not via `git commit` (no source files changed).

Proposed command for the human (nothing requires it, but if they want this handoff file tracked):
```
git add sandbox-handoffs/chem-dy4.md
git commit -m "chem-dy4: handoff -- premise superseded by PR #64, filed chem-8ls"
```
