# chem-0wh — vsdp_hchain.py noise-schedule misalignment for single-value `--bond-dims`

## Verdict: confirmed and fixed.

The review note (PR #55, batch unit 8) was right. `dmrg_upper` in `vsdp_hchain.py` built the
block2 noise schedule as:

```python
noises=(1e-4,) * (len(bond_dims) - 2) + (1e-5, 0.0)
```

For `len(bond_dims) >= 2` this is fine: `(len(bond_dims) - 2)` is `>= 0`, and the result has
exactly `len(bond_dims)` entries, matching `bond_dims` and `thrds=[tol]*len(bond_dims)`.

For `len(bond_dims) == 1` (e.g. `--bond-dims 500`), `len(bond_dims) - 2 == -1`, and
`(1e-4,) * -1` is `()` in Python — but the `+ (1e-5, 0.0)` tail is unconditional, so `noises`
ends up `(1e-5, 0.0)`, **length 2**, one longer than `bond_dims`/`thrds` (length 1). block2's
`DMRGDriver.dmrg` doesn't require these lists to be the same length — it indexes them
independently per sweep, clamping each to its own last entry once the sweep index runs past the
list — so the extra phantom noise entry doesn't raise; it silently shifts which sweep each noise
value lands on. A single-bond-dim ladder (one real "stage") ended up with sweep 0 noise 1e-5
(left over from a nonexistent penultimate stage) and only sweeps 1+ at 0.0, instead of the
intended "last stage gets 0 noise" semantics that every other ladder length gets.

## Reproduction (before fixing, per repo rule)

Monkeypatched `hybrid_quantum_solver.dmrg_reference.dmrg_energy` (and the subprocess/file I/O
around it) to capture the exact `noises`/`bond_dims` kwargs `dmrg_upper` builds, without running
real DMRG:

```
bond_dims=(500,)        (len 1) -> noises=(1e-05, 0.0) (len 2)   MISALIGNED
bond_dims=(250, 500)    (len 2) -> noises=(1e-05, 0.0) (len 2)   ok
bond_dims=(250,500,1000)(len 3) -> noises=(0.0001, 1e-05, 0.0) (len 3)  ok
```

Also confirmed block2 itself tolerates the length mismatch rather than erroring (ruling out "it
would just crash" as a lesser bug): real `dmrg_energy` call on a toy 2-orbital system with
`bond_dims=(50,)`, `noises=(1e-5, 0.0)` (mismatched) ran to completion. On that trivial active
space the extra noise sweep didn't change the converged energy (Hilbert space too small to have a
local minimum to get stuck in), so this bug is a schedule-construction defect, not something
demonstrated here to move a physical answer — the risk is in larger active spaces where the
intended "zero noise on the only/last stage" vs. actual "1e-5 noise on sweep 0" distinction can
matter for convergence to a local minimum. No recorded production number in `data/` depends on a
single-value `--bond-dims` run (default is `--bond-dims 250,500,1000,1500`, length 4); this is a
latent bug for a CLI path that, as far as the repo's tracked files show, hasn't been exercised at
`len==1` for a real result yet.

## Fix

`vsdp_hchain.py` (new `_noise_schedule` helper, used in `dmrg_upper`):

```python
def _noise_schedule(n_stages):
    tail = (1e-5, 0.0)
    if n_stages >= len(tail):
        return (1e-4,) * (n_stages - len(tail)) + tail
    return tail[len(tail) - n_stages:]
```

`_noise_schedule(n)` always returns exactly `n` entries:
- `n=1` -> `(0.0,)` (the lone stage gets the terminal zero-noise sweep, consistent with "last
  stage(s) are noise-free")
- `n=2` -> `(1e-5, 0.0)` (byte-identical to the old formula)
- `n=3` -> `(0.0001, 1e-05, 0.0)` (byte-identical to the old formula)
- `n>=2` in general reproduces the old formula exactly; only `n<2` changes.

Re-ran the same capture harness after the fix: `bond_dims=(500,)` now produces `noises=(0.0,)`
(len 1, matches), and the `len==2`/`len==3` cases are byte-identical to before.

## Gates run

```
$ uv run ruff check .
All checks passed!

$ uv run pytest tests/test_vsdp_hchain_noise_schedule.py -v      # new regression test
4 passed in 0.43s

$ uv run pytest tests/test_hchain_vsdp_bracket_spec.py -q        # the spec gate that imports/exercises vsdp_hchain
7 passed in 85.90s (0:01:25)

$ uv run python -m pytest tests/ --ignore-glob='tests/test_*_spec.py' --ignore=tests/test_dmrg_reference.py -q
60 passed, 4 skipped, 25 warnings in 47.84s
```

I did **not** run the full `tests/test_dmrg_reference.py tests/test_*_spec.py` process (114 spec
files, real DMRG/FCI/benchmark runs, historically hours) — this change touches exactly one helper
used by exactly one CLI path (`vsdp_hchain.py --dmrg-upper`) that no other spec file calls, and
the one spec file that does import `vsdp_hchain` (`test_hchain_vsdp_bracket_spec.py`) is green.
Running the rest would burn the shared 8-CPU/16GB envelope on ~113 files this diff cannot affect.
If a reviewer wants the full gate run for certainty, `scripts/run_gates.sh` (or `make gates`) is
the command; I'd estimate it at well over the 10-minute budget I have per command here, so flagging
rather than running it.

## Hardware / wall time

`uname -a`: `Linux 72e682cc0d98 7.0.14-linuxkit x86_64`, `nproc`: 8. Longest single command was
the vsdp bracket spec gate at 85.9s; nothing else here took more than a minute.

## Pre-registered criteria (repo rule: reproduce before fixing)

| Criterion | Result |
|---|---|
| Reproduce the misalignment before touching code | PASS — captured `noises` length 2 vs. `bond_dims` length 1 at `len(bond_dims)==1`, via a dmrg_energy-call-capturing harness, no fix applied yet |
| Fix generalizes correctly for every `len(bond_dims)` | PASS — `_noise_schedule` tested for `n_stages` 1..5; `n>=2` byte-identical to the old formula (regression-safe for the default `--bond-dims 250,500,1000,1500` path and the existing n=2/n=3 cases), `n<2` now matches length |
| No regression in the spec gate that actually imports `vsdp_hchain` | PASS — `test_hchain_vsdp_bracket_spec.py`, 7/7 |
| Lint clean | PASS — `ruff check .` |

## What I decided not to do, and why

- Did not run the full DMRG-isolated gate process (114 spec files). Scope: this diff changes one
  helper feeding one CLI path (`--dmrg-upper`) that none of those 114 files call except the one I
  did run. Given the repo's own compute-envelope guidance (8 CPU / 16 GB shared, time a mid-size
  point before committing to a large one), running ~113 unrelated spec files — several of which
  are documented multi-hour DMRG ladders — for a change with zero import-time or call-time
  overlap would be guessing at thoroughness rather than buying it. Flagging this explicitly as
  **unverified**: I have not re-run the full `make gates` / `make test` second stage after this
  change.
- Did not attempt to demonstrate the bug moving a real physical answer (e.g., a converged DMRG
  energy on a nontrivial active space). The toy 2-orbital test shows block2 accepts the mismatched
  lengths without erroring; demonstrating an actual wrong answer would need an active space large
  enough to have a real local-minimum risk, which is exactly the "time a mid-size point before
  committing to a large one" territory I didn't have grounds to spend on for a bug that's already
  unambiguous at the list-construction level.
- Did not touch `n_sweeps=4 * len(bond_dims)` or `thrds=[tol]*len(bond_dims)` — both already scale
  correctly with `len(bond_dims)` and were not part of the reported misalignment.

## Files changed

- `vsdp_hchain.py` — added `_noise_schedule`, used it in `dmrg_upper` in place of the inline tuple
  expression.
- `tests/test_vsdp_hchain_noise_schedule.py` — new regression test (4 cases) pinning
  `_noise_schedule`'s length-matching and byte-identical-for-`n>=2` behavior.

## Suggested commands for the human

```bash
git add vsdp_hchain.py tests/test_vsdp_hchain_noise_schedule.py .beads/issues.jsonl sandbox-handoffs/chem-0wh.md
git commit -m "chem-0wh: fix vsdp_hchain noise-schedule misalignment for single-value --bond-dims

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

(`.claude/settings.local.json` is untracked and predates this session — left alone, not part of
this change.)

Bead closed with `bd close chem-0wh` after this handoff was written.
