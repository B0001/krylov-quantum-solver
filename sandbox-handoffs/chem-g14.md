# chem-g14 — `_driver_orbitals` rebuild-reproducibility too strict for conv_tol=1e-9 SCFs

Status: **left open.** This session's job was to close chem-g14 once its two dependencies
(chem-7cx, chem-0m8) landed and `tests/test_eigenstate_reachability_spec.py` is green on Linux.
Both dependencies' fix branches (plus three further branches they themselves spawned:
chem-s6c, chem-0yh, chem-50w) are now merged into `sandbox/chem-g14`, and the file is fully
**deterministic** again — but it is not green: **7/48 still fail**, stably, and root-cause to a
different, still-open bead (chem-8ls) that is outside this one's scope. See "Why this stays open."

Hardware: this container, `Linux f81d0d421506 7.0.14-linuxkit x86_64`, 8 vCPU, 16 GB RAM, pyscf
2.13.1, qiskit-nature 0.8.0. Every run below is single-process (`uv run pytest <file>`), one gate
file at a time, consistent with the block2/pyscf isolation rule — nothing here touches DMRG.

## What I did this session

1. `bd show chem-g14` / `bd show chem-7cx` / `bd show chem-0m8`: both dependencies were marked
   **closed**, but their actual code fixes lived only on their own branches
   (`sandbox/chem-7cx`, `sandbox/chem-0m8`) — neither had been merged into `sandbox/chem-g14` or
   `main`. The bead's own notes said "close once `tests/test_eigenstate_reachability_spec.py` is
   green on Linux," so I merged the dependency branches in to actually get there.
2. Baseline on `sandbox/chem-g14` HEAD (`9a3ac35`, before any merge this session):
   `uv run pytest tests/test_eigenstate_reachability_spec.py -q` → **14 failed, 35 passed**
   (matches the bead's own prior note of "14-15/45 still fail").
3. Merged, one at a time, with `git merge --no-edit`:
   - `sandbox/chem-7cx` (pins `pyscf.lib.num_threads(1)` around `build_molecular_hamiltonian`'s
     free RHF solve — root cause: PySCF's OpenMP integral evaluation is not bit-associative).
   - `sandbox/chem-0m8` (makes `test_G2`/`test_G4` branch-aware instead of asserting a hardcoded
     SCF branch for square H4).
   - `sandbox/chem-s6c` (same thread-pin pattern, applied to `reachability.py`'s own
     `scf_symmetry_status`/`hf_orbital_irreps` RHF calls — a gap chem-7cx's fix didn't cover).
   - `sandbox/chem-0yh` (same pattern, applied to the test file's own `_fci()` helper).
   - `sandbox/chem-50w` (same pattern, applied to `_driver_orbitals`'s own, separately-unpinned
     `run_pyscf()` call).
   All five branch from the same commit (`9a3ac35`), so each merge was close to a fast-forward;
   conflicts were only ever in `.beads/*.jsonl` (append-only interaction logs: unioned; issue
   snapshots: kept whichever side had the more-advanced status, verified with `bd show`) and once
   in `tests/test_eigenstate_reachability_spec.py` and `reachability.py`, both of which
   auto-merged cleanly (no conflict markers). Merge commits:
   `eb4b3ab`, `084c44d`, `6718e4b`, `d4f9ae7`, `71613e4`.
4. Re-ran the gate file after all five merges, three times in a row.

## Evidence

```bash
cd /workspace
uv sync --extra dmrg --extra test --extra dev
uv run ruff check .                                           # -> All checks passed!
uv run pytest tests/test_eigenstate_reachability_spec.py -q   # run 1
uv run pytest tests/test_eigenstate_reachability_spec.py -q   # run 2
uv run pytest tests/test_eigenstate_reachability_spec.py -q   # run 3
```

All three runs: **7 failed, 41 passed, 2 skipped** — identical failing node IDs every time:

```
test_G1_residue_level_identified_warned_and_kept[1.1-1e-06-1e-08]
test_G1_residue_level_identified_warned_and_kept[1.1-1e-06-1e-10]
test_G1_residue_level_identified_warned_and_kept[1.19-1e-06-1e-10]
test_G7_threshold_sites_use_the_population_cut[1.19-1e-06]
test_G9_qksd_from_hf_converges_to_the_residue_level[1.19-1e-06]
test_G9_qksd_from_hf_converges_to_the_residue_level[1.19-None]
test_G9_control_tight_scf_stays_on_the_ag_level[1.1]
```

That exact repeatability (vs. the pre-merge baseline's flicker, and vs. chem-7cx/chem-0m8's own
notes about node IDs "flipping between pass and fail with zero code changes") is itself evidence
that the determinism chain (chem-7cx + chem-s6c + chem-0yh + chem-50w) is complete: every
previously-unpinned RHF/FCI call reachable from this test file is now pinned, so the suite's
result no longer depends on ambient thread count or which run it is.

Related gate files, each run in its own process, to confirm the merged fixes themselves are sound
and to check for collateral regressions:

```bash
uv run pytest tests/test_scf_branch_determinism_spec.py -q            # 4 passed
uv run pytest tests/test_reachability_scf_determinism_spec.py -q      # 5 passed
uv run pytest tests/test_driver_orbitals_scf_determinism_spec.py -q   # 5 passed
uv run pytest tests/test_rebuild_gauge_search_spec.py -q              # 5 passed  (chem-g14's own gate, unaffected)
uv run pytest tests/test_reference_energies.py tests/test_krylov_convergence.py tests/test_noise_resilience.py -q
                                                                        # 20 passed (core gates, unaffected)
uv run pytest tests/test_scf_conv_tol_spec.py -q                      # 4 failed, 10 passed
```

The `test_scf_conv_tol_spec.py` failures (`test_G1[1.2]`, `test_G2[1.35]`, `test_G2[1.19]`,
`test_G3[1.35]`) are **not new** — they match chem-8ls's own enumerated list for that file
(`G1[1.0],G1[1.2], G2[1.35],G2[1.19], G3[1.35]`) and are the same branch-selection phenomenon
discussed below, pre-existing and unrelated to this session's merges.

## Why this stays open

The remaining 7 failures are not a reconciliation-gap or determinism problem — they're a
*which-branch-is-the-witness* problem, which is explicitly chem-8ls's territory (open,
in_progress): "Square-H4 witness specs don't reproduce their recorded Linux numbers; RHF branch
selection is non-deterministic on this container." Direct confirmation:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  uv run python -c "
from reachability import scf_symmetry_status
print('a=1.1 ', scf_symmetry_status('H 0 0 0; H 1.1 0 0; H 1.1 1.1 0; H 0 1.1 0'))
print('a=1.19', scf_symmetry_status('H 0 0 0; H 1.19 0 0; H 1.19 1.19 0; H 0 1.19 0'))
"
# a=1.1  (True, -0.07159850480256758)   <- broken-symmetry branch
# a=1.19 (False, -3.2862601528904634e-14)  <- symmetric branch
```

`a=1.1` landing on `(True, -0.0716...)` matches chem-8ls's own documented finding for that exact
geometry on this exact container ("deterministically lands on the BROKEN-SYMMETRY branch...
instead of the symmetric branch with the 5.07e-10 residue that specs/SPEC_reachability_tolerance.md
...recorded on an earlier Linux run"). `test_G1`/`test_G9_control`/`test_G7`/`test_G9`'s hardcoded
assertions were calibrated against whichever branch an earlier (likely macOS) session landed on,
not this one — the same category of problem chem-0m8 just fixed for `test_G2`/`test_G4` in this
same file, but for a different set of tests that chem-0m8 did not touch.

I did not fix this here: chem-8ls's own text calls this "a scientific decision (which geometries
are stable witnesses here) requiring per-spec re-validation, not a tolerance edit" — exactly the
kind of out-of-scope work this run's instructions say to file, not do. I:

- Added the 7 new/updated node IDs for `test_eigenstate_reachability_spec.py` to chem-8ls's notes
  (`bd update chem-8ls --notes=...`) — its own enumerated list only had `test_G9[1.1]` for this
  file, written before chem-0m8/chem-s6c/chem-0yh/chem-50w's merges changed which tests are
  reachable.
- Wired `chem-g14 → chem-8ls` as a new dependency (`bd dep add chem-g14 chem-8ls`), so this bead
  isn't re-dispatched before chem-8ls actually lands.
- Updated chem-g14's own notes (`bd update chem-g14 --notes=...`) with this session's findings.

## What I decided not to do, and why

- Did not touch chem-8ls itself. It's a pre-existing, separately-tracked, in_progress bead whose
  own text frames the fix as a scientific witness-selection decision across *four* spec files, not
  something to patch inline while closing out chem-g14.
- Did not loosen any of the gate thresholds in `test_G1`/`test_G7`/`test_G9`/`test_G9_control` to
  force a pass. They're not wrong about the physics; they're calibrated to a branch this container
  doesn't reach. Changing the assertion without knowing the "platform-correct" answer is exactly
  what chem-8ls exists to investigate properly.
- Did not re-run the DMRG-isolated gate set (`make gates` / `scripts/run_gates.sh` across all ~117
  spec files) — out of scope for a bead about one test file; ran only the gates that touch the
  merged code paths (listed above) plus the three core validation gates, to bound regression risk
  without spending the whole compute budget on unrelated specs.

## What I could not verify

- Whether chem-8ls's "platform-correct" witness branch for a=1.1/1.19 can be determined without
  a cross-platform SCF audit (e.g., forcing a specific initial guess/DIIS seed and checking it's
  stable on both macOS and this container) — that's the open question chem-8ls itself poses.
- Whether the 5 merged branches' `.beads/issues.jsonl` conflict resolutions (always kept whichever
  side's `status` field was more advanced, confirmed against `bd show <id>` each time) exactly
  match what a `bd dolt pull`-based sync would have produced — I resolved them by hand since all
  five source branches were never pushed/pulled through Dolt in this container.

## Commands to reproduce

```bash
cd /workspace
uv sync --extra dmrg --extra test --extra dev
uv run ruff check .                                                 # All checks passed!
uv run pytest tests/test_eigenstate_reachability_spec.py -q         # 7 failed, 41 passed, 2 skipped (x3 repeats)
uv run pytest tests/test_scf_branch_determinism_spec.py -q          # 4 passed
uv run pytest tests/test_reachability_scf_determinism_spec.py -q    # 5 passed
uv run pytest tests/test_driver_orbitals_scf_determinism_spec.py -q # 5 passed
uv run pytest tests/test_rebuild_gauge_search_spec.py -q            # 5 passed
uv run pytest tests/test_reference_energies.py tests/test_krylov_convergence.py tests/test_noise_resilience.py -q
                                                                      # 20 passed
uv run pytest tests/test_scf_conv_tol_spec.py -q                    # 4 failed, 10 passed (pre-existing, chem-8ls)
```

## Git

On `sandbox/chem-g14`, this session's commits (oldest first):

```
c849ee5  bd: claim chem-g14
eb4b3ab  Merge branch 'sandbox/chem-7cx' into sandbox/chem-g14
084c44d  Merge branch 'sandbox/chem-0m8' into sandbox/chem-g14
6718e4b  Merge branch 'sandbox/chem-s6c' into sandbox/chem-g14
d4f9ae7  Merge branch 'sandbox/chem-0yh' into sandbox/chem-g14
71613e4  Merge branch 'sandbox/chem-50w' into sandbox/chem-g14
20cc7c9  bd: add chem-g14 -> chem-8ls dependency; note chem-8ls's extended Linux failure list
<this handoff + final bd notes commit>
```

Tree otherwise clean except this handoff and the trailing bd-notes export. `.claude/settings.local.json`
(untracked, pre-existing local tool-permission config, unrelated to this task) deliberately not
committed. Per this run's git policy: committed, not pushed, not merged to `main` — a human merges.

## Bead

Left **open**. `bd dep add chem-g14 chem-8ls` run so it depends on chem-8ls going forward. Close
`chem-g14` once `tests/test_eigenstate_reachability_spec.py` is green — which now requires
chem-8ls's witness-geometry work, not more work in chem-g14's own scope.
