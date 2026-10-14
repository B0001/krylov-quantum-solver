# chem-g14 — `_driver_orbitals` rebuild-reproducibility too strict for conv_tol=1e-9 SCFs

Status: **left open, fix implemented and verified, but the bead's own regression gate
(`tests/test_eigenstate_reachability_spec.py` fully green) does not hold** — see "Why this stays
open" below.

## What changed

- `reachability.py` — `_driver_orbitals(mh)`:
  - `_SAME_OPERATOR_ATOL` loosened from `1e-12` to `1e-4` (the old bound rejected every rebuild on
    this host/pyscf version; see the bead's own repro).
  - New energy pre-filter: compares the rebuild's `drv._calc.e_tot` against `mh.hf_energy` to
    `_ENERGY_MATCH_ATOL = 1e-6` Ha before trusting the operator comparison; on mismatch, retries the
    SCF (new `PySCFDriver` run) up to `_MAX_REBUILD_RETRIES = 5` times rather than accepting or
    gauge-searching a genuinely different local minimum.
  - New bounded discrete sign-gauge search: if the energy matches but the unflipped operator
    comparison still fails, searches `2**(k-1)` per-active-orbital sign-flip patterns (one column's
    sign fixed as gauge redundancy) for `k - 1 <= _MAX_SIGN_SEARCH_ORBITALS = 10` active orbitals,
    re-deriving the qubit operator from the sign-flipped `mo_coeff` (cheap: no new SCF) and
    re-comparing.
  - No change to the function's signature or return type.
- `specs/SPEC_rebuild_gauge_search.md` — new spec (goal, root-cause analysis, approach, gates,
  honestly revised status/G2, out-of-scope).
- `tests/test_rebuild_gauge_search_spec.py` — new gate file, G3-G5 (cheap when unneeded, energy
  pre-filter behavior, sign-search bound). 5 tests, all passing.
- Two new beads filed for issues found but out of this bead's scope (see below): chem-7cx, chem-3fd.

## Pre-registered acceptance criteria and pass/fail

The bead (`bd show chem-g14`) frames the fix as closing the rebuild-reproducibility gap so the 17
failures in `tests/test_eigenstate_reachability_spec.py` clear. Breaking that into checkable parts:

1. **The rebuild-reproducibility gap itself is closed** (direct `_driver_orbitals` calls succeed
   where they used to always return `None`). ✅ **Pass.** Calling `reachability._driver_orbitals(mh)`
   directly (bypassing the rest of the test file) for every `(bond length, conv_tol)` combination
   the gate file exercises, 5-6 trials each: **35/36 successes**, vs. an unconditional `None` 5/5
   times at the old `_SAME_OPERATOR_ATOL=1e-12` on the same witness (a=1.2).
2. **No regression — every baseline-passing test still passes.** ✅ **Pass.** Diffed baseline
   `main` (unmodified `reachability.py`, via `git stash`) against the fix across 3 repeated full
   runs of `tests/test_eigenstate_reachability_spec.py`:
   - Baseline: **17 failed, 28 passed** (stable across repeats, though the exact failing IDs
     flicker — this instability is itself the bead's own prior finding).
   - Fixed: **14, 15, 14 failed / 34, 35, 35 passed** across 3 runs.
   - Every test passing on baseline still passes with the fix. Tests that newly pass:
     `test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged[1.4]`,
     `test_G5_stale_build_args_disable_spatial_bits_but_a_shift_does_not`,
     `test_G9_qksd_from_hf_converges_to_the_residue_level[1.19-None]`.
3. **The test suite fully passes (45/45).** ❌ **Fail.** 14-15/45 still fail. This is the reason
   the bead stays open — see below.
4. **The fix's own mechanics are correctly bounded** (G3-G5 in the spec). ✅ **Pass.**
   `uv run pytest tests/test_rebuild_gauge_search_spec.py -v` → **5 passed**.
5. **Lint.** ✅ **Pass.** `uv run ruff check .` → `All checks passed!` (full repo, including the new
   files).

## Why this stays open

Per this run's instructions, chem-g14 closes only when its acceptance criteria are met *and the
test suite passes* — not just "moved the queue." Here, criterion 3 fails: 14-15/45 cases in
`tests/test_eigenstate_reachability_spec.py` still fail after the fix. Root-causing each remaining
failure (one by one, against an unmodified `main` checkout via `git stash` on `reachability.py`
alone) shows **none of them are caused by this fix** — they were previously *masked* by the old
atol always returning `None` (so the downstream code they exercise never ran), and now that the
rebuild succeeds often enough to reach that downstream code, two independent, pre-existing bugs
surface instead:

1. **`build_molecular_hamiltonian`'s own (first, non-rebuilt) free RHF solve is itself
   nondeterministic run-to-run** on square H4 — confirmed directly: 10 fresh builds of the same
   `(atom, conv_tol)` land on two different branches (9/10 at `hf_energy=-1.7825510081` with the
   known residue population at 1.2e-6; 1/10 at `hf_energy=-1.7825511823`, residue at 2.4e-9, below
   the gate's `tol=1e-8` cut). This explains `test_G1_residue_level_identified_warned_and_kept` and,
   conclusively, `test_G9_control_tight_scf_stays_on_the_ag_level` — the latter's own source
   (`tests/test_eigenstate_reachability_spec.py:320-325`) calls `build_molecular_hamiltonian`
   directly and never reaches `_driver_orbitals` at all, so this fix cannot be responsible. **Filed
   as chem-7cx** (already existed as a candidate before this session; updated with this evidence).
2. **`orbital_parity_bounds`/`hf_symmetry_sector`'s square-H4 physics is independently wrong.**
   `test_G2_decision_equals_symmetry_adapted_fci` already fails on unmodified `main` (e.g. a=1.35,
   conv_tol=1e-9: baseline gives a 36-vs-10 "in Ag sector" mismatch against an independent
   symmetry-adapted FCI; with this fix, 20-vs-10 — improved, since the computation now actually
   runs instead of short-circuiting, but still wrong). `test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged[1.2]`
   also fails on baseline (`TypeError: 'NoneType' object is not iterable`, since
   `orbital_parity_bounds` was always `None`); with the fix it fails differently — the rebuilt MOs
   are confirmed to be a *genuine* match to `mh`'s own MOs (residual 1e-6 to 1e-5, deep in this
   module's own noise floor, not a false-positive sign-gauge accept), but one bit's computed `beta`
   comes out ~0.0008 instead of the expected `>= 1.0`. **Filed as chem-3fd.**

Both are confirmed present on baseline `main` (unmodified) too — this fix did not introduce either;
it only stopped silently masking them by making the gap it targets actually close often enough for
the downstream code to run.

### Cross-check against an independent gate file (same session)

Before concluding, I also ran the broader spec suite
(`tests/test_certified_gaps_spec.py tests/test_excited_bounds_spec.py
tests/test_certified_dipole_spec.py tests/test_hf_overlap_certificate_spec.py
tests/test_scf_conv_tol_spec.py`) to check for regressions elsewhere that import `reachability`.
Result: **6 failed, 50 passed**, all 6 in `tests/test_scf_conv_tol_spec.py`
(`test_G1_clean_geometry_does_not_move[1.0/1.2]`,
`test_G2_witness_moves_but_the_eigenvalue_does_not[1.19/1.35]`,
`test_G3_residue_collapses_at_the_tight_tolerance[1.19/1.35]`, the last with
`p0=0.4621025036353699` where `<1e-20` was expected). Traced the call path
(`_measure` → `_dense_hf_projection` / `reachable_eigenpairs` / `exact_reachable_overlap`) and
confirmed **none of it calls `_driver_orbitals`, `hf_symmetry_sector`, or `symmetry_allowed`** —
only `reachable_mask`'s diagnostic warning does, under a broad `try/except` that never changes the
returned mask/measurement. A baseline run of just this file (`git stash` on `reachability.py` alone)
reproduces the same shape of failure (5/14 failed, including `test_G3[1.35]`'s own
`p0` anomaly) with zero changes to this fix. This is the same root cause as (1) above (the original,
non-rebuilt SCF solve's branch nondeterminism), surfacing through a different gate file. Folded into
chem-7cx as further confirming evidence rather than filed as a third bead.

## Additional verification (not required by the acceptance criteria, done for confidence)

- `uv run pytest tests/test_rebuild_gauge_search_spec.py -v` → `5 passed`.
- `uv run ruff check reachability.py tests/test_rebuild_gauge_search_spec.py` → clean, and
  `uv run ruff check .` (full repo) → clean.
- Direct, instrumented reproduction of the sign-gauge search at the a=1.2 witness (5 trials):
  unflipped residuals ranged 4.7e-6 to 1.18e-5; best-over-signs residuals ranged 1.9e-6 to 4.7e-6 —
  both deep in the noise floor, confirming matches found by the search are genuine reconciliations
  of the same MO basis, not false-positive accepts of a different one.

## What I decided not to do, and why

- Did not attempt to fix `build_molecular_hamiltonian`'s own SCF-branch nondeterminism (chem-7cx) or
  `orbital_parity_bounds`/`hf_symmetry_sector`'s square-H4 sector physics (chem-3fd) as part of this
  bead. Both are confirmed pre-existing and independent of `_driver_orbitals`'s rebuild-
  reproducibility (chem-g14's actual, narrow scope); fixing them here would mean fixing problems
  this bead didn't ask me to fix and that would need their own spec/investigation per this repo's
  SDD convention. Filed as their own beads instead, per "work discovered outside this bead's scope
  must be filed as a new bead, not fixed now."
  - Note: both beads were originally filed in the segment before this one; this session's main
    contribution to them is the additional cross-check evidence from `test_scf_conv_tol_spec.py`,
    folded into chem-7cx's notes rather than opening a third bead for the same root cause.
- Did not loosen `_SAME_OPERATOR_ATOL` further, or widen `_ENERGY_MATCH_ATOL`/retry counts, to try to
  force the remaining 14-15 failures closed. They are not reconciliation-gap failures (the rebuild
  itself now matches correctly, per the noise-floor residuals above) — loosening further would only
  risk accepting a false-positive match on some other system, trading a correctness guarantee for a
  cosmetic pass count.
- Did not revise `specs/SPEC_eigenstate_reachability.md` itself (the file whose §9 pre-registered
  this gap). This bead's own spec (`SPEC_rebuild_gauge_search.md`) is the right place to record the
  honest, falsified-then-revised G2 gate; editing the original spec's historical claim seemed out of
  scope and I didn't want to touch a file I wasn't asked to.

## What I could not verify

- Whether the `build_molecular_hamiltonian` branch nondeterminism (chem-7cx) and the
  `orbital_parity_bounds` physics bug (chem-3fd) are platform-specific (this container's PySCF/BLAS
  build) or would reproduce on the host the original `SPEC_eigenstate_reachability.md` baseline
  ("45 passed") was recorded on. Both beads note this as an open question for whoever picks them up.
- Whether `_MAX_SIGN_SEARCH_ORBITALS=10` / `_MAX_REBUILD_RETRIES=5` are the right bounds for systems
  larger than this module's own gate file exercises (max observed need: `k-1=5` at the largest
  gated system) — left as documented, un-exercised headroom per the spec's §7.

## Commands to reproduce

```bash
cd /workspace
uv sync --extra dmrg --extra test --extra dev
uv run ruff check .                                              # -> All checks passed!
uv run pytest tests/test_rebuild_gauge_search_spec.py -v         # -> 5 passed

# Regression check (takes ~2-3 min per run; results flicker run-to-run, see above):
uv run pytest tests/test_eigenstate_reachability_spec.py -q      # -> 14-15 failed, 34-35 passed

# Baseline comparison (reverts only reachability.py):
git stash push -- reachability.py
uv run pytest tests/test_eigenstate_reachability_spec.py -q      # -> 17 failed, 28 passed
git stash pop

# Cross-check gate file (unrelated to this fix, confirms the same upstream nondeterminism):
uv run pytest tests/test_scf_conv_tol_spec.py -q                 # -> ~5-6 failed, ~8-9 passed
```

## Git

Committed on `sandbox/chem-g14` (per this run's git policy, which authorizes committing — but not
pushing or merging — to this branch):

```
<see `git log sandbox/chem-g14` for the commits made alongside this handoff>
```

Tree otherwise clean. `.claude/settings.local.json` (untracked, pre-existing local tool-permission
config, unrelated to this task) deliberately not committed. Per policy, not pushed.

## Bead

Left **open** (`bd update chem-g14 --notes=...`, not closed) because
`tests/test_eigenstate_reachability_spec.py` does not fully pass (criterion 3 above). The fix itself
is implemented, gated, and verified to close chem-g14's actual, narrow claim (rebuild-
reproducibility) without regression; the remaining failures are two separately-filed, pre-existing,
out-of-scope bugs (chem-7cx, chem-3fd) that this bead's own fix does not touch and cannot close.
