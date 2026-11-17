# chem-3ds handoff

## Disposition: fixed

Both leftovers from PR #62 are real and both are now fixed. One claim in the bead text ("odmd.py
**and** trotter_odmd.py set ref = w_eig[0]") turned out to be half stale: `trotter_odmd.py` had
already been corrected (in `d3d28f0`, chem-obf) to use `reach[0]`; only `odmd.py` still had the
`w_eig[0]` bug. Checked that before touching anything (`git log -- odmd.py trotter_odmd.py`,
`grep -rn "w_eig\[0\]"`), so I didn't "fix" a site that wasn't broken.

## What I changed

1. **`odmd.py:67`** — `ODMDProblem.ref` was `float(w_eig[0] - mu)`, the *global* dense-eigh
   ground eigenvalue. `mu`/`width` two lines above are already computed from `reach`
   (`w_eig[reachable_mask(...)]` — the HF-population-filtered band), so `ref` was inconsistent
   with its own frame: on a geometry where the true global ground sits in a spin/spatial sector
   HF cannot reach (the sq-H4 family this repo's reachability work — chem-obf/chem-g14/chem-7cx —
   is all about), `ref` would silently reference an energy ODMD's own HF-overlap signal can never
   see. Changed to `float(reach[0] - mu)`, matching `trotter_odmd.py`'s existing (already-correct)
   convention. `reach` comes from the same ascending `eigh` call `w_eig` does, so `reach[0] ==
   reach.min()`.
   - Updated the `ODMDProblem.ref` docstring/comment (said "exact (FCI) ground energy"; now says
     "exact HF-reachable ground energy ... the global FCI ground only when HF reaches it") and
     the matching line in `trotter_odmd.py`'s `TrotterODMDProblem` docstring, which had the same
     imprecise wording even though its code was already correct.
   - **No observable effect on any currently-gated system** (H2, H4 chain, N2 CAS(6,6)): all have
     HF reaching the true ground state, so `reach[0] == w_eig[0]` there and `ref` is numerically
     identical before/after. The fix only matters for geometries where HF cannot reach the global
     ground — out of this bead's scope to construct a new regression case for (would duplicate
     `test_eigenstate_reachability_spec.py`'s sq-H4 witnesses); flagged as follow-up below.

2. **`rodeo.py`** — `overlap_tol` was a dead parameter on `reference_spectrum` (:49): the function
   never read it, and the two call sites that needed the HF-reachability cutoff
   (`rodeo_ground_energy:71`, the `__main__` demo:87) each re-hardcoded the literal `1e-8`
   instead of using it. Fixed by making `overlap_tol` the single place the cutoff is decided:
   - `reference_spectrum` now clamps `overlaps` at/below `overlap_tol` to exactly `0.0`
     (`np.where(overlaps > overlap_tol, overlaps, 0.0)`), so its own parameter is load-bearing.
   - `rodeo_ground_energy` gained an `overlap_tol: float = 1e-8` parameter, forwarded to
     `reference_spectrum`, and its threshold line is now `w[overlaps > 0].min()` — reading the
     already-clamped array instead of re-guessing the cutoff with a second literal.
   - The `__main__` demo's duplicate `w[ov > 1e-8].min()` likewise became `w[ov > 0].min()`.
   - Updated `specs/SPEC_rodeo.md` §4 to document the new `overlap_tol` kwarg on
     `rodeo_ground_energy` (the spec's own public-interface line was out of sync with the fix).
   - Did **not** touch `tests/test_rodeo_spec.py::test_G3`'s own independent `ov > 1e-8` literal
     — that's the test's own pinned acceptance threshold, not one of the two call sites the bead
     named (`:71`/`:87`), and changing a gate's literal without a reason is exactly what the
     pre-registration rule in the orchestrator prompt says not to do.

## Commands that justify each claim

```
$ git log --oneline -- odmd.py trotter_odmd.py rodeo.py | head -10
d3d28f0 chem-obf: references back on the population cut; symmetry decision demoted to a diagnostic
...
$ git show d3d28f0 -- odmd.py trotter_odmd.py   # confirms trotter_odmd.py already used reach[0]
$ grep -rn "w_eig\[0\]" --include=*.py .        # confirms odmd.py was the only live offender
odmd.py:67: ... ref=float(w_eig[0] - mu))
```

## Gate-run lines (verbatim)

Ran as three separate `uv run pytest -v` invocations (all PySCF/qiskit, no block2 — no process
isolation needed per `tests/`'s own grouping). Host: `Linux 7.0.14-linuxkit x86_64`, 8 vCPUs,
16 GB RAM (`uname -a`, `nproc`, `/proc/meminfo`) — recorded because the longer run exceeded 1
minute by a wide margin.

```
$ uv run pytest tests/test_rodeo_spec.py tests/test_odmd_spec.py tests/test_trotter_odmd_spec.py -v
tests/test_rodeo_spec.py::test_G1_ground_recovery PASSED
tests/test_rodeo_spec.py::test_G2_filter_sharpens_with_cycles PASSED
tests/test_rodeo_spec.py::test_G3_off_resonance_suppression PASSED
tests/test_rodeo_spec.py::test_G4_peak_height_is_reference_overlap PASSED
tests/test_odmd_spec.py::test_G1_overlaps_alone_recover_fci PASSED
tests/test_odmd_spec.py::test_G2_depth_convergence_and_no_variational_floor PASSED
tests/test_odmd_spec.py::test_G3_matched_budget_beats_kqd PASSED
tests/test_odmd_spec.py::test_G4_svd_truncation_is_the_mechanism PASSED
tests/test_trotter_odmd_spec.py::test_G1_step_is_genuinely_trotterized PASSED
tests/test_trotter_odmd_spec.py::test_G2_dmd_adds_nothing PASSED
tests/test_trotter_odmd_spec.py::test_G3_dt_squared_eigenphase_law PASSED
tests/test_trotter_odmd_spec.py::test_G4_richardson_removes_the_bias PASSED
================= 12 passed, 15 warnings in 1026.35s (0:17:06) =================
```

```
$ uv run pytest tests/test_odmd_excited_spec.py tests/test_adaptive_shots_spec.py \
    tests/test_odmd_uq_spec.py tests/test_trotter_restorer_spec.py -v
# every consumer of ODMDProblem.ref/build_odmd_problem besides the three files above
tests/test_trotter_odmd_spec.py::test_G1..G4 PASSED            (re-collected, see note below)
tests/test_odmd_excited_spec.py::test_G1_same_signal_recovers_the_spectrum PASSED
tests/test_odmd_excited_spec.py::test_G2_depth_is_the_excited_resource PASSED
tests/test_odmd_excited_spec.py::test_G3_noisy_gap_beats_qksd_at_matched_noise PASSED
tests/test_odmd_excited_spec.py::test_G4_threshold_mechanism_and_visibility_law PASSED
tests/test_adaptive_shots_spec.py::test_G1_strict_budget_preservation PASSED
tests/test_adaptive_shots_spec.py::test_G2_uniform_is_optimal_on_an_undamped_signal PASSED
tests/test_adaptive_shots_spec.py::test_G2b_geometric_whitening_preserves_eigenphases PASSED
tests/test_adaptive_shots_spec.py::test_G3_optimal_decay_tracks_the_damping_rate PASSED
tests/test_adaptive_shots_spec.py::test_G4_the_win_is_1p3x_in_error_and_zero_in_variance PASSED
tests/test_odmd_uq_spec.py::test_G1_union_coverage_meets_nominal_everywhere PASSED
tests/test_odmd_uq_spec.py::test_G2_width_is_informative PASSED
tests/test_odmd_uq_spec.py::test_G3_union_is_load_bearing PASSED
tests/test_odmd_uq_spec.py::test_G4_bias_is_invisible_the_boundary PASSED
tests/test_trotter_restorer_spec.py::test_G1_there_is_no_variational_floor_leak PASSED
tests/test_trotter_restorer_spec.py::test_G1b_the_extrapolant_is_what_breaks_the_floor PASSED
tests/test_trotter_restorer_spec.py::test_G2_signal_extrapolation_helps_then_degenerates PASSED
tests/test_trotter_restorer_spec.py::test_G3_no_lehmann_bracket_exists_on_this_path PASSED
tests/test_trotter_restorer_spec.py::test_G4_neither_extrapolation_dominates_the_other PASSED
tests/test_trotter_restorer_spec.py::test_G5_the_refinement_parameter_is_reps_not_dt PASSED
================= 23 passed, 14 warnings in 1372.93s (0:22:52) =================
```
(Note: this second invocation actually ran `test_trotter_odmd_spec.py` + the four files above —
I mis-specified the file list on the earlier `nohup` attempt that died silently when its parent
shell context ended; re-ran it correctly as a tracked background task and it is the 23-passed
line above. No test was skipped.)

```
$ uv run ruff check rodeo.py odmd.py trotter_odmd.py
All checks passed!
$ uv run ruff check .      # full-repo lint, confirms no collateral damage
All checks passed!
```

**Total: 35 of 35 tests across every file that imports `rodeo`, `odmd`, or `trotter_odmd` pass.**
I did not run the full `make gates` (≈90 spec files, each its own process — far beyond this
bead's scope and the compute envelope) since nothing outside these three producer files and
their direct consumers was touched.

## Pre-registered criteria vs. outcome

| Criterion (from the bead) | Outcome |
|---|---|
| `rodeo.py:71/87`'s literal `1e-8` replaced by the `overlap_tol` parameter | **Done** — single source of truth now lives in `reference_spectrum`; both call sites read it. |
| `odmd.py`'s `ref = w_eig[0]` fixed to the reachable ground | **Done**. |
| `trotter_odmd.py`'s `ref = w_eig[0]` fixed to the reachable ground | **Not needed** — already `reach[0]` since chem-obf (`d3d28f0`); confirmed, not re-fixed. Docstring wording aligned for consistency instead. |
| "Neither is a symmetry-decision site" — fix or document | **Fixed** (both sites); no symmetry-decision logic was added or implied — `reachable_mask`'s existing population cut is reused as-is. |
| Existing gates stay green | **Done** — 35/35 across all producer + consumer spec files, `ruff check .` clean. |

## What I decided not to do, and why

- **Did not add a new regression test witnessing the `w_eig[0]` vs `reach[0]` divergence on a
  symmetry-forbidden geometry** (e.g. square H4). This bead's scope, per its own text, is the two
  named leftover sites — not re-deriving `SPEC_odmd.md`'s acceptance gates. `H2`/`H4
  chain`/`N2 CAS(6,6)` (the only systems `test_odmd_spec.py`/`test_odmd_excited_spec.py` gate on)
  all have `reach[0] == w_eig[0]`, so no existing gate can distinguish the two formulas — the bug
  was latent, not gated, which is exactly why it survived PR #62. Filing a follow-up bead for a
  square-H4 ODMD witness (mirroring `test_eigenstate_reachability_spec.py`'s G2/G9) rather than
  doing it here, since it's new spec work, not a leftover fix.
- **Did not change `tests/test_rodeo_spec.py::test_G3`'s own `ov > 1e-8` literal.** It's the
  test's pre-registered acceptance threshold, not one of the two call sites the bead named, and
  the orchestrator's pre-registration rule says don't touch gate literals without a measured
  reason.
- **Did not run the full `make gates` suite** (~90 spec files × its own process) — would take
  well over an hour and touches code this bead didn't change. Ran every file that actually
  imports the three touched modules instead (35 tests, listed above).

## What I could not verify

- **No numeric witness that the `odmd.py` fix changes behavior on any real molecule** — by
  construction, every system in this repo's existing ODMD/Trotter-ODMD/excited/UQ/adaptive-shots/
  restorer test suites has `reach[0] == w_eig[0]`, so I cannot show a concrete `prob.ref`
  before/after diff from the current test fixtures alone. The fix is verified *by reading* (mu/
  width/tau all already derive from `reach`, so `ref` deriving from `w_eig[0]` instead was an
  internal inconsistency, not something I can currently force to produce a different float in
  this repo's gated systems) rather than by a numeric regression. Flagged above as a follow-up
  bead rather than claimed as measured.

## Git status at handoff

Dirty, ready to commit on `sandbox/chem-3ds` (cut from `main`). Changed files:
`odmd.py`, `rodeo.py`, `trotter_odmd.py`, `specs/SPEC_rodeo.md`, `.beads/issues.jsonl` (bead
claim/close export), `sandbox-handoffs/chem-3ds.md` (this file).

Proposed commit (not run — git policy for this session is "leave dirty, human commits"):

```
git add odmd.py rodeo.py trotter_odmd.py specs/SPEC_rodeo.md .beads/issues.jsonl sandbox-handoffs/chem-3ds.md
git commit -m "chem-3ds: wire overlap_tol through rodeo.py; odmd.py ref uses the reachable ground"
```
