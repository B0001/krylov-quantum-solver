# chem-j71 handoff — does an HF-weighted band-edge gap predict the Krylov depth M the chem-pm3 η-bound needs?

## Status: DONE on branch `sandbox/chem-j71`. Bead closed with evidence below.

## What changed

- `specs/SPEC_pm3_depth_predictor.md` — new spec, pre-registered 2026-10-14, results filled in
  2026-10-15 (§10).
- `scripts/spec_pm3_depth_predictor.py` — new driver: `find_escape_witnesses`, `weighted_gap`,
  `rank_correlation_report`, `noise_sensitivity_sweep`, `main()`. Reuses
  `scripts/spec_pm3_subspace_eta_bound.py`'s `_chain`, `_asym_chain`, `_build_moment_bound_lp`,
  `run_witness`, `WITNESSES` **unchanged** — no edits to that file.
- `tests/test_pm3_depth_predictor_spec.py` — new gates G1 (regression), G2/G3 (recorded, not
  direction-asserted, per the spec's own design since both could honestly resolve either way).
- `specs/BACKLOG.md` — chem-j71 entry flipped `[ ]` → `[x]` with a one-paragraph result summary.

## The two questions and what the evidence says

**(A) Does WG = min(gap_lo, gap_hi)/η track the bound's difficulty better than the overshoot?**

Command: `uv run python scripts/spec_pm3_depth_predictor.py` (full output in `SPEC_pm3_depth_predictor.md`
§10, raw log was `/tmp/j71_main_run.log`, not committed — rerunning reproduces it; wall time ~45 min
including the witness search).

`find_escape_witnesses()` over 9 pre-registered candidate geometries found 7 new escape witnesses
(R=0.85 and R=1.25 never escaped in M ∈ {6,8,10,12,16,20,24} — recorded as no-witness, not dropped
silently). Combined with chem-pm3's original 3, n=10 witnesses got the full `run_witness` sweep
(M_SWEEP=8,12,16,20,24). Result:

**Spearman ρ(WG, b_min) = −0.867 (p=0.001) vs ρ(overshoot, b_min) = −0.758 (p=0.011).**

Under the pre-registered rule (`|ρ_WG| ≥ 0.7` **and** `> |ρ_overshoot|`): **SUPPORTS** the
weighted-gap hypothesis. Full per-witness table in the spec. This is a genuine (if small-n) finding,
not a wash — the honest caveat is n=10 from one molecule family (linear/asym H₆ chains), so treat the
magnitude of the gap between the two ρ's as suggestive, not a proof of functional form (spec §8).

**(B) How much does the bound actually move under realistic shot noise, measured not just bounded?**

Same run, Part B: perturbed `moments[1:24]` of the linear H₆ R=1.1 witness (M=24, cap=1000,
η_exact=0.0378) with iid complex Gaussian noise at `shot_noise_sigma(shots)` for shots ∈
{1e2,...,1e8}, 40 trials each.

**Every single trial at every shot count — including 1e8 shots — produced an unsound bound** (below
the true η, frequently negative). The empirical mean/std track chem-pm3's own analytic worst-case
slope `(2M−1)·cap·σ` to within roughly a factor of 1-2, meaning that "worst case" estimate is close
to the *typical* case, not a loose upper bound. This is a stronger and more concrete statement than
chem-pm3 §7 could make (it only had the analytic slope, never measured). **This is a load-bearing
caveat for the whole moment-bound-LP approach**: as built (`coeff_cap=1000`, `M=24`), the bound is a
noiseless-only certificate with no margin against even extremely optimistic shot counts. Fixing this
(e.g. a smaller cap, a noise-aware LP) is explicitly out of scope here (spec §7) — recorded as the
finding, not patched.

## Gate-run commands and verbatim pass/fail lines

```
$ uv run pytest tests/test_pm3_depth_predictor_spec.py -v
tests/test_pm3_depth_predictor_spec.py::test_G1_new_witnesses_reuse_the_validated_lp PASSED [ 33%]
tests/test_pm3_depth_predictor_spec.py::test_G2_rank_correlation_is_recorded PASSED [ 66%]
tests/test_pm3_depth_predictor_spec.py::test_G3_noise_sweep_well_formed PASSED [100%]
================= 3 passed, 19 warnings in 1218.16s (0:20:18) ==================
```

Regression check that the shared, *unmodified* chem-pm3 primitives still behave as recorded:
```
$ uv run pytest tests/test_pm3_eta_bound_spec.py tests/test_subspace_floor_resolvability_spec.py -q
............................................... [100%]
47 passed, 8 warnings in 1267.18s (0:21:07)
```
(This ran concurrently with the main driver script, competing for the same 8 cores — hence the
inflated wall time relative to a standalone run; still a clean pass.)

```
$ uv run ruff check .
All checks passed!
```

## Pre-registered criteria, verdict per criterion

| Criterion | Pre-registered rule | Outcome |
|---|---|---|
| G1 | every witness: bound ≥ η, non-increasing in M (true pass/fail) | **PASS** — all 10 witnesses |
| G2 | compute ρ(WG,b_min), ρ(overshoot,b_min); "supports" iff `|ρ_WG|≥0.7` and `>|ρ_overshoot|` | computed: −0.867 vs −0.758 → **SUPPORTS** |
| G3 | measure mean/std/[min,max]/unsound-fraction at every shots value; report both analytic comparators; state plainly if a non-negligible unsound fraction appears at realistic shot counts | computed: **100% unsound at every tested shot count 1e2–1e8**, stated plainly (not softened) |

Done = G1–G3 green (per spec, G2/G3 "green" means computed and recorded, since both could honestly
resolve either way) and §10 filled in. All three hold; §10 is filled in `SPEC_pm3_depth_predictor.md`.

## Decided not to do, and why

- **`m_pass` (full required depth) was not recomputed for the 6 killed new witnesses.** §3 of the
  spec originally sketched an extended sweep to find it; the actual driver run passed `m_extended=()`
  to every `run_witness` call for compute-budget reasons (the base search + 10-witness sweep already
  cost ~45 min on this 8-core/16GB container; an extended sweep to M≈60-64 for 6 witnesses would add
  materially more). G2's *formal* acceptance criterion is stated in terms of `b_min`, which was
  computed in full, so this does not weaken the gate — it only means the richer "required depth"
  comparison sketched in §3's prose was not run. This is recorded as a deviation in §10, not hidden.
  Flagging it here rather than filing a new bead: it is a direct extension of this bead's own driver
  (just a different `m_extended` argument), not new scope.
- **Deriving WG from the LP's actual Markov–Krein/Fejér-kernel extremal theory** — explicitly out of
  scope per spec §7; WG is one cheap, falsifiable candidate operationalization, not a derivation.
- **Fixing the noise-sensitivity problem found in Part B** (smaller cap, noise-aware LP, etc.) — out
  of scope per spec §7; Part B only measures how unsound/volatile the existing bound is.
- **Other molecule families** — out of scope per spec §7, same scope chem-pm3 used (linear/asym H₆).
- **Oracle-free band edges** (`E_d`, `β`, the reachable-level support) — unchanged from chem-pm3,
  still oracle inputs.

## What could not be verified

- **Reproducibility of the exact witness-search outcome across runs is not independently
  cross-checked here against a *third* independent run** — `find_escape_witnesses()` was run twice
  (once inside `main()`, once inside the pytest fixture) and gave the same 7 escapes both times, which
  is some evidence of determinism, but `build_molecular_hamiltonian`'s free-RHF branch is known
  (bd memory `chem-7cx-root-cause`) to be non-bit-associative across PySCF's own OpenMP thread count
  on near-degenerate geometries. `find_escape_witnesses` pins `lib.num_threads(1)` around its own
  SCF/integral build (as chem-7cx-root-cause's fix pattern requires), but the subsequent `run_witness`
  calls (reused unchanged from chem-pm3) and Part B's fresh `build_molecular_hamiltonian` call in
  `main()` do **not** pin thread count — consistent with how chem-pm3's own `run_witness` has always
  been run, but it means a rerun on a container with a different ambient thread count could in
  principle shift results slightly. Not re-verified under a different thread count in this run.
- **The §10 raw run log (`/tmp/j71_main_run.log`, `/tmp/j71_gates.log`) is not vendored into the
  repo** — it lived in `/tmp` on the worker container and is not part of this commit. The numbers in
  `SPEC_pm3_depth_predictor.md` §10 were transcribed directly from it; rerunning
  `uv run python scripts/spec_pm3_depth_predictor.py` reproduces the same table (modulo the thread-count
  caveat above).

## Compute / wall-time summary (for the record)

Container: Linux x86_64, 8 vCPU, 16 GB RAM.
- `find_escape_witnesses()` (9 candidates): ~478-1380s depending on CPU contention (478s standalone
  scratch run; slower inside the committed driver/gates runs because they additionally competed with
  a concurrent regression-test run for part of that window).
- Full `scripts/spec_pm3_depth_predictor.py main()`: ~45 min end-to-end (search + 10×`run_witness` +
  noise sweep), run concurrently with the regression pytest suite for roughly the first half.
- `tests/test_pm3_depth_predictor_spec.py` (re-runs the search + sweep from scratch via its own
  module-scoped fixtures): 1218.16s (20m18s) standalone.
- Regression (`test_pm3_eta_bound_spec.py` + `test_subspace_floor_resolvability_spec.py`): 1267.18s,
  run concurrently with the main driver.

## Suggested next commands for the reviewer

```bash
git log --oneline sandbox/chem-j71 -5
git diff main...sandbox/chem-j71 --stat
uv run pytest tests/test_pm3_depth_predictor_spec.py -v
uv run ruff check .
bd show chem-j71
```

No push performed (per the operating rules for this run — the host pushes, a human merges). No
`git` branch switches, merges, or resets performed.
