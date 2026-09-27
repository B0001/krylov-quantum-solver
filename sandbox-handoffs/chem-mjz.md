# chem-mjz — DMRG stall that regime labels "converged" (Hubbard L=60, U=4)

## Bead (verbatim, `bd show chem-mjz`)

> chem-tjr (PR #49 §10.1): open L=60, U=4, block2's default random MPS: the D=100 stage sat ~13 Ha
> too high and even D=400 was 12.6 mHa high, while every discarded weight was <= 4e-9, so
> `truncation_regime` returned 'converged'. chem-4e9's new checks were designed around a
> 'truncation' ladder (the stalled-stage slope test needs a truncation ladder with >= 3 points);
> its undershoot test runs on converged ladders but won't fire if every stage is high. So this
> stall may still pass. Reproduce it cheaply (small L if possible) and add it as a regime gate; if
> it passes as 'converged', extend the check.

**Acceptance criteria (verbatim):** "The L=60 U=4 stall (or a cheaper reproduction) is classified
'uncontrolled' by truncation_regime alone; all existing regime/extrap/hchain/hubbard gates still
pass."

## What changed

`hybrid_quantum_solver/dmrg_reference.py`:

- New constant `CONVERGED_UNDERSHOOT_FACTOR = 1.0` (next to the existing `UNDERSHOOT_FACTOR =
  10.0`), derived analytically: for `E = E_inf + C/D` with D doubling at every stage (this repo's
  convention), the exact tail past `D_max` telescopes onto `gap_last` itself
  (`gap_last = C/D_n`, since `D_{n-1} = D_n/2`), so `drop/gap_last = 1` for genuine smooth 1/D
  convergence. `UNDERSHOOT_FACTOR=10` was calibrated for the *discarded-weight* axis
  (SPEC_regime_stalled_stage §3) and was never re-derived for the 1/D axis used once the weight
  test already says "converged" — that is the actual gap this bead closes.
- `truncation_regime`'s undershoot check now picks `UNDERSHOOT_FACTOR` on the `"truncation"`
  branch and `CONVERGED_UNDERSHOOT_FACTOR` on the `"converged"` branch (previously always
  `UNDERSHOOT_FACTOR`). No other constant, and no other branch, changed.
- Docstring updated to describe the two-factor check.

`tests/test_hubbard_lieb_wu_spec.py`:

- `test_G7_stalled_ladder_is_labelled_converged_but_fails_the_spread_free_checks` (its own
  docstring: "Pinned so a classifier fix has to update this spec") renamed to
  `test_G7_stalled_ladder_is_now_rejected_by_truncation_regime_alone` and its assertion flipped
  from `== "converged"` to `== "uncontrolled"`, using the exact recorded §10.1 triple, unchanged.

New files:

- `specs/SPEC_regime_converged_undershoot.md` — the spec: goal, derivation, tolerance table,
  rejected alternatives (absolute cap on the fit slope; feeding `stage_dE` into
  `truncation_regime`), and the cheap reproduction.
- `tests/test_regime_converged_undershoot_spec.py` — 7 gates (G1–G4 below).

## Verification commands and their output

**The classifier fix, both real stalls (regenerate anytime, no DMRG needed):**

```
$ uv run python3 -c "
from hybrid_quantum_solver.dmrg_reference import truncation_regime
L60 = [(100, 3.524227915646634e-09, -20.780402942374906),
       (200, 5.688756095253755e-10, -30.024624553972146),
       (400, 1.2913782637611208e-12, -34.04430453150846)]
L50 = [(100, 8.520874863357602e-09, -20.95984405499062),
       (200, 3.640909381391424e-10, -27.53465639784129),
       (400, 1.0316183797650565e-13, -28.320218264676182)]
print(truncation_regime(L60), truncation_regime(L50))
"
uncontrolled uncontrolled
```

Before this change (`UNDERSHOOT_FACTOR` applied to both branches), both printed `converged`
(verified during development, not re-checked after the fix since the old code path no longer
exists — see the diff below and G1/G2 in the new gate file, which pin the *new*, correct output).

**Diff (`git -c safe.directory='*' diff hybrid_quantum_solver/dmrg_reference.py`):** 30
insertions/3 deletions, shown in full at the bottom of this file for a reviewer who wants it
without checking out the branch.

**Gate runs, verbatim, one line per file** (pure gates via `uv run python -m pytest <file> -q`;
DMRG gates via `GATE_JOBS=1 GATE_GLOB='tests/<file>' bash scripts/run_gates.sh`, each in its own
process per the block2/pyscf segfault rule):

| gate file | result | DMRG? |
|---|---|---|
| `tests/test_extrap_regime_spec.py` | `18 passed in 1.72s` | no (unchanged file) |
| `tests/test_regime_stalled_stage_spec.py` | `15 passed in 1.41s` | no (unchanged file) |
| `tests/test_regime_converged_undershoot_spec.py` | `7 passed in 1.33s` | no (new file) |
| `tests/test_hubbard_lieb_wu_spec.py` | `5 passed, 1 warning in 61.66s` | yes (G8 anchor) |
| `tests/test_hubbard_bethe_spec.py` | `4 passed, 4 warnings in 91.44s` | no (PySCF FCI only, unchanged) |
| `tests/test_hchain_largen2_spec.py` | `5 passed in 253.51s (0:04:13)` | yes |
| `tests/test_hchain_tdl_spec.py` | `7 passed in 267.68s (0:04:27)` | yes |

That covers every file matching `regime`, `extrap`, `hchain`, and `hubbard` in `tests/` — the four
gate families named in the acceptance criteria. `uv run ruff check .` → `All checks passed!`.

**Hardware/wall time** (all DMRG runs above and the L=50 reproduction below), per CLAUDE.md's
>1-minute rule: `Linux 6cfb4191db6a 7.0.12-linuxkit x86_64`, `nproc`=8, `MemTotal`≈16 GB.

## The cheap reproduction (L=50)

L=60 (the recorded stall) needed 590 s. L=50, same conditions (open chain, U=4, default random
MPS — no `init_occs` — D=100/200/400, 8 sweeps/stage, 4 threads, 6 GB stack), reproduces the same
failure mode in **206.5 s**:

```
$ uv run python3 -c "
from hybrid_quantum_solver.dmrg_reference import dmrg_energy_extrapolated
from hybrid_quantum_solver.model_hamiltonians import hubbard_chain_integrals
m = hubbard_chain_integrals(50, 4.0, open_chain=True)
r = dmrg_energy_extrapolated(m.h1, m.eri, m.nelec, m.e_core,
                              bond_dims=(100, 200, 400), n_sweeps_per=8,
                              n_threads=4, stack_mem=6*1024**3,
                              scratch='/tmp/repro/.dmrg_tmp_L50')
print(r.per_D); print(r.stage_dE); print(r.regime)
"
per_D    [(100, 8.520874863357602e-09, -20.95984405499062),
          (200, 3.640909381391424e-10, -27.53465639784129),
          (400, 1.0316183797650565e-13, -28.320218264676182)]
stage_dE [1.1970954371986195, 0.4698009917278796, 7.176481631177012e-13]
regime   uncontrolled     # after the fix; was "converged" before it
```

This exact triple is vendored in `tests/test_regime_converged_undershoot_spec.py` (G2,
`L50_REPRODUCTION_PER_D` / `L50_REPRODUCTION_STAGE_DE`) and in `specs/SPEC_regime_converged_undershoot.md`
§5 — a tracked file, not `data/` (which is gitignored). L=10/20/30/40 were tried first and did not
reproduce this specific classifier gap (either converged cleanly, or the old code already flagged
them via the non-variational or stalled-stage checks) — not vendored, since they don't demonstrate
anything the existing gates didn't already cover.

## Pre-registered criteria and pass/fail

| criterion | result |
|---|---|
| L=60 stall (or cheaper repro) classified `"uncontrolled"` by `truncation_regime` alone | **PASS** — both L=60 (recorded) and L=50 (new, cheaper) → `uncontrolled`, using only `(D, dw, E)` triples, no `stage_dE` or other extra input |
| All existing regime/extrap/hchain/hubbard gates still pass | **PASS** — see table above; every file in those four families is green |
| New check is scale-free (not an absolute energy/slope cap) | **PASS** — G3 in the new gate file: the same `drop/gap_last` ratio, energies uniformly scaled by 1e-6, still rejects |
| No false positives on any vendored real ladder | **PASS** — G4: all 6 converged rows in `hchain_tdl_localized_table.csv` and all converged rows in `hubbard_lieb_wu_table.csv` keep `"converged"` |

## What was decided not to do, and why

- **Not adding `stage_dE` as an input to `truncation_regime`.** It's a more direct stall signal
  and `dmrg_energy_extrapolated` already computes it, but the acceptance criterion requires
  classification "by `truncation_regime` alone" from the recorded triples, and the L=60 record's
  `stage_dE` was never captured (it predates that field) — a `stage_dE`-based fix can't reclassify
  the existing pinned fixture. `stage_dE` stays available as an independent, already-wired data
  gate (`hubbard_tdl_analysis.point_check`); this bead doesn't touch that path. Recorded as a valid
  follow-up in SPEC_regime_converged_undershoot.md §7, same as SPEC_hubbard_bethe §10.4 already
  flagged.
- **Not adding an absolute cap on the fit slope `C`.** SPEC_regime_stalled_stage §3 already
  rejected this for the truncation axis (synthetic fixtures in `test_extrap_regime_spec.py` have
  energies deliberately unrelated to weights, implying `C` as large as ~1e6); the same argument
  applies on the 1/D axis, so any absolute cap would also misclassify those fixtures. Kept the fix
  scale-free (`drop/gap_last`), matching the existing style.
- **Not touching `ENERGY_NOISE` or `STAGE_SLOPE_RATIO`.** Both are still pinned by
  `test_regime_stalled_stage_spec.py::test_G5_tolerances_are_pinned`, unaffected by this change.
- **Not vendoring L=10/20/30/40 reproductions.** They don't reproduce this specific gap (see
  above); vendoring them would just be noise.
- **Did not add a BACKLOG.md entry.** SPEC_regime_stalled_stage.md (chem-4e9), the closest
  precedent for this exact kind of classifier-only fix, also has no BACKLOG.md entry — it's a
  narrow bugfix to an existing spec's gate, not a new physics hypothesis.

## What could not be verified / incident during the session

- **A transient hang, self-diagnosed and recovered, not a code issue.** While running the DMRG
  gates, I initially launched `test_hchain_largen2_spec.py`, `test_hchain_tdl_spec.py`, and
  `test_hubbard_lieb_wu_spec.py` as three concurrent background processes to save wall time. All
  three default to the same scratch directory (`dmrg_energy_extrapolated`'s default
  `scratch="./.dmrg_tmp"`), and `test_hchain_largen2_spec.py`'s run hung indefinitely (near-zero
  CPU growth over 30+ minutes, thread state `futex_do_wait`) — almost certainly a scratch-file
  collision between the two DMRG processes that briefly overlapped (`test_hubbard_lieb_wu_spec.py`
  finished cleanly in 63s; `test_hchain_largen2_spec.py` did not). I killed the hung process
  (`kill -9`, confirmed no other work was in flight) and reran it alone; it then passed cleanly in
  254s, in line with its previously recorded 429s. **Lesson for future sessions: never run two
  DMRG-backed spec gates concurrently even briefly, even though `make gates`'s parallelism is
  documented as safe for pyscf/block2 *process* isolation — it does not isolate the *scratch
  directory*, which is a separate, undocumented shared-resource hazard.** This did not affect any
  reported number; the affected process was killed before producing output, and its log
  (`logs/gates/test_hchain_largen2_spec.log`) was empty at kill time.
- `/workspace` was at 99% disk usage (6.2 GB free) for the whole session. Not investigated further
  (out of scope for this bead) and did not appear to cause any observed failure, but is worth
  flagging — a filed-separately concern, not something I fixed or a new bead I'm creating, since I
  didn't reproduce an actual failure caused by it.
- I did not re-verify that the *old*, unfixed code actually returns `"converged"` on the L=50/L=60
  triples after making the fix (the old code path is gone). This was verified during development
  (see Problem Solving in prior session notes) but is not re-checked here as a numbers artifact —
  only the new, correct behavior is pinned in the gates.

## Bead status

Closing `chem-mjz` — acceptance criteria met, evidence above, all named gate families green.

## Git status (left dirty, per instructions — do not commit/push)

```
$ git -c safe.directory='*' status --short
 M hybrid_quantum_solver/dmrg_reference.py
 M tests/test_hubbard_lieb_wu_spec.py
?? specs/SPEC_regime_converged_undershoot.md
?? tests/test_regime_converged_undershoot_spec.py
```

(`.beads/issues.jsonl` and `.claude/settings.local.json` also show as changed/untracked — the
former from `bd update --claim`/`bd close`, unrelated to the code change.)

**Suggested commands for a human to run** (not run here, per git policy):

```
git add hybrid_quantum_solver/dmrg_reference.py tests/test_hubbard_lieb_wu_spec.py \
        specs/SPEC_regime_converged_undershoot.md tests/test_regime_converged_undershoot_spec.py
git commit -m "dmrg_reference: separate undershoot factor for the converged (1/D) regime (chem-mjz)"
git push
bd dolt push
```

## Full diff of `hybrid_quantum_solver/dmrg_reference.py`

```diff
--- a/hybrid_quantum_solver/dmrg_reference.py
+++ b/hybrid_quantum_solver/dmrg_reference.py
@@ -142,6 +142,22 @@ STAGE_SLOPE_RATIO = 1e3
 # bound every gap by ~C*floor ~ 3e-7 Ha, and their recorded drops are 3e-8 to 9e-8 Ha.
 UNDERSHOOT_FACTOR = 10.0
 
+# The 1/D-axis analogue of UNDERSHOOT_FACTOR, for ladders the WEIGHT test alone already calls
+# "converged" (chem-mjz). For E = E_inf + C/D with D doubling at every stage, the exact remaining
+# tail past D_max is gap_last itself: gap_k = C(1/D_k - 1/D_{k+1}) = C/(2*D_k), and the tail
+# C/D_n = 2 * C/(2*D_n) = 2 * gap_{n} ... telescoped from gap_last = C/D_n exactly (D_{n-1} =
+# D_n/2), so drop/gap_last = 1 for a ladder that genuinely follows smooth 1/D convergence.
+# UNDERSHOOT_FACTOR=10 was calibrated for the discarded-weight axis and is far too loose here: the
+# Hubbard L=50 and L=60, U=4 stalls from block2's default random MPS (chem-mjz) both pass the
+# weight floor (every discarded weight <= 4e-9) while a stage is still relaxing out of a metastable
+# charge distribution, giving drop/gap_last = 1.15 (L=60) and 4.2 (L=50) -- both < 10, so
+# UNDERSHOOT_FACTOR let them through as "converged". CONVERGED_UNDERSHOOT_FACTOR = 1 rejects both
+# (drop/gap_last must stay <= 1 + noise/gap_last) while keeping every one of the nine vendored
+# "converged" ladders and the SPEC_extrap_regime synthetic fixtures (drop/gap_last <= 0.5) inside
+# it -- those pass on the ENERGY_NOISE term alone since their gap_last is itself noise-scale, or
+# (the synthetic ones) is a flat 1 mHa/stage regardless of D, unrelated to weight, at 0.5x the bound.
+CONVERGED_UNDERSHOOT_FACTOR = 1.0
+
 # Gate on `regime`, not on `method`. See truncation_regime() for why.
 REGIMES = ("converged", "truncation", "uncontrolled")
 
@@ -215,8 +231,15 @@ def truncation_regime(per_D, *, floor: float = DISCARD_WEIGHT_FLOOR,
       * E(D) rises with D by more than ``energy_noise`` (non-variational);
       * a truncation ladder has a stage gap ``STAGE_SLOPE_RATIO``x larger than another stage's
         slope dE/d(dw) allows (the stalled-stage signature);
-      * the extrapolation falls below E(D_max) by more than ``UNDERSHOOT_FACTOR`` x the last stage
-        gap + ``energy_noise`` -- more than any ladder with a lever arm can support.
+      * the extrapolation falls below E(D_max) by more than the last stage gap times
+        ``UNDERSHOOT_FACTOR`` (truncation ladders) or ``CONVERGED_UNDERSHOOT_FACTOR`` (ladders the
+        weight floor alone already calls converged), plus ``energy_noise`` -- more than any ladder
+        with a lever arm, or genuine 1/D convergence, can support. The floor-only weight test
+        cannot see a stage that stalled in a metastable state with a small but not-yet-converged
+        discarded weight (chem-mjz: Hubbard L=50/60, U=4, default random MPS -- every weight
+        <= 4e-9, but a stage still sat several Ha from the answer); the tighter converged-branch
+        factor catches it because a stalled stage's energy gap is far larger than smooth 1/D
+        convergence with the ladder's own D-doubling would allow.
 
     Order matters: the floor is tested BEFORE weight monotonicity, so a converged ladder whose
     weights wobble by float noise reads as ``"converged"``, not ``"uncontrolled"``.
@@ -237,7 +260,8 @@ def truncation_regime(per_D, *, floor: float = DISCARD_WEIGHT_FLOOR,
     x = dws if regime == "truncation" else 1.0 / Ds
     energy, _ = _linear_extrapolate(x, Es)
     gap_last = max(float(Es[-2] - Es[-1]), 0.0)
-    if Es[-1] - energy > UNDERSHOOT_FACTOR * gap_last + energy_noise:
+    factor = UNDERSHOOT_FACTOR if regime == "truncation" else CONVERGED_UNDERSHOOT_FACTOR
+    if Es[-1] - energy > factor * gap_last + energy_noise:
         return "uncontrolled"
     return regime
```
