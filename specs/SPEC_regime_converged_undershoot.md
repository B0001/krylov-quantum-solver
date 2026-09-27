# SPEC: the converged-regime undershoot check was calibrated for the wrong axis

**Status:** IMPLEMENTED once gates green. Bead `chem-mjz`. Follows SPEC_regime_stalled_stage.md
(bead chem-4e9) and SPEC_hubbard_bethe.md §10.1/§10.4 (bead chem-tjr), which recorded the finding
this spec fixes.

---

## 1. Goal

SPEC_hubbard_bethe §10.1 recorded a genuine DMRG stall that `truncation_regime` mislabels
`"converged"`: an open L=60, U=4 Hubbard chain from block2's **default random MPS** (no
`init_occs`) had every discarded weight ≤ 4e-9 (well under `DISCARD_WEIGHT_FLOOR=1e-8`), so the
weight-only test calls it converged, while the recorded energies were **-20.780 / -30.025 /
-34.044 Ha** at D = 100/200/400 — the D=100 stage sat about 13 Ha above the answer because DMRG was
still moving charge out of a metastable distribution near the Mott gap, not truncating.
SPEC_regime_stalled_stage (chem-4e9) added energy-aware checks on top of the weight test, but its
stalled-stage check only runs on the `"truncation"` branch, and its undershoot check used one
factor (`UNDERSHOOT_FACTOR=10.0`) for both axes. §10.4 flagged this explicitly as unfixed: *"The
classifier itself should learn from [stage_dE]: either take stage convergence as an input, or bound
the undershoot on a converged ladder by ENERGY_NOISE alone."*

**The falsifiable claim:** the existing `(D, dw, E)` triples already contain enough signal to
reject this stall, with no new input (no `stage_dE`, no absolute cap on the fit slope — both
considered and rejected, see §3) — only a tighter, separately-derived undershoot factor on the 1/D
axis.

## 2. Approach

`truncation_regime`'s undershoot check (SPEC_regime_stalled_stage §2.3) compares the fit's
undershoot below `E(D_max)` against `factor · gap_last + ENERGY_NOISE`, where
`gap_last = E(D_{n-1}) - E(D_n)`. It used `UNDERSHOOT_FACTOR = 10.0` unconditionally. That value was
derived from the discarded-weight axis (`drop/gap_last` = 0.07-0.10 on real truncation ladders, ~98x
margin at 10) and never re-derived for the 1/D axis used when the weight test already says
`"converged"`.

**The axes have different natural scales.** On the 1/D axis, if `E = E_inf + C/D` exactly and each
stage doubles D (this repo's convention: 100/200/400/800/...), the tail past `D_n` telescopes
exactly onto the last observed gap:

```
gap_last = E(D_{n-1}) - E(D_n) = C/D_{n-1} - C/D_n = C/D_n         (since D_{n-1} = D_n/2)
tail      = E(D_n) - E_inf     = C/D_n
=> tail / gap_last = 1
```

So a ladder that genuinely converges as 1/D with doubling D should undershoot its last point by
**at most about 1x `gap_last`**, not 10x. `CONVERGED_UNDERSHOOT_FACTOR = 1.0` is used on the 1/D
branch; `UNDERSHOOT_FACTOR = 10.0` is unchanged on the dweight branch. Both real stalls
(L=50 reproduction, §5; L=60 record, §10.1) have `drop/gap_last` of 4.2 and 1.15 respectively —
both comfortably reject at the new 1.0x bound and both would have passed the old 10x bound.

**Rejected: an absolute cap on the fit slope `C`.** SPEC_regime_stalled_stage §3 already rejected
an absolute cap on the truncation-axis slope because `test_extrap_regime_spec.py`'s synthetic
fixtures have energies deliberately unrelated to their weights, implying `C` as large as ~1e6. The
same argument applies here: those same fixtures, reread on the 1/D axis, imply comparably large
`C`, so any absolute cap would also misclassify them. The fix stays scale-free
(`drop/gap_last`), like the rest of §2 of SPEC_regime_stalled_stage.

**Rejected: feeding `stage_dE` into `truncation_regime`.** `dmrg_energy_extrapolated` already
records `stage_dE` (last-sweep energy change per D-stage) as a `dataclass` field, and it is a more
direct stall signal. But the recorded L=60 stall's `stage_dE` was never captured (§10.1 predates
that field), so a `stage_dE`-based fix cannot reclassify the existing recorded triple, and the
chem-mjz acceptance criterion requires the ladder be "classified `uncontrolled` by
`truncation_regime` alone" from its `per_D` triples. `stage_dE` remains a separate, independent
data-gate check (`hubbard_tdl_analysis.point_check`, already wired) — this spec doesn't change that.

## 3. Tolerance derivation and validation

| constant | value | derivation | validated against |
|---|---|---|---|
| `CONVERGED_UNDERSHOOT_FACTOR` | 1.0 | exact tail = gap_last for 1/D convergence with D doubling (§2) | keeps all 9 vendored converged rows in `hchain_tdl_localized_table.csv` and all 30 in `hubbard_lieb_wu_table.csv` converged (their `gap_last` is noise-scale, so `ENERGY_NOISE` alone dominates); keeps every converged fixture in `test_extrap_regime_spec.py` and `test_regime_stalled_stage_spec.py` (`nbn_headline_converged`: drop/gap_last=0.4999; G2 floor and G4 wobble fixtures: same synthetic-E pattern); rejects both real stalls (L=50 reproduction: 4.2x; L=60 record: 1.15x) |

No other constant changes. `ENERGY_NOISE`, `STAGE_SLOPE_RATIO`, and `UNDERSHOOT_FACTOR` are
untouched and still pinned by `test_regime_stalled_stage_spec.py::test_G5_tolerances_are_pinned`.

## 4. Public interface

```
dmrg_reference.CONVERGED_UNDERSHOOT_FACTOR : float               # NEW
dmrg_reference.truncation_regime(per_D, *, floor=..., energy_noise=...) -> str   # now regime-aware on the undershoot check
```

`UNDERSHOOT_FACTOR` keeps its old meaning and value for the `"truncation"` (dweight-axis) branch.

## 5. Reproducing the stall cheaply

L=60 needed 590 s (§10.1). **L=50, open chain, U=4, default random MPS (no `init_occs`), D =
100/200/400, 8 sweeps/stage, 4 threads, 6 GB stack** reproduces the same failure mode in **206.5 s**
on this container (Linux x86_64, 8 vCPU, ~16 GB RAM):

```
per_D:    [(100, 8.520874863357602e-09, -20.95984405499062),
           (200, 3.640909381391424e-10, -27.53465639784129),
           (400, 1.0316183797650565e-13, -28.320218264676182)]
stage_dE: [1.1970954371986195, 0.4698009917278796, 7.176481631177012e-13]
```

Weight-only regime: `"converged"` (max weight 8.5e-9 < floor). Old code
(`UNDERSHOOT_FACTOR=10.0` on both branches): `"converged"` (drop/gap_last = 4.19 < 10). New code
(`CONVERGED_UNDERSHOOT_FACTOR=1.0` on the 1/D branch): `"uncontrolled"` (4.19 > 1). Smaller L (10,
20, 30, 40) either converge cleanly or happen to already fail the existing non-variational/undershoot
checks even under the old factor — L=50 is the cheapest reproduction found that isolates this
specific gap. Command to regenerate: see `tests/test_regime_converged_undershoot_spec.py::G3`
docstring (`hubbard_chain_integrals(50, 4.0, open_chain=True)` into
`dmrg_energy_extrapolated(..., bond_dims=(100,200,400), n_sweeps_per=8)`, no `init_occs`).

## 6. Acceptance gates (`tests/test_regime_converged_undershoot_spec.py`, pure)

- **G1 — the recorded L=60 stall (DEFINITION OF DONE).** The exact §10.1 triple is classified
  `"uncontrolled"` by `truncation_regime` alone.
- **G2 — the cheap L=50 reproduction is classified the same way**, and the recorded numbers
  (`per_D`, `stage_dE`) are pinned so a future change to the driver has to update this spec.
- **G3 — the constant is scale-free, not an absolute cap.** A synthetic fixture with the same
  `drop/gap_last` ratio as the L=50 stall but energies scaled by 1e-6 is still `"uncontrolled"`;
  one scaled so the ratio is under 1.0 is `"converged"`.
- **G4 — no false positives.** Every vendored converged row in both CSV tables
  (`hchain_tdl_localized_table.csv`, `hubbard_lieb_wu_table.csv`) keeps its recorded regime.
- **G5 — existing regressions.** `test_extrap_regime_spec.py` and
  `test_regime_stalled_stage_spec.py` are unmodified and green (checked directly, not re-asserted
  here — see run log in the handoff).

## 7. Out of scope and caveats

- The D-doubling assumption in the tail derivation (§2) is this repo's convention
  (`dmrg_energy_extrapolated`'s default ladder and every vendored table); a ladder with irregular D
  spacing gets a looser or tighter bound than the exact telescoping argument gives, same as
  `UNDERSHOOT_FACTOR` already tolerates on the dweight axis (SPEC_regime_stalled_stage §7).
  Only one data point (`gap_last`) is used, since the fix must work from the ladders it validates
  against, which are all length 3.
- This does not close the general "regime classifier could still learn from `stage_dE`"
  observation in SPEC_hubbard_bethe §10.4 — it closes the specific literal acceptance criterion
  (chem-mjz: classify the recorded stall as uncontrolled using `truncation_regime` alone). A
  `stage_dE`-aware classifier remains a valid, independent follow-up if a future stall evades both
  the weight floor and this undershoot bound.
- Only one new real reproduction (L=50) was run; L=10/20/30/40 were checked and did not reproduce
  this specific classifier gap (see §5), not vendored as they only reproduce cases the old code
  already handled.
