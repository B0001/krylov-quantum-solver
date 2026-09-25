# SPEC: the extrapolation regime cannot see a ramp stage that never converged

**Status:** IMPLEMENTED once gates green. Bead `chem-4e9`. Follows SPEC_extrap_regime.md and the
n=40 finding in SPEC_hchain_largen2.md §11.2.

---

## 1. Goal

`truncation_regime` (`hybrid_quantum_solver/dmrg_reference.py`) classifies a bond-dimension ladder
from its **discarded weights only**. In SPEC_hchain_largen2 §11.2, the n=40 ramp on 100/200/400 had a
D=100 stage (4 sweeps from a random MPS) that stalled **0.71 Ha** above the D=200/400 energies. Its
weights were still small and monotone, so the regime read `"truncation"`. The discarded-weight fit
then landed **0.82 mHa below its own D=400 energy** (stderr 8.0e-4). Gate G4 of that spec checks
`|E_extrap − E(D_max)| < 1e-6` on the vendored table. That is a data gate, not a fix: the
predicate still certifies a stalled ladder.

**The falsifiable claim:** the ladder's own energies are enough to reject a stalled stage. This
works with no reference energy and no new DMRG run, and it rejects none of the repo's recorded
real ladders.

## 2. Approach

`truncation_regime` still runs the SPEC_extrap_regime weight logic. It now also returns
`"uncontrolled"` if any of these holds:

1. **Non-variational.** `E(D_{i+1}) − E(D_i) > ENERGY_NOISE` for some i, so a larger D raised the
   energy. This runs first, for every ladder.
2. **Stalled stage** (truncation ladders with 3 or more points). The linear model `E = E∞ + C·δ`
   says every stage pair has the same slope `C`. For each pair `j`, the slope is bounded above by
   `s_j ≤ (gap_j + ENERGY_NOISE) / Δδ_j`. The ladder is uncontrolled if some other pair `i` has
   `gap_i > STAGE_SLOPE_RATIO · s_j · Δδ_i + ENERGY_NOISE`. The noise term keeps a pair whose gap is
   at noise level from producing a false positive.
3. **Undershoot.** The fit on the ladder's own axis (δ, or 1/D when converged) lies below `E(D_max)`
   by more than `UNDERSHOOT_FACTOR · gap_last + ENERGY_NOISE`, where
   `gap_last = E(D_{n−1}) − E(D_n)`. This also runs on converged ladders.

**The numbers do not change.** `method`, and so the energy and stderr, still come from the weights
alone (`_weight_regime`, the legacy expression). Only `regime` also reads the energies. An
uncontrolled ladder still returns its fit, and the label says not to trust it. The fit is now in a
pure `extrapolate_ladder(per_D)`, which `dmrg_energy_extrapolated` calls. So the same code can
re-classify a recorded ladder without block2.

## 3. Tolerances, and where each comes from

The data are the nine vendored ladders in `specs/hchain_tdl_localized_table.csv`: localized H_n,
n = 8…40, with 3 truncation rows and 6 converged rows.

| constant | value | measured on recorded ladders | margin |
|---|---|---|---|
| `ENERGY_NOISE` | 1e-6 Ha | converged rows move ≤ 5e-10 from D=200 to D=400. The n=40 D=400 energy reproduced to every printed digit, and D=800 moved it by 1e-9. | ≥ 1e3 above noise. 1e2 below the tightest DMRG gate (SPEC_singleramp, 0.1 mHa). Same 1 µHa as SPEC_hchain_largen2 G4. |
| `STAGE_SLOPE_RATIO` | 1e3 | truncation rows n=24/28/32 have pair slopes 18–29 Ha. Within a ladder they differ by ≤ 1.2×. | ~800× |
| `UNDERSHOOT_FACTOR` | 10 | truncation rows have `drop / gap_last` = 0.07–0.10. | ~98× |

- **Why 10 for the undershoot.** For a linear ladder,
  `drop = C·δ_n = gap_last · δ_n / (δ_{n−1} − δ_n)`. A drop above 10·gap_last requires the weights
  to fall by less than 10% over the last stage. With no lever arm left, no extrapolation is
  supportable. On converged (1/D) ladders, the `ENERGY_NOISE` term decides instead. Their weights
  bound every gap by about C·floor ≈ 3e-7 Ha, and their recorded drops are 3e-8 to 9e-8 Ha.
- **Why 1e3 for the slope ratio.** "Far larger than the discarded weight implies" is read as three
  orders of magnitude. **Honest lower bound:** the ratio must also exceed 100, because
  SPEC_extrap_regime G2 pins `[1e-3, 1e-5, 1e-9]` as `"truncation"`. With that gate's synthetic
  energies, the fixture has a slope ratio of about 100. The constraint is disclosed here rather than
  hidden. 1e3 was not tuned to the n=40 case: G1 shows that the undershoot test alone also rejects
  the stall.
- **Not used: an absolute cap on C.** Real slopes are 13–30 Ha, so a cap such as C ≤ 1e3 would be
  physical. But the synthetic fixtures in SPEC_extrap_regime have energies that are unrelated to
  their weights; for example, the G2 floor fixture implies C ≈ 2.5e5. An absolute cap would fail
  them. Only scale-free checks, which compare the ladder with itself, fit alongside those gates
  without editing them.

## 4. Public interface

```
dmrg_reference.ENERGY_NOISE, STAGE_SLOPE_RATIO, UNDERSHOOT_FACTOR : float
dmrg_reference.truncation_regime(per_D, *, floor=..., energy_noise=...) -> str   # now reads E
dmrg_reference.extrapolate_ladder(per_D, *, floor=...) -> ExtrapResult           # NEW, pure
```

## 5. Acceptance gates (`tests/test_regime_stalled_stage_spec.py`, pure)

- **G1 — the stall is uncontrolled (DEFINITION OF DONE).** The fixture uses D = 100/200/400,
  δ = 1e-7 / 2.3e-10 / 1e-12, and E = −20.924 / E400+1e-9 / E400. It is checked first to be
  faithful: the weights-only predicate says `"truncation"`, and the fit reproduces the recorded
  0.82 mHa undershoot and 8.0e-4 stderr. The new regime is `"uncontrolled"`. The gate repeats this
  across five weight scales (δ_min from 3e-8 to 1e-4). It also shows the stage-gap check and the
  undershoot check each reject the stall on their own.
- **G2 — non-variational.** A 10 µHa rise with D gives `uncontrolled`. A 5e-10 rise on a converged
  ladder stays `converged`.
- **G3 — the floor cannot certify a stall.** The same stall with every δ below 1e-8 gives
  `uncontrolled`. The old predicate said `"converged"`.
- **G4 — no false positives.** All nine vendored ladders keep their recorded regime, `method`, and
  energy (to 1e-9, the CSV print precision). The measured margins in §3 are asserted to be at least
  50× on the truncation rows.
- **G5 — the numbers are unchanged.** For the stalled ladder, the energy and stderr are exactly the
  legacy δ-fit. The tolerances are pinned.
- **G6 — the reproduced n=40 row (§8).** The weights-only predicate says `"truncation"`. The new
  regime says `"uncontrolled"`, and `method` stays `dweight`.

## 6. Reproducing the n=40 row

The reproduction repeats the §11.2 settings: localized H₄₀ from `benchmark_hchain_tdl.integrals`,
`protocol="ramp"`, D = 100/200/400, 4 sweeps per stage, 4 threads, 6 GB stack, and no seed. It
takes about 16 minutes on the 4-core container. The result is in §8.

## 7. Out of scope and caveats

- **Two-point ladders** cannot run the stage-gap check. The undershoot check still runs, but on a
  2-point line `drop/gap_last = δ_n/(δ_{n−1} − δ_n)` is exact, so it only catches missing lever arm.
  A stalled D_min stage in a 2-point ladder is only visible if it also breaks variationality.
- **A stall in the D_max stage** is caught only if E(D_max) rises above E(D_{n−1}) (check 1).
- **The constants come from one family.** They were measured on localized H_n (3 truncation rows),
  and the floor's C ≈ 30 came from NbN. They are justified to an order of magnitude, not derived.
- **The fixture weights are reconstructed.** The stalled run's weights were never recorded, only
  its energies, undershoot, and stderr. The G1 fixture is the small, monotone profile that
  reproduces those three numbers. The parametrized scales show the verdict does not depend on it.

## 8. Result (2026-09-25)

**The n=40 stall reproduced**, on the second attempt. The first attempt's process died after about
12 minutes with a flood of `Intel MKL ERROR: Parameter 13 was incorrect on entry to DGEMM` and
produced no row. The cause was not investigated. The second attempt ran for 947 s:

| D | δ | E (Ha) |
|---|---|---|
| 100 | 7.37e-7 | −20.924206023 |
| 200 | 1.68e-9 | −21.634061152 |
| 400 | 2.69e-11 | −21.634061689 |

The D=100 stage stalled 0.71 Ha high, as recorded. The D=400 energy matches the recorded
−21.6340616891 to 1e-10. The δ-fit gives −21.634881430, which is **0.820 mHa below E(D=400)** with
stderr **7.96e-4**. Both match the §11.2 record (0.82 mHa, 8.0e-4). The weights-only predicate
says `"truncation"`. **The new regime says `"uncontrolled"`**, and the triple is vendored as G6.

Each check on its own:

- **Undershoot: rejects with a wide margin.** drop/gap_last = 1527, against a limit of 10.
- **Stage gap: rejects, but only just.** The pair slopes are 9.6e5 and 325 Ha. After the noise
  inflation (§2), the D=100 gap of 0.710 Ha exceeds its allowance of 0.685 Ha by only about 4%.
  The reason is that the D=200 stage was also slightly under-converged: its gap is 5.4e-7, where
  the 200/400/800 rerun has 1.8e-7. That leaves the later-pair slope loose. The undershoot check
  is the robust detector here, and the stage-gap check is a second line.

**Gates run** (each file in its own process, `uv run pytest`):

| gate file | result |
|---|---|
| `tests/test_regime_stalled_stage_spec.py` | 15 passed |
| `tests/test_extrap_regime_spec.py` | 18 passed (unchanged file) |
| `tests/test_hchain_largen2_spec.py` | 5 passed (DMRG, 429 s) |
| `tests/test_hchain_tdl_spec.py` | 7 passed (DMRG, 163 s) |
| `tests/test_singleramp_spec.py` | 3 passed (DMRG, 239 s) |
| `tests/test_dmrg_reference.py` | 2 passed, 1 skipped |
| `tests/test_nbn_dmrg_reference_spec.py` | **3 failed, environmental**: `FileNotFoundError: data/nbn_scf.chk`. That checkpoint is git-ignored and is built from `data/nb_structures/NbN_mp-2634.cif`, which is also not in the repo. The test fails at file load, before any regime code runs, so it cannot pass or fail on this change here. |
