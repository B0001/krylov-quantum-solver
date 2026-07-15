# SPEC: SQD plumbing — promoting the smoke test, and exercising the verdicts it never triggers

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`smoke_test_sqd_plumbing.py` already exercises `run_nbn_sqd_sweep.py`'s real shipped functions with
real assertions — a genuinely good smoke test — but it is a raw script (`assert`-and-print), never
wired into `make gates`, and it only exercises a FEW of the two verdict classifiers' possible
outcomes: `validate_row` can return `PASS`/`FRAME_ERROR`/`ABOVE_TOL`/`SKIP`, but the smoke test only
ever produces `PASS` and one instance of the disjunction `FRAME_ERROR`-or-`ABOVE_TOL`; and
`analyze_sector_trend` can return `CONVERGED`/`CONVERGING`/`STALLED`/`FLAT_SUBSPACE`/`INSUFFICIENT`,
but the smoke test only ever produces `CONVERGED`. This spec promotes the existing smoke test into a
real gate (unchanged logic, tightened assertions where the original was loose), and adds the missing
verdict coverage — cheaply, via synthetic inputs to the classifiers directly, no new SQD runs needed
for the untested branches. False if a classifier's documented verdict is unreachable given its own
logic, or if the promoted smoke-test assertions don't hold.

## 2. Background and honest framing

- `smoke_test_sqd_plumbing.py` and `run_nbn_sqd_sweep.py` are reused unmodified — no new physics,
  only a promotion (script -> gate) and a coverage extension (classifier verdicts never exercised).
- **What you can claim if the gates pass:** the three original smoke-test scenarios hold as real CI
  gates (H2 full-space SQD==CASCI; the deliberate frame-break simulation is caught; H4's
  samples-per-batch sweep stays variational and the trend analyzer reports `CONVERGED` on the real
  swept data, not just printed and eyeballed); AND every verdict either classifier can return is
  reachable and correctly triggered — `validate_row`'s `SKIP` path and `analyze_sector_trend`'s
  `FLAT_SUBSPACE`/`STALLED`/`CONVERGING`/`INSUFFICIENT` paths, none of which the original smoke test
  or any other test in this repo ever exercises.
- **What you cannot claim:** that these verdicts occur naturally at the scale this repo's other
  tests run at (the H2/H4 systems here always land in `PASS`/`CONVERGED`, per the module's own
  docstring note that "small/sparse systems hit the full determinant space immediately") — the
  untested-branch coverage (G3) is deliberately synthetic, constructing `sector_rows`/energy inputs
  directly rather than waiting for a real SQD run to misbehave.
- **Reference:** CASCI (dense diagonalization, the smoke test's own comparison target) for G1/G2/G4;
  direct inspection of `validate_row`/`analyze_sector_trend`'s own documented threshold constants
  (`VALIDATION_TOL_MHA`, `VARIATIONAL_TOL_MHA`, `CONVERGENCE_TOL_MHA`, `TREND_MIN_IMPROVEMENT`) for
  constructing G3's synthetic boundary cases.

## 3. Approach

Reuse `integrals_for_spin`, `generate_bit_array_uniform`, `run_sqd_for_sector`, `validate_row`,
`analyze_sector_trend` from `run_nbn_sqd_sweep.py` unmodified — identical calls to the original
smoke test for G1/G2/G4, seeded (`rng = np.random.default_rng(0)`) for reproducibility. G3 calls
`validate_row`/`analyze_sector_trend` directly on constructed inputs (a `None`/`NaN` energy for
`SKIP`; hand-built `sector_rows` lists with `subspace_dim` sequences designed against the module's
own published tolerance constants for `FLAT_SUBSPACE`/`STALLED`/`CONVERGING`/`INSUFFICIENT`) — no
SQD execution needed for those four verdicts.

## 4. Public interface

No new library code — this spec adds only test-file assertions around `run_nbn_sqd_sweep.py`'s
existing public functions, reused unchanged (mirroring `smoke_test_sqd_plumbing.py`'s own calls).

## 5. Acceptance criteria (validation gates)

- **G1 — H2 CAS(2,2) full space: SQD reproduces CASCI exactly, `validate_row` says `PASS`.**
  `|SQD - CASCI| < 1e-3 mHa`; `SQD >= CASCI - 1e-6` (variational); `validate_row(SQD, CASCI) ==
  "PASS"` exactly. *Measured: |Delta| = 0.0000 mHa.*
- **G2 — the deliberate frame-break simulation is caught, and the EXACT verdict is pinned.**
  Dropping `e_core` from a correct energy (the smoke test's own "507 Ha class of bug" simulation)
  makes the result land BELOW CASCI beyond `VARIATIONAL_TOL_MHA`, so `validate_row` returns
  `"FRAME_ERROR"` specifically — not just "one of FRAME_ERROR/ABOVE_TOL" (the original smoke test's
  looser check; tightened here now that the sign is known: `e_core > 0` on this system, so dropping
  it makes the energy more negative). *Measured: e_core=0.715104 Ha, broken energy well below
  CASCI.*
- **G3 — THE FINDING (definition of done): every verdict either classifier can return is reachable,
  not just the ones real SQD runs happen to produce.** `validate_row(None, CASCI)` and
  `validate_row(float("nan"), CASCI)` both return `"SKIP"`; synthetic `sector_rows` constructed
  against the module's own tolerance constants trigger `"FLAT_SUBSPACE"` (subspace never grows),
  `"STALLED"` (subspace grows, delta improves by less than `TREND_MIN_IMPROVEMENT`),
  `"CONVERGING"` (improves by more, but the endpoint is still above `CONVERGENCE_TOL_MHA`), and
  `"INSUFFICIENT"` (fewer than 2 finite points) — the full five-verdict space of
  `analyze_sector_trend`, only one of which (`CONVERGED`) any real run in this repo has ever
  produced.
- **G4 — H4 CAS(4,4) samples_per_batch sweep: variational at every point, trend analyzer asserted
  (not just printed).** At `samples_per_batch in (10, 30, 80, 200)`: every SQD result is
  `>= CASCI - 1e-6`; `analyze_sector_trend(rows) == "CONVERGED"` on the real swept data — the
  original smoke test computed and printed this value but never asserted on it.

> Definition of done: **G3**. G1/G2/G4 promote what the smoke test already checked (tightened where
> loose); G3 is the new coverage — verdicts that exist in the code but had never been triggered by
> anything in this repo.

## 6. Implementation plan (test-first)

1. Write `tests/test_sqd_plumbing_verdict_coverage_spec.py` encoding G1-G4 (RED in the sense these
   checks are new as pytest gates, even though the underlying calls mirror the existing smoke test).
2. No changes to `run_nbn_sqd_sweep.py` or `smoke_test_sqd_plumbing.py` — every gate calls existing
   public functions.
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- `smoke_test_sqd_plumbing.py` itself is left in place as a standalone script (not deleted or
  redirected to import the new test file) — this spec adds a parallel, gated version rather than
  refactoring the original.
- Larger/Nb3-scale systems `run_nbn_sqd_sweep.py` is ultimately for (this spec stays at the smoke
  test's own H2/H4 scale).
- SQD's own accuracy/convergence behavior beyond what the smoke test already checks (covered
  elsewhere: `SPEC_skqd.md` and related specs already validate SQD-adjacent methods).

## 8. Caveats and risks

- **R1 — G3's synthetic inputs are constructed to hit specific branches given the CURRENT threshold
  constants** (`VALIDATION_TOL_MHA=50.0`, `VARIATIONAL_TOL_MHA=1.0`, `CONVERGENCE_TOL_MHA=10.0`,
  `TREND_MIN_IMPROVEMENT=0.20`); if those constants change, G3's synthetic boundary values would need
  re-deriving (the test reads them from the module rather than hardcoding, where practical, to
  reduce this risk).
- Honest limitation: two systems (H2, H4), matching the original smoke test's own scope exactly.

## 9. Deliverables

- `tests/test_sqd_plumbing_verdict_coverage_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
