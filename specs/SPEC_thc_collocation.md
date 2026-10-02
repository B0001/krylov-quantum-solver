# SPEC: Scoring THC collocation strategies with the native qubitization λ (revisits SPEC_thc_lambda G4)

**Status:** CLOSED — gates G1–G4 PASS (`tests/test_thc_collocation_spec.py`, 2026-10-01, 4 passed in
66.5s). **Finding:** Kill A does not fire — `fit_thc`'s nonlinear collocation lands ≈60–69× below
random collocation in λ, comfortably past the ≥5× bar. Kill B fires at most (not all) seeds —
nonlinear collocation beats `df_lambda` at 4/5 pinned seeds (norb=6) and at the one norb=7 point
checked — but §8 R1's caveat holds exactly as anticipated: every "beats `df_lambda`" run also fails
G2's `<1e-6` reconstruction precondition by 4–5 orders of magnitude, so this is read as "an
error-optimal THC fit can land anywhere in λ-space, including below `df_lambda`, while reproducing
the wrong Hamiltonian" — not as "the ISDF/optimized-collocation route is unnecessary." See §9.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

`SPEC_thc_lambda` locked λ_THC ≈ 62× λ_DF as a finding and put *optimized* (ISDF) collocation out
of scope. But `lambda_ladder.fit_thc` — a nonlinear Levenberg–Marquardt THC fit that predates that
spec — has only ever been scored with the brute-force Pauli λ (`lambda_and_terms`, feasible only to
~4 orbitals and a *different 1-norm convention* than `thc_lambda`/`df_lambda`), never with the
spec's own native `thc_lambda`. This spec scores three collocation strategies — random
(`tensor_hypercontraction`), a new deterministic non-random baseline
(`pair_indicator_collocation`), and `fit_thc`'s nonlinear fit — all under the native `thc_lambda`,
at matched THC rank `M = norb(norb+1)/2`, against `df_lambda` as the reference. Two falsifiable
directions, either of which is an accepted outcome:

- **Kill A (`fit_thc` doesn't actually escape the 62× penalty):** if nonlinear collocation is not
  at least 5× below random collocation's λ.
- **Kill B (`SPEC_thc_lambda` G4 must be revised):** if nonlinear collocation's λ beats `df_lambda`
  outright.

## 2. Background and honest framing

- **What you can claim if the gates pass.** A precise, reproducible measurement of where three
  concrete collocation strategies land in λ-space at matched rank and (where achievable) matched
  reconstruction fidelity, closing the "never scored with the native λ" gap in `fit_thc`.
- **What you cannot claim.** Global optimality of `fit_thc`'s fit (it is a stochastic, local
  Levenberg–Marquardt search — pin the seed, report the spread across restarts/seeds, never claim
  "the best possible nonlinear THC"). Nor does this spec attempt ISDF/optimized collocation — the
  new `pair_indicator_collocation` baseline is a **fixed combinatorial formula with no search and
  no RNG**, deliberately *not* the research-grade optimized route `SPEC_thc_lambda` §7 excludes.
- **A precondition surprise, recorded rather than hidden.** The check as originally posed
  (`specs/BACKLOG.md`) requires "matched reconstruction error < 1e-6" for all three collocations.
  That precondition holds trivially for the two *linear* methods (random and pair-indicator: both
  solve `zeta` by least squares given a full-rank `chi`, exact to round-off by construction) but —
  see §5 G2 — **does not hold for `fit_thc`** within a diligent compute budget. That gap is itself
  load-bearing evidence for the caveat this spec's parent bead already anticipated: `fit_thc`
  minimizes reconstruction error, not λ, and (as measured here) does not even reliably minimize
  reconstruction error to the precision needed for a fully apples-to-apples λ comparison.

## 3. Approach

- **System.** LiH/STO-3G, full space: `norb=6` (the CI-gate cap this bead's cost note requires),
  `M = thc_rank(6) = 21`. DF's own full rank on this system is also 21 (coincidence of
  `norb(norb+1)/2` being the pair-space dimension and this system having no rank deficiency) —
  reported, not assumed.
- **Random collocation** (existing, unchanged): `thc_factorization.tensor_hypercontraction(eri,
  norb, n_thc=M, seed=0)`.
- **Structured (deterministic, non-random) collocation** (new):
  `thc_factorization.pair_indicator_collocation(eri, norb)` — one collocation index per symmetric
  orbital pair (`chi^{(i,i)}=e_i`, `chi^{(i,j)}=e_i+e_j`), `zeta` by the same linear least squares.
  Provably full column rank by construction (not by numerical luck) — see the module docstring for
  why a smooth 1-parameter deterministic curve (Fourier/Chebyshev/Vandermonde collocation, tried
  first and rejected) tops out at rank `2*norb-1 < M` and cannot reconstruct exactly.
- **Nonlinear collocation** (existing, extended): `lambda_ladder.fit_thc(eri, norb, M, seed=0,
  return_factors=True)` now returns the winning restart's `(X, Z)` so it can be scored by
  `thc_lambda` instead of only `lambda_and_terms`.
- **Reference.** `df_factorization.df_lambda` on the DF leaves of the same system, plus
  `thc_factorization.thc_lambda` for all three THC collocations (shared 1-norm convention,
  validated against `df_lambda` already by `SPEC_thc_lambda` G2).

## 4. Public interface

```
thc_factorization.zeta_from_collocation(eri, norb, chi)     -> zeta   (shared LS solve, refactor)
thc_factorization.pair_indicator_collocation(eri, norb)      -> (chi, zeta)
lambda_ladder.fit_thc(..., max_nfev=4000, return_factors=False) -> (X, Z) if return_factors else eri
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_thc_collocation_spec.py`. PySCF/NumPy/SciPy only, `norb=6` (fast: seconds to
low tens-of-seconds per gate at CI defaults — restarts=4, max_nfev=4000, matching `fit_thc`'s
existing defaults).

- **G1 — random and structured collocation reconstruct exactly at matched rank.** At
  `M = thc_rank(norb) = 21`, both `tensor_hypercontraction` and `pair_indicator_collocation` give
  `‖ERI_recon − ERI‖ < 1e-9` (round-off), confirming the "matched rank ⇒ matched reconstruction"
  half of the check is achievable in closed form for linear methods.
- **G2 — recorded precondition failure for `fit_thc` (definition of the honest boundary).** At the
  SAME matched rank `M=21` and `fit_thc`'s own fast defaults (`restarts=4, max_nfev=4000, seed=0`),
  reconstruction error is **> 1e-3** (i.e. it does *not* meet the spec's own `<1e-6` precondition).
  This is a recorded finding, not a bug: dies (would need revision) only if `fit_thc` suddenly *does*
  reach `<1e-6` — in which case the honest-boundary claim below is wrong and must be revised.
- **G3 — λ formula sanity anchors.** `thc_lambda` on the exact structured/random collocations is
  positive, finite, and (structured) within 2× of `df_lambda` — a sanity bound, since
  `pair_indicator_collocation` is not claimed to equal `df_lambda` exactly (unlike
  `thc_from_df`, which already does per `SPEC_thc_lambda` G2).
- **G4 — the two killable directions, evaluated and recorded explicitly (not required to pass in
  a fixed direction).** Compute `ratio_5x = lambda(fit_thc) / lambda(random)` and
  `beats_df = lambda(fit_thc) < df_lambda`. The test asserts these are computed and finite, and
  separately **records** (via the test's own docstring/assert message, and in this spec) which of
  the following occurred at `seed=0`:
    - Kill A (`ratio_5x > 0.2`, i.e. NOT ≥5× below random) — the 62× penalty is not merely
      "nobody optimized the points" if this holds despite `fit_thc` failing G2's precondition; or
    - Kill B (`beats_df` is `True`) — `SPEC_thc_lambda` G4 must be revised.
  Both can be true at once (they are not mutually exclusive): `fit_thc` can simultaneously beat
  random by ≥5× **and** beat `df_lambda` outright, while still failing the reconstruction
  precondition (G2) — which is exactly what this spec measures. See §8 for why "beats λ while
  failing reconstruction" is not read as a clean win for `fit_thc`.

> Definition of done: **G1 + G2** (the reconstruction-fidelity boundary, whichever way it falls,
> is the honest content). G3/G4 report the λ comparison for completeness and are evaluative, not
> pass/fail in a single fixed direction — see the caveat in §8 about reading G4's numbers.

## 6. Implementation plan (test-first)

1. Add `thc_factorization.zeta_from_collocation` (refactor out of `tensor_hypercontraction`, no
   behavior change) and `thc_factorization.pair_indicator_collocation`.
2. Add `lambda_ladder.fit_thc(..., max_nfev=4000, return_factors=False)`, backward compatible.
3. Write `tests/test_thc_collocation_spec.py` encoding G1–G4.
4. Run alone (own process, PySCF/NumPy/SciPy, no block2) via `uv run pytest -q
   tests/test_thc_collocation_spec.py`, then through `scripts/run_gates.sh`.

## 7. Out of scope

- ISDF / optimized collocation (unchanged from `SPEC_thc_lambda` §7) — `pair_indicator_collocation`
  is explicitly NOT this: it is a fixed formula with no search, included only as a non-random
  contrast point to random collocation.
- Diagnosing *why* `fit_thc`'s LM optimizer plateaus short of exact reconstruction (a follow-up, if
  ever wanted, is a new bead — this spec measures the fact, not the cause).
- norb=7 (H₂O full STO-3G, the size this bead's cost note used for the "~600 parameters" estimate)
  is run once outside the CI gate and reported in the PR/handoff, not gated.

## 8. Caveats and risks

- **R1 — λ measured on an imperfect reconstruction is not the same claim as λ of the true THC.**
  If `fit_thc`'s best-of-restarts λ beats `df_lambda` (Kill B) while its reconstruction error is
  ~1e4–1e5× the `<1e-6` bar, that is **not** read here as "nonlinear collocation solves the λ
  problem" — a fit that reproduces the Hamiltonian poorly can trivially have a small 1-norm (e.g.
  the identically-zero operator has λ=0 and infinite reconstruction error). The honest reading is
  reported explicitly in the results, not smoothed into a single verdict.
- **R2 — stochastic optimizer.** `fit_thc` is seeded (`seed=0` pinned) but never claimed globally
  optimal from one run; §5 G2's reconstruction-error finding is corroborated by a larger,
  non-CI-gated sweep (12 seeds × 8000 LM evaluations, ~4.5 CPU-minutes) recorded in the PR/handoff,
  not re-run in the fast CI gate.
- Honest limitation: one small molecule (LiH/STO-3G, norb=6), matching this bead's CI-cost cap; not
  a claim about scaling behavior at dozens of orbitals (where THC's asymptotic case for existing).

## 9. Results (chem-5oj, 2026-10-01)

Gate run: `uv run pytest -q tests/test_thc_collocation_spec.py -v` → **4 passed in 66.53s**
(LiH/STO-3G, norb=6, M=21). `tests/test_thc_lambda_spec.py` (the parent spec's own gates) re-run
alongside and still **4 passed in 28.07s** — this spec's new code does not regress the anchor.

**CI-gate point (seed=0, restarts=4, max_nfev=4000):** `df_lambda=15.4846`,
`random_lambda=1001.6941` (recon err 1.4e-13), `fit_thc`: recon err ≈1.6–1.8e-2 (two runs of the
identical seeded call differed at the 2nd significant figure — BLAS-threaded floating-point
summation order is not bit-reproducible across runs; the qualitative finding is unaffected),
`lambda≈14.62–14.65` → `ratio_to_random≈0.0146` (**≈68× below random ⇒ Kill A does not fire**),
`beats_df=True` (**Kill B fires at this seed**).

**Corroboration sweep, `scripts/thc_collocation_sweep.py`** (not CI-gated; run once this session,
reported here per §8 R2):

1. *12-seed × 8000-LM-eval (restarts=1) reconstruction-error sweep*, LiH norb=6 M=21, 259.0s total:
   errors ranged **3.82e-3 to 2.02e-2** (mean 1.27e-2) across all 12 seeds — **every seed stays
   > 1e-6**, i.e. G2's precondition failure is not a seed=0 fluke.
2. *5-seed beats_df/ratio_to_random spot check* at CI-gate defaults (restarts=4, max_nfev=4000):
   `beats_df=True` at seeds 0, 2, 3, 4 and **False at seed 1** (λ=16.63 > df_lambda=15.48) —
   **4/5, not 5/5**: Kill B is a seed-dependent, reproducible-but-not-universal finding, not a law.
   `ratio_to_random` stayed in **[0.0146, 0.0166]** (≥5× below random) at all 5 seeds — Kill A never
   fires in this spot check either.
3. *Single, not-gated norb=7 (H₂O/STO-3G full space) point*, this bead's cited "~600 LM parameters"
   case (`n_params=602`, confirmed): `df_lambda=86.96`, `random_lambda=5377.30`, `fit_thc`
   (175.2s): recon err **8.95e-2** (even further from `<1e-6` than norb=6), `lambda=85.11` →
   `ratio_to_random=0.0158` (≈63× below random, Kill A does not fire) and `beats_df=True` (Kill B
   fires) — the pattern holds qualitatively one orbital larger, still with the same reconstruction
   caveat.

**Reading (per §8 R1):** `fit_thc` is reliably λ-small relative to random collocation (Kill A never
fires, at any seed or system size checked) but is *not* reliably below `df_lambda` (Kill B fires at
most but not all seeds), and every instance of "beats `df_lambda`" is paired with a reconstruction
error 4–5 orders of magnitude above the `<1e-6` bar the comparison was supposed to be matched at.
The honest finding is not "nonlinear collocation solves the 62× penalty" — it is "an
error-optimal THC fit is λ-small but not λ-controlled, and its apparent wins over `df_lambda` come
from fitting a different (wrong) Hamiltonian, not from a validated small-λ reconstruction of the
real one." This does not revise `SPEC_thc_lambda` G4 (which is about the *exact* random-collocation
THC, re-confirmed unchanged here at ≈62–65× above `df_lambda`) and does not shorten the path to a
genuine λ advantage, which still needs ISDF/optimized collocation that preserves reconstruction
fidelity (out of scope, §7, unchanged).

**Hardware:** `Linux 7843f162c504 7.0.14-linuxkit x86_64`, 8 vCPU, container (no GPU). All of the
above (gate run + 3-part sweep) took ≈13 CPU-minutes wall time.

**Landing re-run (batch/ft-landing, 2026-10-02; macOS 27, Apple M3, 2 BLAS threads).** Gate:
4 passed in 30.4s. Sweep: `uv run python scripts/thc_collocation_sweep.py`, 542 s wall. The
deterministic numbers reproduce exactly: `df_lambda` 15.4846 / 86.9575, `random_lambda` 1001.6941
/ 5377.3005, and `n_params` 602. The LM fits reproduce qualitatively but not bit-for-bit:
- 12-seed reconstruction error is 3.65e-3 to 2.28e-2 (mean 1.40e-2), against 3.82e-3 to 2.02e-2
  above. Every seed is still above 1e-6.
- The 5-seed check still beats `df_lambda` at seeds 0, 2, 3, 4 but not seed 1 (λ = 16.6635), so
  4/5 again. `ratio_to_random` is in [0.0144, 0.0166].
- The norb=7 point has reconstruction error 9.00e-2, λ = 85.1301, ratio 0.0158, and still beats
  `df_lambda`.

Landing corrections to the arithmetic: the λ ratios 0.0144–0.0166 mean fit_thc sits ≈60–69× below
random collocation (this text said ≈65–69×). Random collocation sits 1001.69/15.48 = 64.7× (norb=6)
and 5377.30/86.96 = 61.8× (norb=7) above `df_lambda` (this text said ≈65–68×). Neither change
alters a verdict.

## 10. Deliverables

- `thc_factorization.py` — `zeta_from_collocation`, `pair_indicator_collocation`.
- `lambda_ladder.py` — `fit_thc(..., max_nfev=, return_factors=)`.
- `tests/test_thc_collocation_spec.py` — gates G1–G4.
- `scripts/thc_collocation_sweep.py` — the §9 corroboration sweep (12-seed reconstruction-error
  sweep, 5-seed beats-df spot check, norb=7 report point).
- Results summary (λ values, ratios, reconstruction errors, the diligent-sweep corroboration) in
  §9 above and in the PR description / handoff.
