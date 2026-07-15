# SPEC: lambda_ladder — the docstring's honest caveat, and an unmonotonic accuracy trend it doesn't mention

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`lambda_ladder.py` compares three Hamiltonian factorizations (naive Pauli LCU, double
factorization, tensor hypercontraction) on the qubitization 1-norm lambda — the FT-QPE cost driver
— and its docstring already states an honest caveat: "on a small CAS, a generic THC fit is
comparable to or denser than DF; do not expect a small-system win." That claim has never been gated,
only demonstrated by eye in a `__main__` printout. This spec pins it, and — probing this spec — finds
something the docstring does NOT mention: DF's own rank-truncation accuracy is not monotonic in rank
(a HIGHER rank sometimes gives a WORSE energy than a lower one). False if THC ever beats DF at the
cheapest tested rank on this system, or if DF/THC fail to reconstruct exactly at full rank.

## 2. Background and honest framing

- `lambda_ladder.py` already reuses validated primitives (`qubitization_blueprint`'s exact Pauli
  decomposition, `df_factorization`'s double factorization) — no new physics, only falsifiers around
  claims the module's own docstring and `__main__` output already make but never check.
- **What you can claim if the gates pass:** the docstring's honest caveat is TRUE, not just prose —
  on N2 CAS(3 orbitals, 4 electrons) (the module's own `__main__` system), at the cheapest tested
  rank, THC's 1-norm is measurably LARGER than DF's (denser, not a win) and its LCU term count is
  dramatically larger (THC needs ~5x more terms than DF at comparable rank here); both methods
  reconstruct the Hamiltonian exactly at their own full rank (lambda matches the naive/exact value,
  FCI error is zero).
- **What you cannot claim:** that THC never beats DF (only checked at the cheapest rank on one
  system, matching the docstring's own stated scope — the asymptotic advantage needs "dozens of
  orbitals," out of scope here); that DF's non-monotonic accuracy (G3, the extra finding) is a bug —
  it may be an intrinsic property of rank-truncated double factorization (dropping a mid-magnitude
  eigenvalue can remove a partially-cancelling error term, making an intermediate rank accidentally
  worse than both its neighbors) — recorded as a measured, unexplained phenomenon, not diagnosed.
- **Reference:** the exact CASCI energy (`fci.direct_spin1.kernel`, the module's own `fci_energy_error`
  helper) for the accuracy claims; the naive brute-force Pauli 1-norm (`lambda_and_terms` at the
  exact ERI, the module's own function) for the full-rank reconstruction claims.

## 3. Approach

Reuse `lambda_and_terms`, `fit_thc`, `fci_energy_error` (from `lambda_ladder.py`) and
`double_factorize`, `reconstruct_eri` (from `df_factorization.py`) unmodified — the same building
blocks `lambda_ladder()`'s print loop already calls, but capturing return values directly instead of
parsing printed output. System: N2 CAS(3 orbitals, 4 electrons), the module's own `__main__`
example, DF ranks `1..full_rank`, THC ranks `2..6` (the module's own default `thc_ranks=range(2,7)`).

## 4. Public interface

No new library code — this spec adds only test-file assertions around `lambda_ladder.py`'s and
`df_factorization.py`'s existing functions, reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — full-rank reconstruction is exact for both methods.** DF at its own reported `full_rank`
  and THC at `M=6` (matching `full_rank` on this system) both reproduce the naive lambda to
  `< 1e-6` and give zero FCI error (`< 1e-6` mHa). *Measured: DF R=6 lambda=8.345968 (naive
  8.345968), err 0.0000 mHa; THC M=6 identical to 6 decimal places.*
- **G2 — THE FINDING: the honest caveat is true, not just prose.** At the cheapest tested rank
  (DF R=1 vs THC M=2, the module's own lowest `thc_ranks` value): THC's lambda exceeds DF's by more
  than 20%, and THC's term count exceeds DF's by at least 3x — "comparable to or denser than DF" at
  small CAS is a checked inequality here, not an assertion. *Measured: DF R=1 lambda=8.065, terms=22;
  THC M=2 lambda=10.614 (+32%), terms=117 (5.3x).*
- **G3 — boundary, recorded not smoothed over: DF's rank-truncation accuracy is NOT monotonic.**
  Somewhere in the DF rank sweep (R=1..full_rank), a higher rank gives a WORSE FCI error than a
  lower rank — the docstring never mentions this, and a reader could reasonably assume "more
  factors = strictly better" (as G1's endpoint behavior would suggest). *Measured: R=4 err 12.15
  mHa, R=5 err 43.53 mHa — R=5 is worse than R=4, before R=6 recovers to exact.*
- **G4 — sanity: THC's own accuracy-vs-rank is not similarly broken at the endpoints tested.** THC's
  FCI error at `M=6` (full rank) is exact, and its error at the lowest tested `M=2` is the WORST of
  the swept range (confirms THC's fit isn't accidentally non-monotonic in a way that would confound
  G2's "cheapest rank" comparison). *Measured: THC errors 2622.6, 49.3, 12.15, 43.5, 0.0 mHa for
  M=2..6 — M=2 is the worst, consistent with G2's premise that "cheapest" and "least accurate" align
  for THC here.*

> Definition of done: **G2**. G1 builds the endpoint sanity; G3/G4 keep the surrounding trend honest
> rather than implying a smooth story the data doesn't support.

## 6. Implementation plan (test-first)

1. Write `tests/test_lambda_ladder_honest_caveat_spec.py` encoding G1-G4 (RED in the sense these
   checks are new, even though `lambda_ladder.py`'s functions are not).
2. No changes to `lambda_ladder.py` or `df_factorization.py` — every gate calls their existing
   public functions directly.
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- THC's asymptotic advantage at large active spaces (the docstring's own stated regime, "dozens of
  orbitals") — not reachable with the brute-force Pauli 1-norm this module uses (feasible only for
  `norb <= ~4` per its own docstring).
- Diagnosing WHY DF's rank-5 accuracy dips below rank-4's (G3) — recorded as a measured phenomenon,
  not investigated further.
- Systems beyond N2 CAS(3,4) — a natural follow-up, not attempted here.

## 8. Caveats and risks

- **R1 — one system.** The specific magnitudes (32% lambda gap, 5.3x term-count gap, the R=5 dip)
  are measurements on N2 CAS(3,4); the falsifiable claims (G2's direction, G3's existence of
  non-monotonicity) are not claimed to generalize in magnitude to other systems.
- Honest limitation: brute-force Pauli 1-norm only (feasible at this scale); `fit_thc`'s stochastic
  least-squares fit is seeded (`seed=0` default) for reproducibility, not claimed globally optimal.

## 9. Deliverables

- `tests/test_lambda_ladder_honest_caveat_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
