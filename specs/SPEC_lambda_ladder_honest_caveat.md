# SPEC: lambda_ladder — the docstring's honest caveat, and an orbital-gauge artifact found while gating it

**Status:** ACCEPTED. G1/G2/G4 hold as originally written; G3 was retracted and revised (chem-1yr,
2026-09-26) after its "non-monotonic" claim was found to be an unpinned-orbital-gauge artifact, not
a property of double-factorization rank truncation. G5 was added to keep that retraction's evidence
runnable. See Sec 8 R2.

---

## 1. Goal

`lambda_ladder.py` compares three Hamiltonian factorizations (naive Pauli LCU, double
factorization, tensor hypercontraction) on the qubitization 1-norm lambda — the FT-QPE cost driver
— and its docstring already states an honest caveat: "on a small CAS, a generic THC fit is
comparable to or denser than DF; do not expect a small-system win." That claim has never been gated,
only demonstrated by eye in a `__main__` printout. This spec pins it. An earlier revision of this
spec also claimed DF's own rank-truncation accuracy is not monotonic in rank (a HIGHER rank
sometimes gives a WORSE energy than a lower one) — investigated as bead chem-1yr and RETRACTED: the
non-monotonicity was a symptom of an unpinned orbital gauge in N2's degenerate HOMO pi pair, not a
property of DF itself. See Sec 8 R2 for the finding and G3/G5 for the revised, evidenced claim.
False if THC ever beats DF at the cheapest tested rank on this system, or if DF/THC fail to
reconstruct exactly at full rank.

## 2. Background and honest framing

- `lambda_ladder.py` already reuses validated primitives (`qubitization_blueprint`'s exact Pauli
  decomposition, `df_factorization`'s double factorization) — no new physics, only falsifiers around
  claims the module's own docstring and `__main__` output already make but never check.
- **What you can claim if the gates pass:** the docstring's honest caveat is TRUE, not just prose —
  on N2 CAS(3 orbitals, 4 electrons) (the module's own `__main__` system, gauge-pinned with
  `symmetry='D2h'` — see Sec 8 R2), at the cheapest tested rank, THC's 1-norm is measurably LARGER
  than DF's (denser, not a win) and its LCU term count is dramatically larger (THC needs ~5x more
  terms than DF at comparable rank here); both methods reconstruct the Hamiltonian exactly at their
  own full rank (lambda matches the naive/exact value, FCI error is zero); once the orbital gauge is
  pinned, DF's rank-truncation accuracy is in fact monotonically non-increasing in rank on this
  system (G3, revised).
- **What you cannot claim:** that THC never beats DF (only checked at the cheapest rank on one
  system, matching the docstring's own stated scope — the asymptotic advantage needs "dozens of
  orbitals," out of scope here); that DF's rank-truncation accuracy is monotonic in general — only
  checked on this one gauge-pinned system, and G5 demonstrates that the *unpinned* version of this
  same system reaches non-monotonic error sequences for roughly half of a swept family of equally
  valid, equal-energy orbital gauges, so "monotonic" here is a property of this specific canonical
  gauge choice, not a proof that DF rank truncation is monotonic in general.
- **Reference:** the exact CASCI energy (`fci.direct_spin1.kernel` via `hybrid_quantum_solver.dmrg_reference.fci_energy`,
  through the module's own `fci_energy_error` helper) for the accuracy claims; the naive brute-force
  Pauli 1-norm (`lambda_and_terms` at the exact ERI, the module's own function) for the full-rank
  reconstruction claims.

## 3. Approach

Reuse `lambda_and_terms`, `fit_thc`, `fci_energy_error` (from `lambda_ladder.py`) and
`double_factorize`, `reconstruct_eri` (from `df_factorization.py`) unmodified — the same building
blocks `lambda_ladder()`'s print loop already calls, but capturing return values directly instead of
parsing printed output. System: N2 CAS(3 orbitals, 4 electrons), the module's own `__main__`
example, DF ranks `1..full_rank`, THC ranks `2..6` (the module's own default `thc_ranks=range(2,7)`).
The test fixture builds the molecule with `symmetry='D2h'` (not the module's own unadorned
`gto.M(...)`) to pin the RHF orbital gauge in N2's exactly-degenerate HOMO pi / LUMO pi* pairs — see
Sec 8 R2 for why this is necessary for G3's reproducibility and G2's exact "Measured" digits.

## 4. Public interface

No new library code — this spec adds only test-file assertions around `lambda_ladder.py`'s and
`df_factorization.py`'s existing functions, reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — full-rank reconstruction is exact for both methods.** DF at its own reported `full_rank`
  and THC at `M=6` (matching `full_rank` on this system) both reproduce the naive lambda to
  `< 1e-6` and give zero FCI error (`< 1e-6` mHa). *Measured (symmetry='D2h' gauge): naive
  lambda=8.4228763, terms=34; DF R=6 (full_rank=6) lambda=8.4228763, err 0.0000 mHa; THC M=6
  lambda=8.4228763, err 0.0000 mHa.*
- **G2 — THE FINDING: the honest caveat is true, not just prose.** At the cheapest tested rank
  (DF R=1 vs THC M=2, the module's own lowest `thc_ranks` value): THC's lambda exceeds DF's by more
  than 20%, and THC's term count exceeds DF's by at least 3x — "comparable to or denser than DF" at
  small CAS is a checked inequality here, not an assertion. *Measured: DF R=1 lambda=7.871482,
  terms=22; THC M=2 lambda=10.696709 (+35.9%), terms=109 (4.95x).*
- **G3 — REVISED (chem-1yr): DF's rank-truncation accuracy is monotonic once the orbital gauge is
  pinned.** Originally this gate asserted the opposite (a higher rank gives a WORSE FCI error than a
  lower rank somewhere in the sweep) and was measured true under an *unpinned* gauge. Investigation
  (bead chem-1yr, 2026-09-26) found that result was an artifact of RHF's arbitrary orthonormal choice
  within N2's exactly-degenerate HOMO pi pair, which CASCI(3,4)'s core/active freeze turns into a
  physically consequential (not null) decision — see Sec 8 R2. Under the canonical `symmetry='D2h'`
  gauge the error sequence is monotonically non-increasing every time (10/10 fresh-process runs).
  *Measured: err_mHa for R=1..6 = [141.673, 104.597, 52.431, 47.649, 5.716, 0.000] — strictly
  decreasing.* G5 keeps the retracted claim's gauge-freedom evidence runnable rather than just prose.
- **G4 — sanity: THC's own accuracy-vs-rank is not similarly broken at the endpoints tested.** THC's
  FCI error at `M=6` (full rank) is exact, and its error at the lowest tested `M=2` is the WORST of
  the swept range (confirms THC's fit isn't accidentally non-monotonic in a way that would confound
  G2's "cheapest rank" comparison). *Measured: THC err_mHa for M=2..6 = [2606.032, 69.313, 47.649,
  5.716, 0.000] — M=2 is the worst, consistent with G2's premise that "cheapest" and "least
  accurate" align for THC here.*
- **G5 — NEW (chem-1yr): the retracted G3 claim was an unpinned-orbital-gauge artifact, evidenced not
  asserted.** Build the *unpinned* system (plain `gto.M`, no symmetry), confirm mo indices 4 and 5
  are exactly degenerate, then sweep an explicit 2x2 rotation of that pair through 13 angles in
  `[0, pi]` (a full period of the gauge freedom) and rerun CASCI + the DF rank sweep at each angle.
  Assert the swept family contains BOTH a monotonic and a non-monotonic outcome — proof that the
  same physical system (same RHF energy at every angle) can produce either verdict depending purely
  on which equally-valid gauge is selected. *Measured: across 10 repeated fresh-process runs, the
  non-monotonic count varied 5-7 out of 13 angles, but every run produced a mix of both outcomes (the
  test only asserts the mix, not the exact count, since the count depends on the unpinned baseline's
  arbitrary starting phase).*

> Definition of done: **G2**. G1 builds the endpoint sanity; G3 (revised) and G4 keep the
> surrounding trend honest; G5 is the falsifiable record of why G3 needed revising.

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
- Explaining WHY the specific unpinned-gauge rotation angles that give a non-monotonic DF error
  sequence do so (i.e. the mechanism within double factorization) — G5 establishes THAT gauge
  determines the outcome, not WHY particular rotations break monotonicity.
- Systems beyond N2 CAS(3,4) — a natural follow-up, not attempted here.
- Fixing the underlying gauge-dependence of DF-truncated accuracy in general (e.g. by canonicalizing
  degenerate subspaces inside `df_factorization.py` itself) — out of scope for this spec, which only
  needed a reproducible test fixture; a real fix (if one is wanted) belongs in a new spec.

## 8. Caveats and risks

- **R1 — one system.** The specific magnitudes (35.9% lambda gap, 4.95x term-count gap, the R=1..6
  error sequence) are measurements on N2 CAS(3,4); the falsifiable claims (G2's direction, G3's
  monotonicity) are not claimed to generalize in magnitude, or in the monotonic/non-monotonic
  verdict itself, to other systems or other active-space choices.
- **R2 — orbital gauge changes DF-truncated accuracy when CASCI freezes part of a degenerate
  manifold (chem-1yr, 2026-09-26).** N2/sto-3g's HOMO pi pair (mo indices 4,5) and LUMO pi* pair
  (7,8) are exactly degenerate under D∞h. RHF's SCF solver returns SOME orthonormal basis for each
  degenerate pair, and without a symmetry constraint that specific basis is arbitrary and varies
  run to run (BLAS/LAPACK `eigh` tie-breaking, thread-count dependent). This is normally a null
  gauge freedom — any rotation gives the same total energy — but CASCI(norb=3, ne=4)'s default
  energy-ordered active-space selection freezes only ONE member of the HOMO pair (orbital 4) into
  "core" while keeping the other (orbital 5) active. Freezing breaks the pair's rotational
  invariance, so which specific rotation RHF happened to return becomes a physically consequential
  choice: it changes the frozen-core Fock contribution, hence `get_h1eff()`/`get_h2eff()`, hence
  CASCI's total energy and every downstream DF/THC number. This was originally discovered as gate
  flakiness (`test_G3` passing or failing depending on process-level BLAS nondeterminism, tracked as
  bead chem-1yr) and confirmed as a genuine gauge effect — not incidental noise — by explicitly
  sweeping a 2x2 rotation of the degenerate pair through a full period and observing the DF error
  sequence's monotonicity flip along with CASCI's total energy itself (G5). Two fixes are available:
  (a) pin the gauge with `symmetry='D2h'` so PySCF solves each real abelian irrep block
  independently and reproducibly (used here, for the G1-G4 fixture), or (b) accept the
  gauge-dependence as a real property of this construction and test for its EXISTENCE rather than a
  specific direction (used for G5). Loosening G3's threshold after seeing this data would have
  violated the repo's "never loosen thresholds after seeing data" rule; retracting and replacing the
  claim, with the replacement backed by a new falsifier (G5), is the sanctioned alternative
  (specs/README.md: "if a gate proves unsatisfiable, revise the spec and record why").
- Honest limitation: brute-force Pauli 1-norm only (feasible at this scale); `fit_thc`'s stochastic
  least-squares fit is seeded (`seed=0` default) for reproducibility, not claimed globally optimal.

## 9. Deliverables

- `tests/test_lambda_ladder_honest_caveat_spec.py` — gates G1-G4 (unchanged in substance; fixture
  gained `symmetry='D2h'`) plus revised G3 and new G5 (chem-1yr). No library code changes.
- Results summary (with the R1/R2 caveats) in the PR description / BACKLOG entry.
