# SPEC: Certified two-sided brackets on a spectral gap (no FCI oracle)

**Status:** IMPLEMENTED — gates G1–G4 green. The definition-of-done gate (G1) is zero escapes at M≥6.

---

## 1. Goal

The fundamental gap Δ = E₁ − E₀ of the HF-reachable sector can be **bracketed**, from Krylov data
alone (no FCI), by an interval [Δ_lo, Δ_hi] that provably contains the exact reachable gap once the
Krylov space resolves the excited state (M ≥ 6). Falsifiable: a single depth at which the exact
reachable gap falls outside the certified interval (in the M ≥ 6 regime) kills the claim.

## 2. Background and honest framing

`temple_bracket` gave every QKSD ground-state *energy* a certified two-sided bracket. Spectroscopy
measures *gaps*, not total energies, and those were still point estimates. This spec extends the
Temple/Weinstein machinery to the gap, at the cost of one extra ⟨H²⟩ on the first-excited Ritz state
(no new measurements beyond what the excited solve already needs).

- **What we can claim:** a certified error bar on an excitation gap computed without the exact
  answer, valid under the same checkable premise as the ground-state bracket, demonstrated with zero
  escapes on H₄ / LiH / N₂ CAS(6,6) and closing to sub-mHa (H₄) / tens of mHa (N₂).
- **What we cannot claim:** unconditional two-sided rigor. The lower certificate rests on the
  Weinstein premise ε₁ ≤ E₁ (θ₁ resolves E₁), which a real-time Krylov space cannot self-verify and
  which fails at M ≤ 4 for multireference cases — the recorded boundary. Sector-restricted (two
  lowest HF-reachable levels). Exact statevector: the shot cost of ⟨H²⟩ (~λ²-sized Pauli expansion)
  is not modeled. On systems small enough to diagonalize (the 4-qubit Nb₃X₈ dimer) certification is
  pointless — the value is where FCI is out of reach.

## 3. Approach

**Reference = the exact reachable gap** (dense diagonalization + HF overlap, `reachable_gap`, used
only to *check* the bracket in tests — never fed to the estimator). The bracket composes:

- Cauchy interlacing: θ₁ ≥ E₁ (second Ritz value bounds E₁ from above).
- Temple (1928): τ₀ ≤ E₀ (ground lower bound, using ε₁ as the E₁ floor).
- Weinstein / self-ε: ε₁ = θ₁ − σ₁ ≤ E₁ under the premise that θ₁ resolves E₁.

With θ₀ ≥ E₀ (variational):

```
Δ_hi = θ₁ − τ₀                >= E₁ − E₀ = Δ      (upper certificate: interlacing + Temple)
Δ_lo = (θ₁ − σ₁) − θ₀         <= Δ                (lower certificate: premise-gated)
```

**Numeric result (self-ε mode, no oracle):** exact reachable gap inside [Δ_lo, Δ_hi] at every M ≥ 6
for all three systems; width closes H₄ 342 → 0.7 mHa, LiH 181 → 3.3 mHa, N₂ 391 → 37 mHa over
M = 6…24. At M = 4 the premise fails for H₄/N₂ and the lower certificate escapes.

## 4. Public interface

Reuses `QuantumKrylovSolver.eigenstates`, `temple_bounds._mean_and_variance`.

```
certified_gaps.gap_bracket(mh, m, e1=None, solver=None) -> GapBracket
certified_gaps.gap_bracket_ladder(mh, dims, e1=None, solver=None) -> list[GapBracket]
certified_gaps.reachable_gap(mh) -> float          # dense, reference/validation only
GapBracket(m, gap_lower, gap_upper, width, theta0, theta1, sigma1, eps1, eps1_source)
certified_gaps (CLI)                               # per-system M-ladder + inside? column
```

## 5. Acceptance criteria (validation gates)

`tests/test_certified_gaps_spec.py` (test-first).

- **G1 — zero escapes (DEFINITION OF DONE).** For every system and every M ≥ 6, the exact reachable
  gap lies in [Δ_lo, Δ_hi] (self-ε mode). One escape kills the claim.
- **G2 — the bracket closes.** width(M=24) < ½·width(M=6) and < 50 mHa for every system — a genuine
  error bar, not a vacuous one; finite throughout the certified regime.
- **G3 — the boundary (finding).** At M = 4 the premise ε₁ ≤ E₁ fails for the multireference cases
  (H₄, N₂; checked against the oracle E₁ the live path never has) and the lower certificate escapes
  — certification requires M ≥ 6. Failure need only *exist*, not be universal (LiH may pass at M=4).
- **G4 — upper certificate is robust + scope.** The upper certificate brackets the gap from above at
  every tested depth including M = 4 (the premise-sensitive side is the lower one); the reachable
  gap is positive; oracle mode (feeding exact E₁) never escapes in the certified regime.

## 6. Implementation plan (test-first)

1. `tests/test_certified_gaps_spec.py` encoding G1–G4 (initially failing — no module).
2. `certified_gaps.py` composing interlacing + Temple + Weinstein from the validated primitives.
3. `make gates` (own process; no block2).

## 7. Out of scope

- Excited states beyond E₁ (interior levels need lower bounds on E₂, E₃… — a follow-up).
- Self-verification of the Weinstein premise (fundamentally needs a lower bound on E₂).
- Hardware ⟨H²⟩ shot-cost modeling.
- Certification on exactly-diagonalizable toys (Nb₃X₈ dimer) — pointless by construction.

## 8. Caveats and risks

- **R1 — premise unverifiable without an oracle.** Mitigated by inheriting the temple_bracket M ≥ 6
  boundary and gating the M = 4 failure explicitly (G3), so the interval is never quoted below the
  depth where it is trustworthy.
- The lower certificate can be negative (vacuous but valid) at shallow depth; quote a gap ± only
  once the width is finite and small and M ≥ 6.

## 9. Deliverables

- `certified_gaps.py` — gap-bracket estimator + CLI.
- `tests/test_certified_gaps_spec.py` — gates G1–G4.
- `specs/SPEC_certified_gaps.md` — this spec.
- `specs/BACKLOG.md` — entry with the finding recorded.
