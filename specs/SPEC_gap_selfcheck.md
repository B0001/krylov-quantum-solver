# SPEC: Oracle-free trustworthiness certificate for the certified gap bracket

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G2 (oracle-free interval
covers the truth).

---

## 1. Goal

The certified gap bracket (`certified_gaps`) can **self-verify without an oracle**: a bracket at
Krylov depth M is *corroborated* iff it is consistent (overlaps) the converged brackets at greater
depth, and this adaptively catches the premise-failure regime — rejecting the shallow bracket only
when it is genuinely inconsistent. Falsifiable: if a premise-failed bracket (M=4 on H₄/N₂) were
*not* rejected, or a valid one (M=4 on LiH, M≥6 anywhere) *were* rejected, the certificate is wrong.

## 2. Background and honest framing

`certified_gaps` recorded an open limitation: its lower certificate rests on the premise ε₁ ≤ E₁,
which fails at shallow depth and is verifiable only against an oracle (exact E₁). This spec removes
the oracle by cross-depth consistency: a premise failure inflates/shifts a bracket so it no longer
overlaps the deep ones, and is caught with no FCI.

- **What we can claim:** an adaptive, oracle-free flag that distinguishes real premise failures from
  merely-shallow-but-valid brackets, and an oracle-free gap interval validated (against FCI, in
  tests only) to contain the exact reachable gap and to exclude the shallow outliers.
- **What we cannot claim:** sufficiency. The anchor — the deepest brackets — is taken on trust; the
  check cannot certify its own convergence, and cannot see a consistently-biased sequence shrinking
  toward a wrong value (the model-misspecification blind spot of any single-run resampling, cf.
  `SPEC_odmd_uq`). It catches the *known* failure mode and pairs with a depth-convergence check.

## 3. Approach

**Reference (tests only) = the exact reachable gap** (`certified_gaps.reachable_gap`), used to
confirm the oracle-free interval covers the truth. The certificate itself uses only the bracket
ladder:

```
anchor          = intersection of the deepest k brackets           (convergence reference)
corroborated(M) = bracket(M) overlaps the anchor
self-checked    = intersection of the corroborated brackets        (oracle-free interval)
```

**Numeric result:** on H₄ / N₂ the M=4 bracket (premise failed, escapes) does not overlap the deep
anchor → rejected; on LiH the M=4 bracket overlaps → accepted (premise holds there); every M≥6 is
corroborated on all systems. The self-checked interval contains the exact gap for all three
(H₄ [624.5, 625.3] ⊇ 625.3; LiH [130.1, 133.4] ⊇ 133.4; N₂ [673.4, 709.9] ⊇ 708.0 mHa).

## 4. Public interface

Reuses `certified_gaps.gap_bracket_ladder`, `GapBracket`.

```
gap_selfcheck.anchor_interval(brackets, k=2) -> (lo, hi)
gap_selfcheck.corroborated_flags(brackets, k=2) -> list[bool]
gap_selfcheck.self_checked_gap(brackets, k=2) -> (lo, hi)
gap_selfcheck.self_checked_gap_from(mh, dims, k=2, solver=None) -> ((lo, hi), flags)
gap_selfcheck (CLI)                                             # per-system interval + rejected M
```

## 5. Acceptance criteria (validation gates)

`tests/test_gap_selfcheck_spec.py` (test-first).

- **G1 — adaptive rejection.** Rejects M=4 for H₄, N₂ (premise failures) but accepts M=4 for LiH
  (premise holds); accepts every M≥6 on all systems. Oracle-free.
- **G2 — oracle-free interval covers the truth (DEFINITION OF DONE).** The self-checked interval is
  non-empty, finite, and contains the exact reachable gap for every system.
- **G3 — the self-check repairs the naive estimate.** Naively intersecting *all* brackets (incl. the
  M=4 outlier) is empty for the premise-failure cases; dropping the uncorroborated M=4 restores a
  non-empty interval containing the gap.
- **G4 — anchor robustness + minimal premise.** Flags identical for anchor depth k=2 vs k=3; the
  deep anchor is itself self-consistent (non-empty). The anchor is trusted, not proven (§2).

## 6. Implementation plan (test-first)

1. `tests/test_gap_selfcheck_spec.py` encoding G1–G4 (initially failing — no module).
2. `gap_selfcheck.py` composing the `certified_gaps` ladder with interval consistency.
3. `make gates` (own process; no block2).

## 7. Out of scope

- Certifying the convergence of the deepest bracket (fundamentally needs an external anchor).
- Detecting a consistently-biased (model-misspecified) sequence — the recorded blind spot.
- Shot-noise statistics on the flags (the brackets here are exact statevector).

## 8. Caveats and risks

- **R1 — necessary, not sufficient.** Stated up front; the certificate pairs with a
  depth-convergence check, never replaces it.
- The anchor depth k is a small heuristic; G4 gates flag-stability across k to bound the risk.

## 9. Deliverables

- `gap_selfcheck.py` — corroboration certificate + CLI.
- `tests/test_gap_selfcheck_spec.py` — gates G1–G4.
- `specs/SPEC_gap_selfcheck.md` — this spec.
- `specs/BACKLOG.md` — entry with the finding recorded.
