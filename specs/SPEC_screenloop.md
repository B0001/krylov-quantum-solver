# SPEC: ScreenLoop — bracket-aware screening with zero false eliminations (#19)

**Slug:** `screenloop` · **Depends on:** CertChem-M1 certified brackets (via a `BoundedOracle`
adapter; the core needs only intervals). Tasks `specs/tasks/14-screening-loop.md`.

## Goal

Screen a candidate space for those whose (unknown) property lands in a target region, eliminating
candidates by **certified interval dominance** with a *provable* guarantee of **zero false
eliminations** — something a point-estimate screener cannot promise. Demonstrated on a
50-candidate toy space with exhaustive ground truth.

## The rule (sound by construction)

A candidate carries a certified bracket `[lo, hi]` that (ADR-0001 containment) always encloses its
true value. Target region `[t_lo, t_hi]`. **Eliminate iff the bracket strictly excludes the
target:** `hi < t_lo` or `lo > t_hi`. Because the true value lies in `[lo, hi]`, an eliminated
candidate's true value cannot lie in `[t_lo, t_hi]` — so a genuine hit is *never* eliminated
(zero false eliminations). Boundary contact is resolved conservatively (non-strict overlap keeps
the candidate). A survivor whose bracket is fully inside the target (`t_lo <= lo, hi <= t_hi`) is a
**confirmed** hit; otherwise it is **undecided** and gets refined.

## Interface

```python
from screenloop import (
    Interval, Verdict, classify,             # the dominance rule (pure)
    BoundedOracle, SyntheticOracle,          # interval source protocol + a test oracle
    screen,                                   # the loop (v1 prune-only, v2 acquisition)
    point_estimate_screen,                    # the unsound baseline (for contrast)
)
```

## Acceptance gates (`tests/test_screenloop_spec.py`)

1. **Dominance rule.** `classify` eliminates iff strictly disjoint; boundary contact keeps;
   full-containment → confirmed. Unit-tested on hand cases.
2. **Zero false eliminations (the invariant).** On a 50-candidate synthetic space with known
   ground truth (~10 true hits), `screen` eliminates **no** true hit — at pilot precision and
   through refinement. This is the claim.
3. **Acquisition is cheaper.** v2 (overlap-ordered refinement) uses fewer oracle evaluations than
   exhaustive-at-full-precision, with the identical final hit set.
4. **Baseline is unsound.** `point_estimate_screen` (eliminate by best-estimate) makes **> 0**
   false eliminations on the same space where `screen` makes 0 — the demonstration.
5. **Oracle is generic.** A second `BoundedOracle` (a synthetic conformal-ML mock with different
   interval widths) runs through `screen` unchanged, still zero false eliminations.

## Out of scope / caveats

- The guarantee is conditional on the oracle's brackets actually containing truth (CertChem's
  contract for the certified oracle; assumed for mock oracles). Garbage brackets → garbage
  guarantee; the soundness is *relative to* bracket validity.
- A certchem-backed oracle adapter is provided and smoke-tested, but the headline gates run on
  fast synthetic oracles so the invariant is checked over a large space deterministically.
