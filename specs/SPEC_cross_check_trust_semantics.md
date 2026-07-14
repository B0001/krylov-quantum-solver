# SPEC: cross_check's own trust semantics — what "reference" means when CASCI can't run

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`cross_check.py` is "the capstone for the near-term half of the stack" — four independent solvers
(CASCI, Krylov, ADAPT, SQD) agreeing is what justifies trusting a number before spending QPU time on
it. Its `__main__` informally asserts agreement, never CI-gated. But the harness's OWN correctness
has a subtler question its docstring doesn't address: `reference = available.get("CASCI",
next(iter(available.values())) ...)` — when CASCI is out of reach (a real, expected case; its own
docstring says CASCI runs "up to a large determinant dimension" but is capped), the "reference"
silently becomes whichever OTHER method happened to run first in Python dict insertion order
(Krylov, then ADAPT, then SQD) — not a principled "most trustworthy remaining method" choice. This
spec asks: is that fallback order actually sound, and does the harness know when it has degraded
from "checked against the exact answer" to "three approximate methods agreeing with a fourth
approximate method"? False if the fallback ever silently promotes the least reliable available
method (SQD) ahead of a more reliable one (Krylov/ADAPT), or if the harness cannot be made to report
an honest "no reference available" state when it should.

## 2. Background and honest framing

- `cross_check.py` already reuses every validated solver in this repo unmodified (`krylov_ground_state`,
  `adapt_ground_state`, `fci.direct_spin1`, `qiskit_addon_sqd`) — this spec adds no new physics, only
  falsifiers around the HARNESS's own trust logic, not any individual solver.
- **What you can claim if the gates pass:** with all four solvers reachable, they agree within the
  harness's own 5 mHa tolerance (pinning the informal `__main__` assertion); the reference fallback,
  while undocumented, is at least ordered sensibly — when CASCI is capped out, the reference
  degrades to Krylov (not silently to SQD, the more approximate configuration-sampling method); when
  CASCI and Krylov are BOTH capped out, it degrades to ADAPT before SQD — the priority order (CASCI
  > Krylov > ADAPT > SQD) an insertion-order accident happens to match a sensible general
  reliability ranking, at least in the direction of never promoting the noisiest method over a more
  precise one still running.
- **What you cannot claim:** that this priority order is PRINCIPLED (it is literally Python dict
  insertion order in the source, not a documented or configurable trust ranking) — a future edit
  that reorders the `res[...] = ...` assignments would silently change which method becomes the
  reference, with no test catching it except this one; that SQD is always the least accurate method
  in general (measured on H4 CAS(4,4) here, SQD's deviation from CASCI, 1.0e-9 mHa, was actually
  TIGHTER than Krylov's or ADAPT's — a system-specific observation, not a general ranking).
- **Reference:** the exact CASCI energy (dense diagonalization) as ground truth for G1; the harness's
  own returned `reference`/`results` fields, compared for exact equality against the individual
  solver values, for G2-G4 (a structural/mechanism check, not an accuracy check).

## 3. Approach

Reuse `cross_check` unmodified, calling it with different `fci_max_dim`/`krylov_max_dim`/
`qubit_dense_max_orb` cost caps (already public keyword arguments) to force specific solvers to be
skipped, and inspect the returned dict's `reference`/`results`/`skipped` fields — no new library
code, no monkeypatching, no internal access.

## 4. Public interface

No new library code — this spec adds only test-file assertions around `cross_check.py`'s existing
`cross_check` function, reused unchanged, driven entirely through its public cost-cap parameters.

## 5. Acceptance criteria (validation gates)

- **G1 — baseline agreement, pinning the informal `__main__` assertion.** On H2 CAS(2,2) and H4
  CAS(4,4) with all four solvers reachable: `out["agree"]` is `True` and `out["max_dev_mHa"] <=
  5.0`. *Measured H4: max deviation 0.0047 mHa (Krylov), well inside tolerance.*
- **G2 — THE FINDING: the reference fallback is exact insertion-order, not a documented ranking.**
  With `fci_max_dim=0` (CASCI forced unreachable) on H4 CAS(4,4): `out["reference"]` equals
  `out["results"]["Krylov"][0]` EXACTLY (not merely close) — proving the fallback is a literal
  first-available pick, not a recomputed "best" estimate.
- **G3 — the accidental priority order at least never promotes the noisiest method early.** With
  BOTH `fci_max_dim=0` and `krylov_max_dim=0` (CASCI and Krylov both forced unreachable): the
  reference becomes ADAPT (`out["reference"] == out["results"]["ADAPT"][0]` exactly), NOT SQD, even
  though SQD is also available — the insertion order (`CASCI, Krylov, ADAPT, SQD` in the source)
  happens to keep the configuration-sampling method last in line.
- **G4 — boundary, recorded not smoothed over: SQD has no cost-cap knob and cannot be suppressed
  through the public API.** With ALL THREE other caps forced to zero
  (`fci_max_dim=krylov_max_dim=qubit_dense_max_orb=0`): CASCI, Krylov, and ADAPT are all in
  `out["skipped"]`, but SQD still ran and became the sole reference
  (`out["reference"] == out["results"]["SQD"][0]`) — there is no way to drive `cross_check` into a
  fully-empty "no reference available" state through its public cost caps alone; SQD is
  unconditionally attempted, a real asymmetry in the four solvers' treatment worth knowing about
  before assuming the caps make all four symmetric knobs.

> Definition of done: **G2**. G1 alone would miss the whole point — the interesting question is not
> "do all four agree when everything runs" but "what does the harness actually DO when they can't
> all run," which is exactly the regime the docstring's own cost-cap design anticipates.

## 6. Implementation plan (test-first)

1. Write `tests/test_cross_check_trust_semantics_spec.py` encoding G1-G4 (RED in the sense these
   checks are new, even though `cross_check.py`'s function is not).
2. No changes to `cross_check.py` — every gate is driven through its existing public keyword
   arguments (`fci_max_dim`, `krylov_max_dim`, `qubit_dense_max_orb`).
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- Fixing the reference-fallback logic to be a documented, principled ranking instead of insertion
  order (a follow-up, if the maintainer wants one — this spec measures the current behavior).
- A general claim about which of Krylov/ADAPT/SQD is "most accurate" (system-specific, per §2).
- Any change to `cross_check.py`'s SQD-cannot-be-capped asymmetry (recorded as a boundary, not
  patched).

## 8. Caveats and risks

- **R1 — G2/G3's forced-unavailability scenarios are constructed via cost caps, not naturally
  occurring at this system size.** On H4 CAS(4,4) all four solvers are normally reachable (G1); the
  caps are a deliberate stress test of the fallback logic, matching how a LARGER system (where CASCI
  genuinely exceeds `fci_max_dim`) would actually behave — the mechanism, not the trigger condition,
  is what's being verified.
- Honest limitation: one system (H4 CAS(4,4)) for the fallback-order gates; the insertion-order
  finding is a fact about the CURRENT source text, not something that would be caught if the
  function were refactored without re-running this test (which is exactly why it's now gated).

## 9. Deliverables

- `tests/test_cross_check_trust_semantics_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
