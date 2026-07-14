# SPEC: Gap self-check under shot noise — intersection concentrates noise where composition diluted it

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`gap_selfcheck.py`'s oracle-free trustworthiness certificate (a bracket is "corroborated" iff it
overlaps the intersection of the deepest brackets; the self-checked interval is the intersection of
the corroborated ones) was validated only exact-statevector. This spec asks the noise question
`certified_thermochem_noise.md` already asked of a *difference* composition, but here of an
*intersection* composition: does the self-checked interval retain useful coverage of the true gap
under i.i.d. shot noise, and does the existing z=2 single-bracket rule (`certified_noise`) restore
it? The claim is false if intersecting several independently-noisy brackets turns out to need
*less* inflation than a single bracket (the same direction `certified_thermochem_noise` found) —
the interesting content here is the contrast.

## 2. Background and honest framing

- Builds on `certified_noise.md` (single-bracket coverage collapse + z-inflation fix),
  `certified_gaps.md`/`gap_selfcheck.md` (the exact-statevector self-check), and is the intersection
  counterpart to `certified_thermochem_noise.md`'s difference composition (same repo, two composition
  operators, opposite findings).
- **What you can claim if the gates pass:** composing certified brackets by AND/intersection
  (self-check's mechanism) concentrates noise-induced miscalibration rather than diluting it —
  the opposite of composing by difference (thermochem). This is a genuine, actionable caveat for
  anyone tempted to run `gap_selfcheck` on noisy/hardware data expecting the same z=2 rule to work.
- **What you cannot claim:** a closed-form derivation of the required z (empirical, per-system, like
  `certified_thermochem_noise`'s z*); a fix for the miscalibration (out of scope — this spec measures
  the problem, a follow-up would design the repair); non-i.i.d. / grouped-Pauli noise; more than the
  two systems tested.
- **Reference:** the exact reachable gap (`certified_gaps.reachable_gap`, dense diagonalization) —
  the same reference `gap_selfcheck`'s own exact-statevector gates already use.

## 3. Approach

Reuse `gap_selfcheck`'s corroboration/intersection logic (`self_checked_gap`, unmodified) and
`certified_noise`'s i.i.d. shot-noise model and inflation convention (`certified_half_width =
z*lambda_H/sqrt(shots)`), composed the simplest, least ambiguous way: build ONE noisy realization of
each depth's bracket (single noisy `theta0, var0, theta1, sigma1` per depth, matching
`certified_gaps.gap_bracket`'s exact formula — no internal directional surgery, which the probe for
this spec found is genuinely ambiguous for a two-eigenstate gap bracket, unlike the ground-state
Temple case), then pad each depth's `[gap_lower, gap_upper]` post-hoc by
`+/- certified_half_width(lambda_H, shots, z)` before handing the padded ladder to
`gap_selfcheck.self_checked_gap` unchanged. Monte Carlo (>= 3000 trials, seeded) over noise
realizations measures coverage of the exact gap on H4 and LiH (the same two systems
`gap_selfcheck`'s own exact-statevector example already uses, at `dims=(6,8,12,16,20,24)`).

## 4. Public interface

```
gap_selfcheck_noise.self_check_noise_coverage(
    mh, dims, shots, z=2.0, trials=4000, seed=0, k=2,
) -> dict   # coverage, frac_empty (no bracket corroborated), lam_h, lam_h2
gap_selfcheck_noise.minimal_z_for_selfcheck_coverage(
    mh, dims, shots, target=0.9, resolution=0.25, z_max=6.0, ...
) -> float
```

## 5. Acceptance criteria (validation gates)

- **G1 — the self-checked interval inherits the collapse, and it is shot-count-independent.** Raw
  (z=0) coverage on H4 and LiH is < 0.3 at shots in {1e4, 1e5, 1e6} (dims as above), and the spread
  across those three shot counts is < 0.1 per system. *Measured: H4 0.145/0.148/0.198 (spread
  0.053); LiH 0.081/0.056/0.077 (spread 0.025).*
- **G2 — THE CONTRAST: the existing z=2 rule does NOT fix it here.** At shots=1e5, z=2.0 inflation
  (the rule that restores `certified_noise`'s single bracket to ~0.98 and is already MORE than
  enough for `certified_thermochem_noise`'s composed bracket) leaves coverage < 0.85 on BOTH H4 and
  LiH. *Measured: H4 0.702, LiH 0.527.*
- **G3 — THE FINDING (definition of done): the self-check needs MORE inflation than a single
  bracket, not less.** The minimal z (grid search, 0.25 resolution) reaching 90% coverage at
  shots=1e5 is > 2.0 on both systems — the opposite direction from `certified_thermochem_noise`'s
  composed bracket (which needed *less* than z=2). *Measured: z\*_H4 = 3.25, z\*_LiH = 4.00 —
  roughly 1.6-2x the single-bracket rule, not a fraction of it.*
- **G4 — honesty diagnostic.** The self-checked interval is empty (no bracket corroborates —
  `gap_selfcheck._intersect` on zero kept brackets) in a non-trivial fraction of raw (z=0) trials
  (> 0.05 on at least one system) but that fraction drops below 0.02 once z=2.0 padding is applied —
  padding fixes "inconclusive" well before it fixes "correct." *Measured: frac_empty z=0: H4 0.11-0.13,
  LiH 0.09-0.13; z=2.0: < 0.005 both.*

> Definition of done: **G3**. If a future system is found where the self-check needs *less*
> inflation than a single bracket (contradicting the direction found here), that is the boundary and
> must be recorded, not silently averaged away.

## 6. Implementation plan (test-first)

1. Write `tests/test_gap_selfcheck_noise_spec.py` encoding G1-G4 (RED — module doesn't exist).
2. Add `gap_selfcheck_noise.py`: `self_check_noise_coverage` (per-depth single noisy realization +
   post-hoc `certified_half_width` padding, reusing `gap_selfcheck.self_checked_gap` unmodified) and
   `minimal_z_for_selfcheck_coverage` (grid search over z).
3. `make gates` / targeted pytest to green; ruff clean.

## 7. Out of scope

- A repair for the miscalibration (a follow-up spec — e.g. a Bonferroni-style per-depth confidence
  correction was probed informally and did not fully fix it either; that is itself worth a follow-up,
  not asserted here).
- A closed-form reason *why* intersection needs more z while difference needs less (empirical only).
- Oracle e1 mode, more than two composed geometries/systems, non-Gaussian or grouped-Pauli noise.

## 8. Caveats and risks

- **R1 — the specific z\* numbers are system-dependent (like `certified_thermochem_noise`'s R1).**
  H4 and LiH give different z* (3.25 vs 4.00); a system with a wider corroboration margin at depth
  might need less. The falsifiable claim (G3) is the *direction* (self-check needs more, not less),
  kept system-relative rather than a universal constant.
- **R2 — coverage at fixed z degrades further as shots grow** (the absolute pad
  `z*lambda_H/sqrt(shots)` shrinks with more shots, same as `certified_noise`'s half-width, so a
  fixed z inevitably drifts toward the (broken) raw coverage as shots -> infinity). This is expected,
  not a new pathology, but means z must be re-derived per shot budget, not fixed once.
- Honest limitation: two systems, one padding convention (simplest defensible choice — post-hoc
  half-width, not an internally-optimized directional inflation, which the probe found is genuinely
  ambiguous for a two-eigenstate gap bracket).

## 9. Deliverables

- `gap_selfcheck_noise.py` — `self_check_noise_coverage`, `minimal_z_for_selfcheck_coverage`.
- `tests/test_gap_selfcheck_noise_spec.py` — gates G1-G4.
- Results summary (with R1/R2 caveats) in the PR description / BACKLOG entry.
