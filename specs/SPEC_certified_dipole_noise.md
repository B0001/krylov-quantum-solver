# SPEC: Certified dipole under shot noise — inflation cannibalizes the gap margin it needs to stay finite

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`certified_dipole.py`'s property bracket rests on `sin(theta) <= sigma_0 / Delta_lo` (Davis-Kahan),
so its half-width is finite only while `Delta_lo` (the certified gap lower bound) stays positive and
exceeds `sigma_0`. `certified_noise.md`/`certified_thermochem_noise.md`/`gap_selfcheck_noise.md`
already showed z-inflation is the standard fix for shot-noise coverage collapse, in two different
compositions (difference: less z needed; intersection: more z needed) — both *monotonic* in z. This
spec asks the same question of the property bracket, which composes THREE noisy quantities
(`sigma_0`, `Delta_lo`, and the dipole's own `mu`/`sigma_A`) through the `s = sigma_0/Delta_lo < 1`
gate. The claim: because `Delta_lo` is both (a) the resource the certificate is built on and (b) the
thing inflation must conservatively shrink, padding to fix coverage *cannibalizes the same margin
that keeps the bracket finite at all* — producing a non-monotonic (rise-then-fall) coverage curve in
z, a genuinely different failure shape from the two prior noise specs. False if coverage turns out
monotonic in z here too (same as the prior two), or if the ceiling never actually bites at realistic
shot budgets.

## 2. Background and honest framing

- Builds on `certified_dipole.md`'s Davis-Kahan bracket (out-of-scope there: "shot-noise statistics
  on sigma_0, sigma_A"), and directly extends `certified_noise.md`'s single-bracket noise model to a
  THIRD, more heavily composed quantity, after `certified_thermochem_noise` (difference of two) and
  `gap_selfcheck_noise` (intersection of many).
- **What you can claim if the gates pass:** on a system where `Delta_lo` has healthy margin (HeH+),
  a moderate z (z~1) restores >=0.9 coverage — inflation still *works* there — but pushing z further
  measurably *reduces* both the finite-bracket rate and the achieved coverage at low shot budgets,
  an inflation ceiling neither prior noise spec showed (both were monotonic in z, verified up to
  z=5-6). On a system where `Delta_lo` is already fragile (LiH, at this Krylov depth/shot range) the
  bracket is mostly vacuous regardless of z — inflation cannot rescue what it has no margin to work
  with.
- **What you cannot claim:** a repair (a follow-up would need a bracket construction that doesn't
  spend the same margin twice — e.g. deeper Krylov depth first, then inflate); the exact z at which
  the ceiling bites (system/shot-budget-dependent, R1); non-i.i.d. noise; active-space dipoles.
- **Reference:** the exact `mu_exact` used for coverage is this spec's own noiseless M=16 Krylov
  point estimate (not the dense-diagonalization FCI value `certified_dipole` ultimately targets) —
  isolating the *noise* question from the *depth-convergence* question `certified_dipole`'s own gates
  already own; a cheap, honest simplification stated up front, not discovered after the fact.

## 3. Approach

Reuse `certified_dipole.spectral_width` (unchanged) and `certified_noise`'s i.i.d. shot-noise model
and `certified_half_width(lambda, shots, z) = z*lambda/sqrt(shots)` padding convention, applied
post-hoc (the same safe, unambiguous recipe `gap_selfcheck_noise` settled on, after an earlier probe
for *that* spec found internal directional perturbation ambiguous for composed quantities — the same
risk applies here, with three composed quantities instead of two).

One noisy realization per trial of `(theta0, var0, theta1, var1)` (self-mode, exactly
`certified_gaps.gap_bracket`'s formula) gives raw `sigma_0`, `Delta_lo`; one noisy realization of
`(mu, <A^2>)` (using the dipole operator's own Pauli 1-norms `lambda_A`, `lambda_A2`, computed the
same way `certified_noise.hamiltonian_one_norms` computes `lambda_H`) gives raw `mu`, `sigma_A`.
Pad `Delta_lo` down and `sigma_0`, `sigma_A` up by `certified_half_width(lambda_H_or_A, shots, z)`;
if the padded `s = sigma_0/Delta_lo >= 1` the trial is vacuous (scored as not-covered, tracked
separately as `finite_frac` — the same honest convention `gap_selfcheck_noise.frac_empty` uses).
Report `half = 2*sigma_A*s + W*s^2 + hw_A` (the last term pads the point estimate `mu` itself).
Monte Carlo (>= 6000 trials, seeded) on HeH+ and LiH (the same two systems `certified_dipole`'s own
example already uses), M=16.

## 4. Public interface

```
certified_dipole_noise.operator_one_norms(op: SparsePauliOp) -> (lambda, lambda2)
certified_dipole_noise.dipole_noise_coverage(
    mh, a_op, m, shots, z=2.0, trials=6000, seed=0,
) -> dict   # coverage, finite_frac, mu_exact, lam_h, lam_a, ...
```

## 5. Acceptance criteria (validation gates)

- **G1 — raw (z=0) coverage is broken on both systems.** At shots in {1e4, 1e5, 1e6}, M=16: HeH+
  coverage < 0.6, LiH coverage < 0.05. *Measured: HeH+ 0.294/0.491/0.505; LiH 0.001/0.002/0.009.*
- **G2 — inflation DOES work on the healthy-margin system, at moderate z.** At shots in {1e5, 1e6},
  z=1.0 restores coverage >= 0.9 on HeH+. *Measured: 0.981 / 1.000.*
- **G3 — THE FINDING (definition of done): an inflation ceiling, not seen in either prior noise
  spec.** At shots=1e4 (the tight-budget regime where `Delta_lo`'s margin is itself noisy), HeH+'s
  finite-bracket rate is STRICTLY LOWER at z=3 than at z=1 — more inflation shrinks the fraction of
  trials where a bracket can be constructed at all, the opposite of `certified_thermochem_noise`
  and `gap_selfcheck_noise`, both monotonically non-decreasing in z. *Measured: finite_frac
  z=1: 0.728, z=3: 0.554 (and coverage falls in step: 0.728 -> 0.554).*
- **G4 — LiH boundary, recorded not fixed.** Across every tested z in {0, 1, 2, 3} and every shot
  count, LiH's finite-bracket rate stays < 0.3 — inflation cannot rescue a system whose `Delta_lo`
  margin is already this thin at M=16; a system-dependent boundary, exactly the R1 pattern of the
  two prior noise specs. *Measured: max finite_frac over the whole (z, shots) grid = 0.260.*

> Definition of done: **G3**. If a system is later found where the dipole certificate stays
> monotonic in z (no ceiling), that contradicts the generality of the mechanism and must be
> recorded as a boundary on G3, not folded in silently.

## 6. Implementation plan (test-first)

1. Write `tests/test_certified_dipole_noise_spec.py` encoding G1-G4 (RED — module doesn't exist).
2. Add `certified_dipole_noise.py`: `operator_one_norms` (mirrors `hamiltonian_one_norms`),
   `dipole_noise_coverage` (single noisy realization + post-hoc padding, reusing
   `certified_dipole.spectral_width` and `certified_gaps`'s self-mode formula unmodified).
3. `make gates` / targeted pytest to green; ruff clean.

## 7. Out of scope

- A bracket construction that avoids double-spending the `Delta_lo` margin (the natural follow-up
  the ceiling motivates).
- Comparison against the dense-diagonalization FCI dipole (uses the M=16 Krylov point estimate as
  ground truth instead — isolates noise from depth-convergence, stated up front in §2).
- Oracle e1, more than two systems, active-space dipole convention, non-Gaussian noise.

## 8. Caveats and risks

- **R1 — the ceiling's exact location is shot-budget- and system-dependent** (it appears at
  shots=1e4 but not yet within z<=3 at shots=1e5-1e6 for HeH+; LiH never leaves the vacuous regime
  in the tested range). The falsifiable claim (G3) is the *existence* of a non-monotonic regime at
  a stated, reproduced shot budget, not a universal optimal-z formula.
- **R2 — three composed noisy quantities is the most fragile construction in this noise-spec
  series**; an earlier probe using internal directional perturbation (rather than post-hoc padding)
  produced a spurious *strictly decreasing* coverage curve from a sign error (theta0 is negative,
  so naive "shift down is conservative" reasoning inverted) — caught before it became a claim. The
  post-hoc padding recipe used here is deliberately the simplest, least assumption-laden choice.

## 9. Deliverables

- `certified_dipole_noise.py` — `operator_one_norms`, `dipole_noise_coverage`.
- `tests/test_certified_dipole_noise_spec.py` — gates G1-G4.
- Results summary (with R1/R2 caveats) in the PR description / BACKLOG entry.
