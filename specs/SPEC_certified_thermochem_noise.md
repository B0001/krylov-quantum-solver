# SPEC: The certified relative-energy bracket under shot noise — composition beats either endpoint alone

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`certified_thermochem.py` composes two exact-statevector Temple/Ritz brackets (one per geometry)
into a rigorous interval on a relative energy Delta = E(B) - E(A). `certified_noise.py` showed that
a *single* such bracket's raw (uninflated) coverage of the true energy collapses to ~0.40 under shot
noise — a coin flip, not a bound — because the Ritz value sits so close to the true energy that
symmetric noise pushes it below E0 half the time. This spec asks the composed question: **does the
same collapse happen to the difference of two independently-noisy brackets, and does the same z=2
inflation rule (calibrated on one bracket) still have to pay for it?** The claim is false if the
composed bracket needs the *same or more* inflation than a single endpoint to reach 90% coverage.

## 2. Background and honest framing

- Builds on `certified_thermochem.md` (exact-statevector composition, explicitly out-of-scope:
  "shot-noise statistics on the brackets") and `certified_noise.md` (single-bracket coverage under
  i.i.d. shot noise). This is the direct intersection of the two, reusing both without new physics.
- **What you can claim if the gates pass:** relative energies (reaction/dissociation energies —
  chemistry's actual currency) inherit the certified_noise breakage, but composition partially
  self-corrects it: two independent one-sided coin-flips combined into a two-sided difference
  interval are measurably better calibrated than either coin-flip alone, so restoring 90% coverage
  costs *less* inflation than the single-bracket z=2 rule.
- **What you cannot claim:** a derivation of *why* composition helps beyond the empirical
  measurement (no closed-form proof of the reduced-z number); a rigorous (non-probabilistic) bound
  under sampling (impossible from finite samples, same as `certified_noise`); real grouped-Pauli
  measurement covariances (i.i.d. lambda-1-norm idealization, same as `certified_noise`); systematic
  (Trotter/basis) bias (statistical noise model only).
- **Reference:** the exact in-basis relative energy `mh_b.ground_state_energy() -
  mh_a.ground_state_energy()` (dense diagonalization) — the same reference `certified_thermochem`
  already validates its noiseless bracket against.

## 3. Approach

Reuse `certified_thermochem`'s two-geometry composition and `certified_noise`'s i.i.d. shot-noise
model, applied independently at each endpoint (different geometries -> independent noise draws,
different lambda_H/lambda_H2 1-norms). For each geometry, sample noisy `(rho, tau)` exactly as
`certified_noise.shot_noise_coverage` does (`th ~ N(th0, se_h)`, `h2 ~ N(h2_exact, se_h2)`,
`tau = th - var/(eps1-th)` when `eps1 > th` else `-inf`), each geometry's total energy including its
own `energy_offset` (the offsets differ across geometries and do **not** cancel in Delta, unlike in
`certified_gaps`/`certified_noise` where a single-geometry gap cancels the offset — this repo's
first cross-geometry noise composition, so this is a real new seam, not a copy-paste).

Compose: `Delta_lower = tau_B - rho_A`, `Delta_upper = rho_B - tau_A` (noiseless composition already
validated by `certified_thermochem`); apply inflation `z*se` at each endpoint independently before
composing. Reference: exact relative energy on the H4 symmetric stretch (0.9 -> 2.3 A), the same
system `certified_thermochem` uses, oracle `E1` per geometry (isolates the noise-composition
question from the self-mode Temple-premise question `certified_gaps`/`gap_selfcheck` already own).

Monte Carlo over noise realizations (>= 20000 trials, seeded) measures coverage of the exact Delta,
same style as `certified_noise.shot_noise_coverage`.

## 4. Public interface

```
certified_thermochem_noise.thermochem_noise_coverage(
    mh_a, mh_b, m, shots, e1_a=None, e1_b=None, trials=20000, z=2.0, seed=0,
) -> dict   # cov_raw, cov_inflated, min_z_for_90pct, lam_h_a, lam_h_b, ...
certified_thermochem_noise.minimal_z_for_coverage(mh_a, mh_b, m, shots, target=0.9, ...) -> float
```

## 5. Acceptance criteria (validation gates)

- **G1 — composition inherits the collapse, and it is shot-count-independent.** Raw (z=0) coverage
  of Delta on the H4 stretch is < 0.85 (broken) at shots in {1e4, 1e5, 1e6}, all M=16, oracle E1;
  and the spread across those three shot counts is < 0.05 (N-independence, mirroring
  `certified_noise`'s core finding). *Measured: raw coverage 0.727 / 0.700 / 0.688 (spread 0.039).*
- **G2 — the existing single-bracket rule (z=2) still works when reused as-is.** Applying the
  `certified_noise` z=2 inflation independently at each endpoint before composing restores coverage
  >= 0.9 (sanity: reuse doesn't break). *Measured: 0.999-1.000 at all three shot counts.*
- **G3 — THE FINDING (definition of done): composition needs LESS inflation than either endpoint
  alone.** The minimal z (grid search, 0.05 resolution) that restores composed coverage >= 0.90 is
  strictly less than 2.0 (the single-bracket number) at shots=1e5, M=16. *Measured: z* = 0.60
  (shots=1e4: z*=0.55) -- roughly 3-4x less inflation than the single-endpoint rule.*
- **G4 — mechanism check.** At z=0, the composed bracket's raw coverage exceeds EITHER single
  endpoint's raw coverage (from `certified_noise.shot_noise_coverage` on `mh_a` and `mh_b`
  separately, same M/shots/seed) by a wide margin, not a rounding-level one: composed >= single + 0.2.
  *Measured: composed 0.700 vs endpoint-A 0.403 / endpoint-B 0.413 at shots=1e5 -- margin ~0.29.*

> Definition of done: **G3**. If bisection ever lands z* >= 2.0 at some tested condition the
> composition-helps claim is false there and must be recorded as a boundary, not silently dropped.

## 6. Implementation plan (test-first)

1. Write `tests/test_certified_thermochem_noise_spec.py` encoding G1-G4 (RED — module doesn't exist).
2. Add `certified_thermochem_noise.py`: `thermochem_noise_coverage` (mirrors
   `certified_noise.shot_noise_coverage`'s noise model, applied per-endpoint with
   `mh.energy_offset` added before composing) and `minimal_z_for_coverage` (bisection over z).
3. `make gates` / targeted pytest to green; ruff clean.

## 7. Out of scope

- A closed-form explanation of *why* the required z shrinks (empirical only).
- Non-oracle (self) eps1 mode — that premise question belongs to `certified_gaps`/`gap_selfcheck`;
  this spec isolates the noise-composition question with an oracle E1 per geometry.
- Grouped-Pauli measurement covariances; systematic bias; more than two composed geometries.
- Any claim about `certified_dipole`/`certified_gaps` under noise (separate compositions).

## 8. Caveats and risks

- **R1 — the "composition helps" number may be specific to this system.** The H4 stretch has very
  different lambda_H/lambda_H2 at each endpoint (10.26/62.9 eq vs 6.84/29.9 stretched) and a
  near-degenerate stretched-endpoint gap; a system with symmetric endpoints might show a different
  margin. Mitigation: G4 is the mechanism gate, kept system-relative rather than an absolute z*
  number, so the claim (composition > either alone) is the falsifiable part; the specific z* ~ 0.5
  is reported as a measurement, not a universal constant.
- Honest limitation: two geometries only, oracle gap, i.i.d. Gaussian shot noise — same idealization
  chain as `certified_noise`.

## 9. Deliverables

- `certified_thermochem_noise.py` — `thermochem_noise_coverage`, `minimal_z_for_coverage`.
- `tests/test_certified_thermochem_noise_spec.py` — gates G1-G4.
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
