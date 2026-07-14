# SPEC: Iterative QPE — a fat-tailed bit-flip cascade the median error hides

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`iterative_qpe.py`'s docstring claims "precision ~ 1/2^bits," never gated, and — unlike
`qpe_walk_readout.py`'s noiseless exact-simulation oracle — this module's per-bit measurement IS
genuinely stochastic (a majority vote over `shots_per_bit` Bernoulli trials). Iterative QPE reads
bits least-significant-first and feeds each one back as a rotation that cancels its contamination
for every subsequent (coarser) bit — so a single flipped LOW bit corrupts the feedback angle for
every bit measured afterward. This spec asks whether that architecture creates a fat tail: does the
MEDIAN error improving with more bits (as the docstring implies) hide a MAXIMUM (worst-case) error
that is far worse, at low per-bit shot counts? False if the median/max ratio stays close to 1 at
every tested `(n_bits, shots_per_bit)` combination — i.e., if there is no cascade risk to find.

## 2. Background and honest framing

- `iterative_qpe.py` already reuses `qubitization_blueprint`'s exact JW Hamiltonian for the
  validation-oracle eigenphase — no new physics, only falsifiers around a genuinely stochastic
  simulation the module runs (median-behavior only, in its `__main__`) but never characterizes by
  its tail.
- **What you can claim if the gates pass:** precision does improve monotonically (median error
  non-increasing) with more bits at the module's own default `shots_per_bit=15`; but at LOW
  `shots_per_bit` (3, tested) and LOW `n_bits` (8, tested) the worst-case error over seeds is
  **~10x the median** — a real fat tail from cascading LSB bit-flips, not visible if you only look
  at a typical run; and that tail VANISHES (max == median, fully deterministic) once
  `shots_per_bit` crosses a measured threshold (>=3 at `n_bits=12` on this system) — cheap
  insurance against the cascade, quantified rather than assumed.
- **What you cannot claim:** that the `shots_per_bit>=3` threshold generalizes to other systems or
  larger `lambda` (it depends on how close intermediate rotation angles land to the ambiguous
  `sin^2(theta/2)=0.5` point, a system-specific quantity); that this reproduces a real device's
  measurement statistics (each "shot" here is an exact-`p1` Bernoulli draw, no readout error, no
  decoherence during the feedback loop).
- **Reference:** the exact ground eigenphase (`np.linalg.eigvalsh`, dense) — the same reference
  `iqpe_ground_energy` itself already uses as its comparison oracle.

## 3. Approach

Reuse `iqpe_ground_energy` unmodified. For each `(n_bits, shots_per_bit)` combination, run >= 30
independently-seeded trials and record the median and maximum absolute error against CASCI. Compare
median trends across `n_bits` at fixed `shots_per_bit=15` (the module's own default); compare
median-vs-max at fixed `n_bits=8` across `shots_per_bit in {3, 15}`; sweep `shots_per_bit` at fixed
`n_bits=12` to find where max collapses onto median.

## 4. Public interface

No new library code — this spec adds only test-file assertions around `iterative_qpe.py`'s existing
`iqpe_ground_energy`, reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — precision improves monotonically with more bits, at the module's own default shots.** On
  H2 CAS(2,2), `shots_per_bit=15`, median error (30 seeds) is non-increasing across
  `n_bits = 4, 6, 8, 10, 12, 14, 16, 18, 20`. *Measured: 281.6 -> 85.2 -> 12.9 -> 0.66 -> 0.66 ->
  0.106 -> 0.085 -> 0.0105 -> 0.0015 mHa — a staircase, same character as `qpe_walk_readout`'s
  precision law, never increasing.*
- **G2 — THE FINDING (definition of done): a fat-tailed cascade at low shots_per_bit.** At
  `n_bits=8`, `shots_per_bit=3` (40 seeds): `max_err / median_err > 5`. At the SAME `n_bits=8` but
  the module's own default `shots_per_bit=15`: the ratio is exactly 1 (fully deterministic across
  all 40 seeds — the default already avoids the pathology). *Measured: ratio 10.39 at spb=3;
  median==max at spb=15 (12.9325 mHa both).*
- **G3 — the tail has a measured threshold, not an open-ended risk.** At `n_bits=12`,
  `shots_per_bit=1`: max/median > 5 (cascade present). At `shots_per_bit >= 3`: max == median
  exactly across 40 seeds (deterministic, tail gone). *Measured: spb=1 ratio 10.29; spb=3,5,7,9 all
  deterministic.*
- **G4 — boundary, recorded not smoothed over: the module's own smallest demonstrated bit counts
  never reach chemical accuracy.** At the module's own default `shots_per_bit=15`, `n_bits=8`
  median error exceeds 10 mHa (chemical accuracy) on every one of 40 seeds — the docstring's
  `~1/2^bits` claim is true in direction but the `__main__` demo's own printed `n_bits=4,6,8` rows
  are honestly still in the gross-error regime; only `n_bits >= 10` crosses into accuracy here.

> Definition of done: **G2**. G1 alone (median improves) is exactly the kind of check that would
> miss the tail — G2 is what actually falsifies "1/2^bits is the whole precision story."

## 6. Implementation plan (test-first)

1. Write `tests/test_iqpe_cascade_risk_spec.py` encoding G1-G4 (RED in the sense that these checks
   are new, even though `iterative_qpe.py`'s functions are not).
2. No changes to `iterative_qpe.py` — a genuine external characterization of behavior the module
   already exhibits (stochastically, per-run) but never summarizes across seeds.
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- A fix for the cascade (e.g. bit-flip-tolerant feedback, repetition-coded bit reads) — this spec
  measures the risk, a follow-up would need to design and gate a mitigation.
- Systems beyond H2 CAS(2,2), or a device noise model (readout error, decoherence during the
  feedback loop) beyond the module's own exact-`p1`-Bernoulli abstraction.
- Deriving the `shots_per_bit` threshold analytically (measured only, per R1).

## 8. Caveats and risks

- **R1 — the measured `shots_per_bit` threshold (3, at n_bits=12) is a property of THIS system's
  rotation-angle sequence**, not a universal constant; a system whose intermediate angles land
  closer to `sin^2(theta/2)=0.5` would need more shots per bit for the same margin. The falsifiable
  claim (G2/G3) is that the fat tail EXISTS at low shots and VANISHES above a measured threshold on
  at least one system, not a portable formula.
- Honest limitation: exact-`p1` Bernoulli simulation, no device noise beyond that; one system.

## 9. Deliverables

- `tests/test_iqpe_cascade_risk_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
