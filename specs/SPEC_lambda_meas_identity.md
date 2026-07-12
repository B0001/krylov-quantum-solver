# SPEC: λ_meas must exclude the identity — re-scoring the published crossovers honestly

**Status:** IMPLEMENTED once gates green. Follow-up promised by
[`SPEC_shift_both_sides.md`](SPEC_shift_both_sides.md) §3; supersedes its "left in place" note.

---

## 1. Goal

`precision_cost.measurement_lambda` summed **every** Pauli coefficient of the qubit Hamiltonian,
identity included. The identity term is a constant of **zero variance** — it costs zero shots — so
the near-term shot count N = (z·λ/ε)² was overstated for *raw and shifted Hamiltonians alike*, and
every published `precision_cost`/`cost_advisor` crossover number inherited the bias. Fix the metric
in place (exclude identity by default, `include_identity=True` retained for archaeology) and
re-score the published findings. Falsifiable: if excluding the identity did not move the crossovers
by the identity fraction, or if any *structural* finding of `SPEC_precision_cost` broke, the
correction (or the original) is wrong.

## 2. The numbers (H₂ / H₂O(4,3) / N₂(6,6), CASCI sto-3g, ε = 1.6 mHa, z = 2)

| quantity | H₂ | H₂O(4,3) | N₂(6,6) |
|---|---|---|---|
| λ_meas, identity included (old) | 2.699 | 9.925 | 22.844 |
| λ_meas, identity excluded (honest) | 1.887 | 5.400 | 14.297 |
| identity fraction of the 1-norm | 30.1% | 45.6% | 37.4% |
| resource ratio R at 1.6 mHa (old → honest) | 1.20e4 → 5.84e3 | 8.58e4 → 2.54e4 | 2.08e5 → 8.14e4 |
| crossover ε* (old → honest, ρ=1) | 19.1 → 9.34 | 137 → 40.6 | 332 → 130 |

## 3. THE FINDING — one qualitative verdict flips

Re-scoring is not just a constant rescale:

- **H₂ at ρ = 10⁴ flips FT → near-term.** With the inflated λ_meas the advisor called FT cheaper
  at chemical accuracy for a per-query cost 10⁴× a shot; honestly scored, near-term wins there.
  The identity bias was *hiding a regime where the near-term method is the right answer.*
  (H₂O and N₂ verdicts are stable at ρ ∈ {1, 10², 10⁴, 10⁶}.)
- **Every structural finding of `SPEC_precision_cost` survives** re-scoring: the exponent gap is
  ε-scaling (identity-free by construction); raw λ_DF still exceeds the honest λ_meas for N₂
  (24.94 > 14.30 — *stronger* now); the shifted λ_DF still beats λ_meas for every molecule with a
  margin that grows with size — but the margin is **1.94× → 2.95× → 3.58×**, not the published
  2.8×/5.4×/5.7×. `test_precision_cost_spec` G3's `> 5.0` bar is revised to `> 3.0` and this spec
  records why — the old bar measured identity mass, not FT advantage.
- The BACKLOG conjecture "no qualitative verdict flips" was **wrong** (H₂ above) — recorded here,
  gated in G3.

## 4. Public interface

```
precision_cost.measurement_lambda(mh, include_identity=False) -> float   # honest by default now
```

(Signature change is the fix; `shift_both_sides.shot_lambda` remains the standalone equivalent and
`SPEC_shift_both_sides` G5 still pins the shift-scoring artifact independently.)

## 5. Acceptance criteria (validation gates)

`tests/test_lambda_meas_identity_spec.py`.

- **G1 — the metric is honest by default.** `measurement_lambda(mh)` excludes the identity
  (matches `shot_lambda`); `include_identity=True` reproduces the old value; the identity fraction
  is material (> 25% for all three molecules — this was never a rounding error).
- **G2 — the crossovers move by the identity fraction.** ε* and R scale as the λ ratio (algebra:
  ε* ∝ λ_meas², R ∝ λ_meas²) — old/honest ratios match (λ_old/λ_new)² to 1e-9.
- **G3 — the verdict flip (the finding).** H₂ at ρ=10⁴, ε=1.6 mHa: `advise` says FT with the old
  metric and near-term with the honest one. N₂ and H₂O verdicts are unchanged at
  ρ ∈ {1, 10², 10⁴, 10⁶}.
- **G4 — the structural findings survive.** Under the honest metric: raw λ_DF > λ_meas for N₂
  still holds; shifted λ_DF < λ_meas for all three; the shift margin still grows with size
  (H₂ < H₂O < N₂, N₂ > 3.0).

## 6. Out of scope

- Grouped/commuting measurement schemes (the 1-norm is the ungrouped shot model, as before).
- The z-score and i.i.d. shot model themselves (unchanged from `certified_noise`).

## 7. Caveats and risks

- **R1 — committed literals move.** `test_precision_cost_spec` G3's N₂ bar 5.0 → 3.0 (recorded in
  §3); the SPEC_precision_cost numeric block gains a revision note. No other committed gate pins
  the inflated values (`cost_advisor`'s 22.844 is a hard-coded algebra fixture, still valid as a
  fixture).
- **R2 — history.** `SPEC_shift_both_sides` measured its shift gains with `shot_lambda` from the
  start, so its gates and numbers are untouched by this fix.

## 8. Deliverables

- `precision_cost.py` metric fix + docstring re-score; revised `tests/test_precision_cost_spec.py`
  G3; `tests/test_lambda_meas_identity_spec.py`; SPEC_precision_cost revision note; updated
  `SPEC_shift_both_sides.md` §3 pointer; this spec; BACKLOG close-out.
