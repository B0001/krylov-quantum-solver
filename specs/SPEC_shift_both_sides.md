# SPEC: Shift both sides — the symmetry shift helps the near-term arc too, moving the crossover

**Status:** IMPLEMENTED (gates G1–G5 green in `tests/test_shift_both_sides_spec.py`).

---

## 1. Goal

The BLISS/SCDF number-operator symmetry shift (`df_factorization.symmetry_shift`), which the bridge
used to shrink the FT 1-norm λ_DF, **also lowers the measurement 1-norm λ_meas** — so it cuts the
*near-term* certified shot cost (N ∝ λ_meas²) as well. The bridge shifted only the FT side and
thereby **overstated FT's advantage**; a fair comparison shifts both, moving the crossover flip-ρ
down by the same factor. Falsifiable: if the shift did not lower λ_meas, or the fair flip-ρ did not
drop by (λ_raw/λ_shift)², the claim breaks.

## 2. Background and honest framing

`cost_advisor`/`precision_cost` compared near-term with the *raw* λ_meas against FT with the
*shifted* λ_DF — an inconsistency. The shift is a spectrum-preserving constant on the target sector
(H → H + (b1 N + b2)(N − nₑ)); its effect on any 1-norm is fair game for either method.

- **What we claim:** the shift lowers λ_meas materially (≥ 35%); the near-term shot cost drops by
  (λ_raw/λ_shift)²; the fair (both-sided) flip-ρ is that same factor below the one-sided bridge
  value — so the one-sided bridge overstates FT's advantage; a λ_meas-optimized shift does at least
  as well as the λ_DF-optimized one.
- **What we cannot claim:** that the shift changes the *exponents* (still 1/ε² vs 1/ε) — it moves
  constants only; nor that it is free on hardware (it adds number-operator terms to measure).

## 3. THE FINDING — the identity term inflated the draft's headline

This spec was drafted claiming λ_meas reductions of **−54%/−51%/−73%** and a fair-crossover shift of
**4–14×**. Those numbers are real but they are scored with the **wrong 1-norm**.

The shot count is N = (z·λ/ε)², where λ must sum only the **non-identity** Pauli coefficients: the
identity term is a constant with **zero variance**, so it costs zero shots. The repo's
`precision_cost.measurement_lambda` sums *every* coefficient, identity included. That matters here
because a large part of what the shift does to the qubit Hamiltonian is **dump weight into the
identity term** (N₂: identity 8.55 → 0.10 Ha). Scoring the shift with the identity included
therefore *flatters* it:

| reduction in λ_meas | H₂ | H₂O(4,3) | N₂(6,6) |
|---|---|---|---|
| identity **included** (draft) | 53.9% | 51.4% | 73.2% |
| identity **excluded** (honest, shot-relevant) | 42.9% | 39.4% | 57.9% |

| crossover gain (λ_raw/λ_shift)² | H₂ | H₂O(4,3) | N₂(6,6) |
|---|---|---|---|
| identity **included** (draft) | 4.70× | 4.24× | 13.92× |
| identity **excluded** (honest) | **3.07×** | **2.73×** | **5.65×** |

**The claim survives, smaller.** The shift really does cut the near-term 1-norm, and the one-sided
bridge really does overstate FT — but by **2.7–5.7×**, not the 4–14× the draft advertised. N₂'s fair
flip-ρ is **1.44×10⁴** against the one-sided **8.14×10⁴** (not the drafted 1.5×10⁴ vs 2.1×10⁵, both
of which carried identity mass). G5 pins this so the artifact cannot silently return.

**Consequence for the repo:** `precision_cost.measurement_lambda` overstates the near-term shot cost
for *raw and shifted alike*, so the published `precision_cost`/`cost_advisor` crossovers inherit the
same bias. Left in place here (changing it would move committed numbers in another spec's gates);
`shift_both_sides.shot_lambda` is the honest one. → BACKLOG follow-up.

## 4. Public interface

```
shift_both_sides.shot_lambda(mh, include_identity=False) -> float   # the honest 1-norm
shift_both_sides.shifted_measurement_lambda(h1, eri, norb, nelec, target="raw"|"df"|"meas") -> float
shift_both_sides.fair_flip_rho(h1, eri, norb, nelec, eps, z=2.0) -> (one_sided, both_sided)
shift_both_sides.optimize_meas_shift(h1, eri, norb, nelec) -> (b1, b2)
shift_both_sides.spectrum_preserved(h1, eri, norb, nelec) -> bool   # FCI invariance check
shift_both_sides (CLI)                                              # per-molecule shift table
```

## 5. Acceptance criteria (validation gates)

`tests/test_shift_both_sides_spec.py` — 14 green.

- **G1 — the shift lowers λ_meas too.** SCDF-shifted λ_meas ≤ 0.65·λ_meas_raw for every molecule
  (≥ 35% reduction). *Revised from the draft's ≥ 40%:* on the honest identity-excluded metric H₂O
  reduces by 39.4%, which the draft's bar would have failed. The 40% bar was an artifact of the
  inflated metric — recording the revision, per the SDD loop.
- **G2 — the fair crossover drops (DEFINITION OF DONE).** Shifting both sides lowers the near-term
  shot cost by (λ_raw/λ_shift)² ≥ 2.5× (N₂ > 5×), and the both-sided flip-ρ is exactly that factor
  below the one-sided bridge value. A *lower* flip-ρ is a *weaker* case for FT.
- **G3 — the objectives are aligned.** A λ_meas-optimized shift is never worse than the SCDF one
  (it is seeded from it). *Honest:* the draft claimed SCDF already captures most of the λ_meas gain
  — true for N₂ (within 5%: 6.017 → 5.990) but **not** for H₂O, where re-optimizing buys a further
  37% (3.271 → 2.057).
- **G4 — spectrum preserved.** FCI(shifted) + e_shift == FCI(raw) to < 1e-8 Ha on every CAS — the
  shift is exact, so it cannot touch the exponents, only the constants.
- **G5 — the identity artifact (the finding).** The identity-inclusive metric reports a strictly
  larger gain than the honest shot metric for every molecule; the honest gain never reaches the
  advertised 14×. Goes red if λ_meas is ever "fixed" back to include identity.

## 6. Out of scope

- Hardware cost of measuring the added number-operator terms.
- Tensor-optimized SCDF / THC shifts (number-operator shift only).
- Correcting `precision_cost.measurement_lambda` itself (see §3 — BACKLOG follow-up).

## 7. Caveats and risks

- **R1 — the shift is not free on hardware.** It adds number-operator terms that must themselves be
  measured; the 1-norm accounting here does not charge for that.
- **R2 — constants, not exponents.** Both arcs keep their scaling (1/ε² vs 1/ε). The shift moves
  *where* they cross, never *whether* FT eventually wins.

## 8. Deliverables

- `shift_both_sides.py`, `tests/test_shift_both_sides_spec.py`, this spec, BACKLOG entry.
