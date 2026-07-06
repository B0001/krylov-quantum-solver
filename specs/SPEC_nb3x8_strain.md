# SPEC: Nb3X8 strain / pressure response — Grüneisen parameters in the hopping

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G3 (the sharp prediction).

---

## 1. Goal

With inter-layer hopping |t| as the compression knob (uniaxial strain increases the dimerization
overlap at ~fixed U₀, Us), the leading strain response d ln X/d ln|t| of the Nb₃X₈ dimer observables
is: (a) the **spin gap stiffens everywhere** (γ_J > 0), running monotonically from the atomic-limit
value 2 (Nb₃F₈) toward 1 (Nb₃I₈); (b) the **charge gap is non-monotonic** (a minimum at |t\*|) and
the halide family **straddles** it, so γ_gap changes sign F/Cl (<0) → Br/I (>0); (c) **Nb₃Cl₈ sits at
its charge-gap minimum**, giving a spin-charge-decoupled strain response. Falsifiable: a negative
γ_J, a monotonic charge gap, or a Nb₃Cl₈ with strongly strain-dependent charge gap would break it.

## 2. Background and honest framing

Strain-tunable magnetism in Nb₃Cl₈ is an active experimental target (e.g. arXiv:2601.14524), so the
signs and magnitudes here are predictions against the lab, not just internal numbers. The exact
dimer (cRPA params of arXiv:2501.10320) supplies γ = d ln X/d ln|t| in closed form (spin) and by ED
(charge).

- **What we can claim:** the spin-gap strain exponent and its atomic→metallic run (2→1), the
  charge-gap minimum and the halide straddle, and the Nb₃Cl₈ spin-charge decoupling — falsifiable,
  experimentally-relevant, and numbers the source paper did not report.
- **What we cannot claim:** |t| is the *sole* strain proxy — real strain also shifts U₀, Us,
  geometry, and in-plane couplings this isolated dimer omits. A linear-response (Grüneisen)
  statement, not a strain phase diagram; density-density only. Nb₃F₈'s γ_J = 2 is the analytic
  atomic-limit exponent (J ∝ t²), not a magnitude claim (its J ~ 0.05 meV is below neglected terms).

## 3. Approach

**References:** the closed-form dJ/dt (exact) checked against central finite differences; the exact
charge gap Δc = E(3)+E(1)−2E(2) by ED (`exact_charge_gap`); the analytic small-/large-t limits
γ_J → 2 / 1.

**Numeric result (LT-bulk):** γ_J = 2.000/1.888/1.780/1.518 (F/Cl/Br/I), closed form == finite
difference; charge-gap minimum |t\*| = 271/152/122/76 meV (falls F→I); γ_gap = −0.004/−0.017/+0.078/
+0.368 (sign straddle); Nb₃Cl₈ |t| = 136 ≈ |t\*| = 152, so γ_J/γ_gap ≈ 111 (decoupled). χ-max and
Schottky Grüneisen equal γ_J to < 0.02% for Cl/Br, deviate 3.4%/5.9% for the iodide (the E_s/J
charge-contamination boundary, consistent with the rest of the thread).

## 4. Public interface

Reuses `exact_charge_gap`, `dimer_exchange_analytic`, `chi_max_temperature`, `schottky_peak_temperature`.

```
nb3x8_strain.spin_gap_gruneisen(U0, t, Us) -> float          # gamma_J, closed form
nb3x8_strain.charge_gap_gruneisen(U0, t, Us, rel=1e-4) -> float   # gamma_gap, finite diff
nb3x8_strain.charge_gap_min_hopping(U0, Us, t_hi=400) -> float    # |t*| (meV)
nb3x8_strain (CLI)                                           # family Grüneisen table + finding
```

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_strain_spec.py` (test-first).

- **G1 — spin gap stiffens, runs 2→1.** γ_J > 0 for all; closed form == central FD (< 1e-3); γ_J
  strictly decreasing F→I; γ_J(F) = 2 (atomic limit), γ_J(I) ∈ (1, 1.6).
- **G2 — charge-gap minimum + sign straddle.** Δc(t) has an interior minimum at |t\*| (below both
  flanks); |t\*| falls F→I; γ_gap > 0.1 (I), > 0 (Br), < 0.05 (Cl, F) — the family straddles |t\*|.
- **G3 — Nb₃Cl₈ spin-charge decoupled (DEFINITION OF DONE).** |t| within 20% of |t\*|; γ_J > 1.5,
  |γ_gap| < 0.05, |γ_J/γ_gap| > 30 — strain moves the spin gap, not the Mott gap.
- **G4 — J-scale observables share γ_J.** χ-max and Schottky Grüneisen == γ_J to < 0.1% for Cl/Br;
  the iodide deviates > 1% and ≫ the chloride (the E_s/J charge-contamination boundary).

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_strain_spec.py` encoding G1–G4 (initially failing — no module).
2. `nb3x8_strain.py` composing the closed-form dJ/dt + ED charge gap into Grüneisen parameters.
3. `make gates` (own process; no block2).

## 7. Out of scope

- Full strain response including dU₀/dstrain, dUs/dstrain, geometry, in-plane couplings.
- A strain phase diagram / quantum critical point at |t\*| (this is linear response).
- Nb₃F₈ spin-gap magnitude (J below neglected terms).

## 8. Caveats and risks

- **R1 — |t|-only strain.** Mitigated by framing as the leading (dominant overlap) response and
  stating the omitted channels; the qualitative signs (spin stiffens, charge straddles) are robust.
- Reproduction-adjacent machinery (Grüneisen analysis); the value is the ab-initio signs/magnitudes
  and the Nb₃Cl₈ decoupling prediction, reported with the |t|-only caveat.

## 9. Deliverables

- `nb3x8_strain.py` — Grüneisen tools + CLI.
- `tests/test_nb3x8_strain_spec.py` — gates G1–G4.
- `specs/SPEC_nb3x8_strain.md` — this spec.
- `specs/BACKLOG.md` — entry with the finding recorded.
