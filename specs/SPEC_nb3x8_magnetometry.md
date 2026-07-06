# SPEC: Nb3X8 interlayer exchange vs measured magnetometry — a parameter-free prediction

**Status:** IMPLEMENTED — gates G1–G4 green. The definition-of-done gate (G3) records the finding.

---

## 1. Goal

The ab-initio interlayer singlet–triplet gap J of the Nb₃X₈ bilayer dimers (extracted in
`odmd_spin` / `nb3x8_susceptibility` from the cRPA parameters of `arXiv:2501.10320`), with **no fit
parameters**, predicts the *scale* and *ordering* of the measured magnetic-singlet transitions of
the real solids Nb₃Cl₈ (~90 K) and Nb₃Br₈ (~382 K); and where it misses, the miss is quantified and
localizes the missing physics. Falsifiable: a wrong J, or a J unrelated to the transition, would
produce the wrong ordering or a scale off by orders of magnitude.

## 2. Background and honest framing

The low-temperature nonmagnetic state of Nb₃Cl₈/Nb₃Br₈ **is** the interlayer dimerization of
adjacent-layer Nb₃ (S=½) clusters into singlets (Sheckelton et al., *Inorg. Chem. Front.* **4**, 481
(2017), `arXiv:1701.05528`; Haraguchi et al., *Inorg. Chem.* **56**, 3483 (2017)). So the coupling J
that `odmd_spin` computes from the downfolded bilayer is exactly the coupling that sets the
singlet-formation temperature. A Heisenberg/Bleaney–Bowers dimer's susceptibility peaks at
k_BT_max ≈ 0.625 J; that T_max is the ab-initio prediction of the transition scale.

- **What we can claim:** a parameter-free, experiment-referenced prediction of the *scale and
  ordering* of the Nb₃Cl₈/Nb₃Br₈ transitions from the downfolded interlayer J, plus a quantified,
  physically-interpreted account of where it breaks (a number the cluster papers did not report).
- **What we cannot claim:** the isolated bilayer dimer has no in-plane kagome exchange, no phonons,
  and no cooperative/first-order structural transition — it sets *scales*, it does not reproduce a
  first-order transition. Nb₃I₈ has no interlayer-singlet transition (a moment-retaining ground
  state) and is excluded from the comparison. This is a comparison against measured values, **not a
  fit**; the experimental numbers are quoted with primary-source citations.

## 3. Approach

**Reference = experiment** (the measured transition temperatures and Weiss temperature), plus the
analytic Bleaney–Bowers dimer for the machinery anchor (G1).

For each halide, compute J with the exact `dimer_exchange_analytic`, the exact dimer χ(T) with
`susceptibility` (`nb3x8_susceptibility`), and T_max as the maximizer of χ(T). Convert to kelvin and
compare to the measured Tc; compare θ_CW = −J/4 to the measured Curie–Weiss θ_W.

**Numeric result (LT-bulk cRPA parameters):**

| system | J (meV) | pred. χ-max (K) | obs. Tc (K) | overpred. | −J/4 (K) vs obs. θ_W |
|--------|--------:|----------------:|------------:|----------:|---------------------:|
| Nb₃Cl₈ | 66.2 | 479 | 90 | 5.3× | −192 vs −13.1 (15×) |
| Nb₃Br₈ | 119.1 | 862 | 382 | 2.3× | — |
| Nb₃I₈ | 245.9 | 1756 | (no transition) | — | — |

## 4. Public interface

Reuses validated primitives (`susceptibility`, `curie_weiss_theta`, `dimer_exchange_analytic`,
`NB3X8_LT_BULK`); the only new code is the predictor + the cited experimental table.

```
nb3x8_magnetometry.chi_max_temperature(U0, t, Us) -> float        # meV, exact chi(T) maximum
nb3x8_magnetometry.predicted_transition_K(U0, t, Us) -> float     # K, predicted transition scale
nb3x8_magnetometry.overprediction_factor(name) -> float           # pred / measured Tc
nb3x8_magnetometry.theta_over_measured(name) -> float             # (-J/4) / measured theta_W
nb3x8_magnetometry.EXPERIMENT                                      # cited measured references
nb3x8_magnetometry (CLI)                                          # family table + finding
```

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_magnetometry_spec.py` (test-first).

- **G1 — machinery anchored.** χ·T → ½ (S=½ pair Curie constant) deep in the spin window,
  θ_CW = −J/4 exactly, and the exact χ maximum sits within 5% of 0.625 J for every J-resolvable
  member — so `chi_max_temperature` is a faithful transition-scale estimator.
- **G2 — measured scale and ordering (the win).** For Nb₃Cl₈ and Nb₃Br₈ the parameter-free
  prediction lands within an order of magnitude of the measured Tc (`1 < pred/obs < 10`) and
  reproduces the observed Cl < Br ordering.
- **G3 — the finding (DEFINITION OF DONE).** (a) The isolated dimer overpredicts Tc for both halides
  (factor > 2) and the overprediction weakens monotonically Cl → Br (the isolated-cluster →
  cooperative-lattice renormalization). (b) −J/4 overshoots the measured θ_W of Nb₃Cl₈ by > 5×, so
  the interlayer J that sets Tc is a *different* coupling from the weak in-plane exchange that sets
  θ_W.
- **G4 — honest boundary.** Nb₃I₈ (no interlayer-singlet transition) is excluded from the
  experimental set; the predictor is a positive, monotone-in-J *scale* estimator, never a claim to
  reproduce a first-order cooperative transition.

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_magnetometry_spec.py` encoding G1–G4 (initially failing — no module).
2. `nb3x8_magnetometry.py` composing the validated χ(T)/J primitives + the cited experimental table.
3. `make gates` (own process; no block2).

## 7. Out of scope

- In-plane kagome exchange, the cooperative/first-order structural transition, phonons.
- Nb₃I₈ magnetometry (different ground state).
- Any fit to the experimental data (this is a prediction/comparison).
- A quantitative theory of the isolated→solid renormalization factor (a follow-up: coordination /
  mean-field reduction of T_max, cf. the `nb3x8_gaps` coordination-scan treatment of the charge gap).

## 8. Caveats and risks

- **R1 — experimental provenance.** The measured Tc/θ_W are quoted from the cited primary sources
  (Sheckelton 2017, Haraguchi 2017); Nb₃Br₈'s ~382 K is the bulk value. The gates use bounded
  ranges (order-of-magnitude scale, factor > 2 overprediction), not exact hits, so they are robust
  to modest revision of the experimental numbers while still being falsifiable.
- The overprediction is expected physics (a single dimer over-counts the coupling relative to the
  cooperative lattice transition), reported honestly — the value is the *quantification* and the
  two-coupling separation, not a claim of quantitative agreement.

## 9. Deliverables

- `nb3x8_magnetometry.py` — predictor + cited experimental table + CLI.
- `tests/test_nb3x8_magnetometry_spec.py` — gates G1–G4.
- `specs/SPEC_nb3x8_magnetometry.md` — this spec.
- `specs/BACKLOG.md` — entry moved to done with the finding recorded.
