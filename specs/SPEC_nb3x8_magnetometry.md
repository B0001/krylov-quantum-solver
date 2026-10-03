# SPEC: Nb3X8 interlayer exchange vs measured magnetometry — a parameter-free prediction

**Status:** IMPLEMENTED — gates G1–G5 green. The definition-of-done gate (G3) records the finding;
G5 (bd chem-jiy) records a primary-source correction to the θ_W input and kills a proposed
phase-assignment correction to it.

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
| Nb₃Cl₈ | 66.2 | 479 | 90 | 5.3× | −192 vs −51.2 (3.75×) |
| Nb₃Br₈ | 119.1 | 862 | 382 | 2.3× | — |
| Nb₃I₈ | 245.9 | 1756 | (no transition) | — | — |

**Correction (bd chem-jiy, 2026-09-30):** θ_W was originally quoted as −13.1 K (15× miss). That
number does not appear anywhere in the cited primary source and reached this repo via an unverified
search-engine summary. A primary-source read of `arXiv:1701.05528` p.4 gives: *"an analysis of the
inverse susceptibility data for T > 140 K yields a Curie constant C = 0.484 emu·K·mol⁻¹·Oe⁻¹
(p_eff = 1.97, consistent with S_eff = 1/2) and a Weiss temperature of θ = −51.2 K."* C and μ_eff
were already correctly quoted; only θ_W was wrong. The verified miss is 3.75×, not 15× — real, but
smaller than originally published. See G5 below.

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
- **G5 — the phase-correction hypothesis, tested and killed (bd chem-jiy).**
  `specs/BACKLOG.md`'s Nb₃X₈ entry conjectured that the G3(b) miss is a phase-assignment artifact:
  θ_W is fitted in the HT (undimerized) phase, but G3(b) compared it to a LT-phase J; swapping in
  the HT-phase J (`nb3x8_gaps.NB3X8_HT_BULK`, Table IV, already unused) was claimed to invert the
  miss to a clean mean-field z_eff = 4.1 (θ = −z·J/4), using the (wrong) θ_W = −13.1 K.
  Pre-registered kill conditions: z_eff unphysical (<1 or >12 for a layered stack), or χ_HT(200 K)
  misses the verified Curie constant by >2×. Recomputed with the **verified** θ_W = −51.2 K:
  `z_eff_ht("Nb3Cl8")` = 16.1 — **outside** the [1, 12] bound. The hypothesis is **killed**: it was
  an artifact of the wrong, unverified −13.1 K datum, not a real phase-assignment fix. (χ_HT(200 K)
  alone is within the 2× sanity bound — 1.53× — so it alone would not have killed the hypothesis;
  the z_eff bound is what kills it.)

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
  **Checked and FALSIFIED (chem-g78, `tests/test_nb3x8_gaps_spec.py::test_G8_*`):** running the
  identical coordination machinery in the spin channel (`nb3x8_gaps.coordination_spin_gap`) gives
  J_eff = 66.20 (z=0) → 71.08 (z=1) → 66.55 (z=2) → 71.12 meV (z=3) for Nb₃Cl₈ and 119.11 → 126.50
  → 120.11 → 126.63 meV for Nb₃Br₈ (Br measured at landing; the original chem-g78 run covered Cl
  only), i.e. it oscillates and is *larger*, not smaller, at the largest cluster reached — nowhere
  near the J₀/3 ≈ 22 meV (Cl) / J₀/2.26 ≈ 53 meV (Br) a rescuing reduction would need — while the SAME
  machinery's charge-channel control on the SAME clusters gives the expected monotonic softening
  (33.5% Cl, 30.1% Br at z=3). So coordination/mean-field reduction, as modelled by finite open
  L≤8 clusters (no 3-D triplon band), is **not** the mechanism behind the 5.3×/2.3× Tc
  overprediction. What is left open: the cooperative structural transition, in-plane kagome
  exchange, or a direct lattice renormalization of t_s⊥ (candidates, none tested here) — this
  follow-up is closed as a negative result, not merely deferred.

## 8. Caveats and risks

- **R1 — experimental provenance.** The measured Tc/θ_W are quoted from the cited primary sources
  (Sheckelton 2017, Haraguchi 2017); Nb₃Br₈'s ~382 K is the bulk value. The gates use bounded
  ranges (order-of-magnitude scale, factor > 2 overprediction), not exact hits, so they are robust
  to modest revision of the experimental numbers while still being falsifiable. **This risk was
  realized**: θ_W = −13.1 K (an earlier version of this spec/module) was not actually in the primary
  source — see G5. Fixed 2026-09-30 by reading `arXiv:1701.05528` p.4 directly; C and μ_eff were
  already correct.
- **R2 — z_eff is extracted, not predicted (G5).** `z_eff_ht` solves θ = −z·J/4 for z given a
  measured θ_W — it is a diagnostic, gated only against a physical bound (1–12), never reported as
  an independent prediction. A z_eff that happens to land inside the bound would not, by itself,
  validate the phase-assignment hypothesis; it would only fail to falsify it.
- The overprediction is expected physics (a single dimer over-counts the coupling relative to the
  cooperative lattice transition), reported honestly — the value is the *quantification* and the
  two-coupling separation, not a claim of quantitative agreement.

## 9. Deliverables

- `nb3x8_magnetometry.py` — predictor + cited experimental table + CLI + `z_eff_ht`/`chi_ht_curie_ratio`.
- `tests/test_nb3x8_magnetometry_spec.py` — gates G1–G5.
- `specs/SPEC_nb3x8_magnetometry.md` — this spec.
- `specs/BACKLOG.md` — entry moved to done with the finding recorded (chem-jiy: hypothesis killed).
- `nb3x8_gaps.py` — `NB3X8_HT_BULK` (Cl/Br HT-phase cRPA parameters, Table IV), factored out for
  reuse (previously only embedded inside `NB3X8_CLUSTERS`).
