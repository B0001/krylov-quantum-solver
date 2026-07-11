# Task Breakdown 10 — #15 Shots-to-Certainty Planner
Goal: two-stage (pilot → allocate) measurement budgeter validated on the golden suite. Depends: CertChem-M1 (uses the estimator + error bars).

1. **Sensitivity map** — from the single-signal error-bar machinery, derive/implement d(bracket width)/d(variance at time point t_k) — which samples buy the most certainty.
   ✓ Synthetic-signal test: sensitivity ranking matches brute-force leave-one-out. (M)
2. **Forecast model** — width(shots per point) predictor from pilot-run noise estimates; state assumptions (shot-noise-dominated regime).
   ✓ On simulated data: predicted vs achieved width within stated tolerance across a 10× shot-budget range. (M)
3. **Two-stage planner** — API: `plan(target_width, pilot_data) -> allocation`; hard-restrict v1 to two-stage (no mid-run adaptivity) so bracket validity is untouched — document WHY (adaptive bias risk) in the module docstring.
   ✓ Allocation sums to budget; monotone: tighter target → more shots. (S)
4. **Validity check** — statistical test that two-stage allocation preserves bracket coverage: 500 simulated repetitions, count bracket-containment frequency vs nominal.
   ✓ Coverage ≥ nominal within binomial error. This is the certifying experiment. (L)
5. **Golden-suite validation** — for each golden system: request widths {10, 5, 2, 1.6, 1} mHa; plot predicted vs achieved shots.
   ✓ Honest plot published, including where the forecast breaks down. (M)
6. **Integrate + document** — expose as `certified_energy(..., target_width=)` convenience; costs auto-planned.
   ✓ One-line UX; docs show a before/after shot-count saving on N₂. (S)
