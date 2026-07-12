# Task Breakdown 2 — ChemCheck M1–M2: Tiers + Mode A Scorer
Goal: frozen tier definitions and a spec-sheet → scorecard pipeline calibrated against noisy simulation.

1. **Freeze tier registry** — `tiers.py`: for T0–T3, serialize Hamiltonian (canonical term stream + SHA-256), FCI reference, accuracy thresholds (1.6/16 mHa), resource counts from `benchmark_resources.py`. Version string `chemcheck-2026.1`.
   ✓ Registry loads; hashes stable across machines; T4 stub marked aspirational. (M)
2. **Submission validator** — implement `chemcheck-submission.schema.json` validation (jsonschema lib); reject with pointer-precise errors.
   ✓ 3 valid + 6 invalid fixture files behave correctly. (S)
3. **Routing/overhead model** — map each tier circuit onto connectivity classes (all-to-all, heavy-hex, grid, linear): SWAP-overhead multiplier per class. Start with published lookup factors; document as v1 crudeness.
   ✓ Overhead(all_to_all)=1.0; heavy-hex > grid > all-to-all ordering holds; factors documented with sources. (M)
4. **Error-budget model** — expected total error = f(two-qubit count × routed overhead × error rate, decoherence vs depth/T₂). Depolarizing-only v1, stated loudly.
   ✓ Pure function, unit-tested against hand-computed cases. (M)
5. **Headroom computation** — invert the budget: required two-qubit error to hit PASS threshold per tier; headroom = current/required.
   ✓ Monotonicity tests: better error → smaller headroom; headroom=1.0 exactly at threshold. (S)
6. **Calibration harness (the crux)** — run tier circuits through the repo's simulator with injected depolarizing noise at swept rates; find empirical pass/fail transition; compare against Mode A prediction.
   ✓ Predicted transition within stated model uncertainty for T0–T2; discrepancy for T3 documented, not hidden. (L)
7. **Floor-detector test** — feed the Mode B scorer skeleton the old codebase's known-bad energies.
   ✓ 100% flagged UNPHYSICAL; zero false positives on golden results. (S)
8. **Scorecard emitter** — JSON per `chemcheck-scorecard.schema.json` + rendered Markdown table with the `classically_simulable` disclaimer on T0–T2.
   ✓ Schema-validates; renders readable on GitHub. (S)
9. **Launch artifact** — Mode-A-score 3 real devices from published spec sheets; write "State of Quantum Chemistry Hardware" report with headroom factors as headline numbers.
   ✓ Report cross-checked: no device certified for a task published results show it can't do. (M)
