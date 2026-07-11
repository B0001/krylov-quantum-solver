# Task Breakdown 11 — #8 Co-Design Advisor (SEQUENCED: after ChemCheck task 6 calibration)
Goal: molecule + CAS in → device-requirement table out.

1. **Requirement inversion API** — reuse ChemCheck's error-budget model inverted (headroom machinery already computes required error); generalize from fixed tiers to arbitrary cap-compliant molecules via `benchmark_resources.py` counts.
   ✓ For the four golden systems, output matches ChemCheck tier thresholds exactly (consistency check). (M)
2. **Connectivity scenarios** — requirement table rows per connectivity class; columns: qubits, required 2q error, depth, with/without-ShallowForge (consume `baselines.json` + best-stack manifests).
   ✓ ShallowForge column populated from real manifests, never hand-entered. (S)
3. **Vendor-distance view** — join against the ChemCheck device registry: factor-to-requirement per known device.
   ✓ Reuses scorecard data; no new claims beyond calibrated model. (S)
4. **Uncertainty honesty** — every requirement carries the Mode-A model uncertainty band; CLI refuses to run if calibration data older than the model version.
   ✓ Stale-calibration refusal tested. (S)
5. **CLI + report** — `codesign advise mol.xyz --cas 6,6` → Markdown table.
   ✓ Golden-system outputs reviewed against ChemCheck launch report for contradiction — zero tolerated. (M)
