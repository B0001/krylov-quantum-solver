# Task Breakdown 3 — ShallowForge M1–M2: Metric Harness + First Rungs
Goal: CX@ε measurement infrastructure, frozen baselines, then R1 (grouping/ordering) and R2 (Trotter/randomization) with go/no-go readouts.

1. **Term-stream IR** — canonical representation of a Hamiltonian as ordered Pauli terms; SHA-256 of canonical form = `hamiltonian_hash`. Transforms are pure functions IR→IR.
   ✓ Hash invariant under term reordering of the *input*; round-trip serialization. (M)
2. **Manifest emitter** — every transform appends its entry per `compiler-manifest.schema.json` (params, predicted ε, lossless flag); totals assembled at export.
   ✓ Schema-validates; a lossless transform emits ε=0. (S)
3. **QASM3 export + verifier** — compile IR → circuit; verify on T0–T2 via statevector against exact `exp(-iHt)`.
   ✓ Fidelity ≥ 1−1e-8 per step at reference step size, all three systems. (M)
4. **CX@ε harness** — end-to-end: compiled evolution → ODMD → final energy error; binary-search step-count/params to the cheapest config with error ≤1.6 mHa on the golden suite. THIS is the only reporting path — no raw counts leave without ε.
   ✓ Harness reproduces the ~6,500 CX/step N₂ baseline figure; baselines frozen to `baselines.json`. (L)
5. **R1a Pauli grouping** — commuting-set grouping (greedy coloring v1).
   ✓ CX@ε improves on ≥3/4 golden systems → stays in default stack; multiplier recorded. (M)
6. **R1b ordering + cancellation** — magnitude vs lexicographic ordering; adjacent-exponential gate cancellation pass.
   ✓ Same go/no-go readout; interaction with R1a measured (stacked run). (M)
7. **R2a second-order Suzuki** — implement; predicted-ε model for Trotter order must bound observed ε.
   ✓ Prediction bounds observation on all four systems, or error model fixed before proceeding. (M)
8. **R2b randomized compiling (priority experiment)** — qDRIFT-style sampling; measure whether ODMD's noise tolerance absorbs the stochastic error more cheaply than deterministic formulas.
   ✓ Report CX@ε for randomized vs deterministic at matched ε — win or lose, this readout is publishable. (L)
9. **Interim report** — multiplier table (individual + best stacked), N₂ before/after, frontier plot v0.
   ✓ Every number carries its ε; manifests attached for every claimed circuit. (S)
