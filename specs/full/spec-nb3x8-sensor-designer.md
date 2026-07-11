# Spec: Nb₃X₈ Quantum Sensor Design Tool

**Working name:** SenseForge (rename freely)
**Status:** Draft v0.1
**Depends on:** krylov-quantum-solver (ODMD suite: `nb3x8_strain.py`, `nb3x8_magnetometry.py`, `nb3x8_susceptibility.py`, certified-bracket modules)

---

## 1. Goal

A computational design pipeline that screens Nb₃X₈ cluster variants (X = Cl, Br, I) across strain and magnetic-field regimes to identify the composition and operating point with the **highest spin-gap sensitivity**, i.e. the best candidate for a strain gauge or magnetometer based on a 2D magnet.

The deliverable of a pipeline run is not an energy — it is a **ranked list of sensor design candidates**, each with a figure of merit and a certified error bracket.

## 2. What counts as the invention

The invention is the *design output*, not the solver: "A strain/magnetic-field sensor comprising an Nb₃X₈ layer with halide X and pre-strain ε₀, operated at the point of maximal |d(gap)/dε| (or |d(gap)/dB|), as identified by certified quantum-chemical screening."

A defensible result = a specific (X, ε₀, B₀) operating point plus predicted sensitivity with rigorous bounds.

## 3. Inputs

| Input | Format | Notes |
|---|---|---|
| Cluster geometry | CIF → finite dimer cluster (existing path) | Document cluster size used |
| Halide choice | Cl / Br / I | Iodide flagged: spin–orbit coupling caveat (see §8) |
| Strain sweep | ε ∈ [−2%, +2%], step 0.25% (default) | Uniaxial and biaxial modes |
| Field sweep | B ∈ [0, 10] T, step 0.5 T (default) | Zeeman term added to Hamiltonian |
| Active space | (n_elec, n_orb) per system | Keep ≤ existing validated sizes |
| Method config | ODMD settings, Krylov depth, shots (if sampled) | Reuse existing defaults |

## 4. Core computation (per grid point)

1. Build clamped/strained cluster geometry; generate active-space Hamiltonian.
2. Run ODMD spin-gap extraction (singlet–triplet or relevant gap).
3. Attach **certified two-sided bracket** to the gap (`certified_gaps.py` / Temple-style bounds).
4. Store: gap, bracket width, convergence diagnostics, variational-floor sanity check (reject any energy below the floor — hard fail, never silently keep).

## 5. Derived quantities

- **Strain sensitivity:** S_ε = d(gap)/dε via central finite differences over the sweep, with error propagated from bracket widths.
- **Field sensitivity:** S_B = d(gap)/dB likewise.
- **Figure of merit (FoM):** |S| / (bracket width at that point). A huge slope with a huge error bar is worthless; this ratio rewards *credible* sensitivity.
- **Operating-point stability:** second derivative — prefer plateaus of high |S| over knife-edge points.

## 6. Outputs

- `results/{X}/gap_vs_strain.csv`, `gap_vs_field.csv` (gap, lower, upper per point)
- Sensitivity curves (matplotlib PNG + CSV)
- `candidates.md`: ranked table — halide, operating point, S_ε, S_B, FoM, bracket width, convergence flags
- A per-candidate one-page "design card" suitable for a paper figure or provisional patent appendix

## 7. Validation plan (gates, in order)

1. **Regression gate:** existing test suite green; H₂/H₄/LiH/N₂ references still reproduce to documented tolerance.
2. **Cluster-size convergence:** repeat the top candidate at ≥2 cluster sizes; the qualitative ranking must survive. If it doesn't, the result is a cluster artifact — report, don't publish.
3. **Cross-method check:** top 1–2 candidates re-computed with DMRG (or AFQMC) on the same active space. Agreement within combined error bars required.
4. **Bracket honesty audit:** for every system with an exact FCI reference available, the certified bracket must contain the FCI value 100% of the time. One violation = stop and debug the bounds, nothing else matters until fixed.

## 8. Known risks and honest caveats

- **Finite cluster ≠ periodic solid.** All claims must be phrased as cluster-model predictions until periodic validation exists. This goes in every output file header automatically.
- **Iodide needs spin–orbit coupling.** Without SOC, Nb₃I₈ numbers are indicative only. Either add SOC (large scope increase) or restrict launch claims to Cl/Br.
- **Active-space sensitivity.** Run the top candidate at two active-space sizes; report both.
- **Sensitivity ≠ device.** Readout mechanism (optical? transport?) is out of scope for v1; the deliverable is the material/operating-point design.

## 9. Milestones

1. **M1 — Sweep harness (1–2 weeks):** driver script that runs strain/field grids over existing modules, writes CSVs. No new physics.
2. **M2 — Certified derivatives (1 week):** finite-difference sensitivities with propagated brackets; FoM ranking.
3. **M3 — Validation pass (2–3 weeks):** cluster-size and DMRG cross-checks on top candidates.
4. **M4 — Design cards + writeup:** candidate report; decide paper vs provisional patent filing.

## 10. Success criteria

- At least one (X, ε₀ or B₀) candidate whose FoM survives M3 validation.
- Predicted gap shift measurable at experimentally reachable strain (|ε| ≤ 2%) or field (≤ 10 T).
- Full pipeline reproducible from one config file + one command.
