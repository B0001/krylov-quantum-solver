# Spec: Shallow-Evolution Circuit Compiler for Chemistry Hamiltonians

**Working name:** ShallowForge (rename freely)
**Status:** Draft v0.1
**Depends on:** krylov-quantum-solver (`qubitization_blueprint.py`, `thc_factorization.py`, `benchmark_resources.py`, ODMD pipeline as the downstream consumer)

---

## 1. Goal

Cut the two-qubit gate count per time-evolution step for molecular Hamiltonians by a **target factor of 5–10×** relative to the repo's current baseline (first-order Trotter, naive Pauli grouping — the ~6,500 CX/step N₂ figure), while keeping the evolution accurate enough that the downstream ODMD/Krylov energy estimate still reaches chemical accuracy.

This is the research-lever project: unlike the other three specs, success is not guaranteed. The spec is therefore structured as a sequence of measurable experiments, each independently publishable, rather than one product.

## 2. What counts as the invention

Any of: (a) a compilation pipeline (technique stack + ordering heuristics) achieving ≥5× CX reduction at fixed downstream energy error, (b) a new combination of THC factorization + qubitization tailored to Krylov/ODMD workloads specifically, or (c) a negative-but-rigorous result mapping the true Pareto frontier of depth vs. accuracy for these systems (valuable to the field, publishable, feeds the benchmark project's resource model).

## 3. The metric that matters (define once, use everywhere)

**CX@ε:** two-qubit gate count per evolution step such that the final ODMD ground-state energy error stays ≤ ε = 1.6 mHa on the golden suite. Comparing raw gate counts without fixing downstream ε is how compiler papers mislead; this project never reports a gate count without its ε.

Secondary metrics: total circuit depth (accounts for parallelism), ancilla count, classical compile time.

## 4. Technique ladder (implement and measure in this order)

Each rung is a self-contained experiment with a go/no-go readout. Gains multiply only if errors compose benignly — measure the stack, not just the rungs.

1. **R1 — Term ordering & Pauli grouping (cheap wins first).** Commuting-set grouping, lexicographic vs. magnitude ordering, gate cancellation across adjacent exponentials. Literature suggests 1.5–3× is available here. Baseline infrastructure for everything after.
2. **R2 — Higher-order / randomized Trotter.** Second-order Suzuki, randomized compiling (qDRIFT-style) tuned to ODMD's noise tolerance — ODMD's robustness to stochastic error is a genuine structural advantage: it may tolerate cheap randomized formulas that a plain phase-estimation pipeline cannot. This interaction is the most original angle in the project.
3. **R3 — THC-compressed Hamiltonians.** Use `thc_factorization.py` to reduce the term count before compilation; measure THC rank vs. downstream ε tradeoff explicitly (rank sweep, not a single point).
4. **R4 — Qubitization / block-encoding path.** Flesh out `qubitization_blueprint.py` into a counted, simulable circuit; compare CX@ε against the best Trotter stack. Qubitization trades depth for ancillas — report both.
5. **R5 — Numerical circuit optimization.** Peephole + template matching + (if time) small-block resynthesis on the output of the best stack. Bounded expectations: 10–30%.

## 5. Architecture

- Input: active-space Hamiltonian (same objects the solver already produces)
- IR: a term-stream representation that every rung transforms; each transform is a pure function with a recorded provenance tag (so any compiled circuit can be audited — this also satisfies the benchmark project's compilation-audit requirement)
- Output: circuit in QASM3 + a manifest (technique stack, parameters, predicted ε contribution per stage)
- Verification harness: statevector simulation on T0–T2 systems comparing compiled evolution against exact `exp(-iHt)`; operator-norm / spectral checks where statevector is too big

## 6. Experimental protocol (per rung)

1. Compile golden-suite Hamiltonians (H₂, H₄, LiH, N₂ CAS) with the rung enabled vs. disabled.
2. Verify unitary fidelity on simulable sizes; record CX count, depth, ancillas.
3. Run the **full ODMD pipeline** on the compiled circuits (noiseless sim) → measure final energy error.
4. Report CX@1.6mHa and the rung's multiplier. Go/no-go: rung stays in the default stack only if it improves CX@ε on ≥3 of 4 golden systems.

## 7. Validation plan

- **Correctness gate:** every compiled circuit's action verified against exact evolution on T0–T2 (fidelity ≥ 1 − 10⁻⁸ per step at the chosen step size) before any gate-count claim is made.
- **No-silent-approximation rule:** every accuracy-losing transform (THC rank truncation, Trotter order, randomization) must emit its predicted ε contribution into the manifest; total predicted ε must bound observed ε on the golden suite. If prediction and observation disagree, the error model is wrong — fix it before proceeding.
- **End-to-end proof point:** reproduce the N₂ CAS(6,6) energy to chemical accuracy through the full compiled pipeline, and state the new CX/step figure next to the ~6,500 baseline. That single before/after number is the headline result.

## 8. Risks

- **Gains don't stack.** R1–R3 may partially overlap (grouping gains shrink after THC compression). Mitigation: the protocol measures the stack combinatorially on at least the best-2 combinations, not just individually.
- **Qubitization ancilla cost eats the win** on near-term devices. Report honestly; it may still win in the fault-tolerant cost model — keep both cost models in the report.
- **5–10× is not reachable.** Then outcome (c) from §2 — a rigorous frontier map — is the deliverable, and it directly upgrades ChemCheck's Mode A resource model. No wasted work.
- **Scope creep into general-purpose compiler.** Hard scope fence: chemistry Hamiltonians consumed by ODMD/Krylov workloads only. General compilation is a different (crowded) field.

## 9. Milestones

1. **M1 — Metric + harness (1–2 weeks):** CX@ε measurement pipeline, baseline numbers frozen for all four golden systems.
2. **M2 — R1+R2 (2–3 weeks):** cheap wins + Trotter/randomization study; first multiplier report. The R2×ODMD-robustness interaction is the priority experiment.
3. **M3 — R3+R4 (3–4 weeks):** THC rank sweep; qubitization circuit counted and simulated.
4. **M4 — Stacked result + writeup:** best stack, N₂ before/after headline, frontier plot. Decide: paper, or fold into the solver as the default compilation path (ideally both).

## 10. Success criteria

- Primary: ≥5× reduction in CX@1.6mHa on N₂ CAS(6,6) vs. frozen baseline, verified end-to-end through ODMD.
- Fallback (still success): a published Pareto frontier with a validated error model, adopted as ChemCheck's resource model.
- All compiled circuits carry audit manifests; zero correctness-gate failures in the final report.
