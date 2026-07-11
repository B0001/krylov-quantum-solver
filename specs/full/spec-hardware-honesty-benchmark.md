# Spec: Quantum Hardware Chemistry Honesty Benchmark

**Working name:** ChemCheck (rename freely)
**Status:** Draft v0.1
**Depends on:** krylov-quantum-solver (`benchmark_resources.py`, variational floor, golden reference suite, ODMD noise-robustness machinery)

---

## 1. Goal

A standardized, reproducible scorecard answering one question: **"Can this quantum computer do real chemistry yet — and if not, how far off is it?"** Vendors report gate fidelities and qubit counts; neither tells a chemist whether a device can produce a usable N₂ dissociation curve. This benchmark converts hardware specs and (optionally) real device runs into a single honest verdict per molecule tier.

## 2. What counts as the invention

The invention is the **scoring methodology**: a fixed suite of chemistry problems with exact references, a resource model that maps circuit requirements onto device error budgets, and a hard "nonsense detector" (variational floor + bracket containment) that catches results which look plausible but are physically impossible. Comparable in spirit to MLPerf, but for quantum chemistry — the value is the neutral, gameable-resistant protocol.

## 3. Benchmark tiers

| Tier | System | Why it's there |
|---|---|---|
| T0 | H₂ / STO-3G (4 spin orbitals) | Sanity floor — any credible device or simulator must pass |
| T1 | H₄ chain (8 spin orbitals) | First multireference stress; correlation matters |
| T2 | LiH (validated CAS) | Realistic small molecule; ionic + covalent character |
| T3 | N₂ stretched, CAS(6,6) | Strongly multireference; HF fails by ~0.5 Ha; ~6,500 two-qubit gates per Trotter step — the current hardware wall |
| T4 (aspirational) | Nb₃X₈-class cluster active space | Materials-relevant target; defines "chemistry-useful" |

Each tier ships with: Hamiltonian spec, exact FCI reference energy, required accuracy (chemical accuracy: 1.6 mHa), and the resource count from `benchmark_resources.py` (qubits, two-qubit gates per Trotter step, total depth for the reference protocol).

## 4. Two evaluation modes

### Mode A — Paper score (no hardware access needed)
Input: a device spec sheet (qubit count, two-qubit gate error, T₁/T₂, connectivity, native gate set).
Computation: map each tier's circuit onto the device (routing overhead included — all-to-all vs heavy-hex matters enormously), compute expected total error via a standard error-budget model, output pass/fail per tier with a **fidelity headroom number**: "you need two-qubit error ≤ X to pass T3; you are at Y — a factor of Z away."
This mode makes the benchmark useful immediately, before anyone lends you a QPU.

### Mode B — Live score (device runs)
Input: measurement data from actually executing the tier circuits (or a vendor-submitted results file in the benchmark's schema).
Scoring per tier:
1. Energy estimate extracted via the ODMD pipeline (noise-robust by design — this is a fair-to-hardware choice, not a gotcha).
2. **Floor check:** any energy below the variational floor = automatic FAIL flagged `UNPHYSICAL`, regardless of how close it looks to the reference. This is the anti-fraud core.
3. **Accuracy check:** |E − FCI| ≤ 1.6 mHa → PASS; ≤ 16 mHa → MARGINAL; else FAIL.
4. **Reproducibility check:** ≥3 independent runs; spread must be consistent with reported shot noise.

## 5. Anti-gaming provisions

- Circuits are fixed and published; compilation is allowed but the compiled circuit must be submitted and is re-verified against the Hamiltonian (no "optimizing" to a different, easier problem).
- Error mitigation is allowed but must be declared; mitigated and raw scores are both reported.
- Classical shadow / post-processing budgets are capped and disclosed (otherwise a classical computer does the chemistry and the QPU is decoration).
- A "classically spoofable" flag: T0–T2 are trivially classically simulable, and the scorecard says so — passing them proves correctness of the stack, not quantum advantage. Only honest framing survives scrutiny.

## 6. Outputs

- Per-device scorecard: table of tiers × (paper score, live score, headroom factor), one summary verdict line
- Public leaderboard-ready JSON schema + a rendered Markdown/HTML report
- Reproducibility bundle: exact Hamiltonians, circuits, seeds, and scoring code pinned by version hash

## 7. Validation plan

1. **Self-test:** run Mode B against the repo's own noiseless simulator — every tier must PASS. Then against simulators with injected depolarizing noise at known rates — Mode A's predictions must match Mode B's observed pass/fail transitions within stated tolerance. This calibrates the paper-score model against ground truth.
2. **Floor-detector test:** feed the scorer the old codebase's known-bad outputs (energies far below the floor); it must flag every one UNPHYSICAL.
3. **Compilation audit test:** submit a deliberately altered circuit; the Hamiltonian re-verification must reject it.
4. **External sanity:** score 2–3 real devices from published spec sheets (Mode A) and check the verdicts against published experimental chemistry results on those devices — the benchmark must not certify a device that published results show can't do the task, and vice versa.

## 8. Risks

- **Vendor pushback / fairness disputes.** Mitigation: publish the full methodology, accept methodology PRs, version the benchmark (ChemCheck-2026.1) so scores are comparable within a version.
- **Error-budget model too crude.** Depolarizing-only models mispredict real devices. Start simple, state assumptions, refine with Mode B calibration data; report model uncertainty on headroom factors.
- **Benchmark rot.** Hardware improves; tiers must be extensible upward (T5+) without invalidating old scores — hence versioning.
- **Being ignored.** Distribution matters as much as correctness: the launch artifact should be a readable report scoring today's flagship devices, not a bare repo.

## 9. Milestones

1. **M1 — Tier definitions frozen (1 week):** Hamiltonians, references, resource counts, accuracy thresholds published in one document.
2. **M2 — Mode A scorer (2 weeks):** spec-sheet → scorecard, calibrated against noisy simulation.
3. **M3 — Mode B scorer + anti-gaming checks (2–3 weeks):** results-file schema, floor detector, compilation audit.
4. **M4 — Launch report:** score 3–5 current devices in Mode A; publish "State of Quantum Chemistry Hardware 2026" with the headroom factors as the headline numbers.

## 10. Success criteria

- Mode A pass/fail predictions match noisy-simulation ground truth across the calibration sweep.
- Floor detector catches 100% of injected unphysical results, zero false positives on the golden suite.
- The launch report produces one memorable, defensible number per flagship device ("factor of ~40 in two-qubit error from passing T3") that others cite.
