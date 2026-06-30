# SPEC: Classical shadows estimate the energy unbiasedly, with the shadow-norm variance bound

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); `classical_shadows.py` merged. Unbiased on HF and
FCI states (within 4·stderr), ~1/√shots convergence, empirical single-shot variance 1.95 ≤ shadow
norm 3.09. Finding: weight-≥3 terms carry ≈22% of the norm (3^w cost; max weight 4 on H₂).

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

Estimate energies and Pauli observables from randomized single-qubit (random-Pauli) measurements —
the classical-shadows protocol. Claim: the random-Pauli shadow energy estimator is **unbiased**
(its mean converges to the exact `⟨ψ|H|ψ⟩` for any state, HF or correlated FCI), its error shrinks
as `1/√shots`, and its single-shot variance is bounded by the **shadow norm** `Σ_k |c_k|² 3^{w_k}`.
The claim is false if the estimator is biased, if it fails to converge, or if the empirical variance
exceeds the shadow norm.

## 2. Background and honest framing

- **Prior art / reference.** Huang, Kueng & Preskill, *Predicting many properties of a quantum system
  from very few measurements*, Nature Physics 16, 1050 (2020), arXiv:2002.08953. A measurement /
  estimation primitive, distinct from the energy *methods* (Krylov, PDS, …) in this repo; it
  connects to the RDM-measurement theme of the QKSD-properties and PDS rungs.
- **Ground truth.** The exact expectation `⟨ψ|H|ψ⟩` (dense) for the HF and FCI states of the same
  qubit Hamiltonian.
- **What we can claim if gates pass.** Random-Pauli classical shadows give an unbiased, `1/√shots`
  energy estimator with the HKP shadow-norm variance bound, validated on HF and FCI states.
- **What we cannot claim (stated up front).** (a) No quantum advantage — exact statevector
  simulation of the measurement, tiny systems. (b) **The 3^{w} cost is the honest finding:** the
  shadow norm grows as `3^{weight}`, so high-weight Pauli terms are sample-expensive — random-Pauli
  shadows are not the cheapest route for molecular Hamiltonians; grouped / derandomized shadows
  mitigate this and are **out of scope**. (c) Only the energy / Pauli expectations here, not the full
  many-property prediction or fidelity estimation.

## 3. Approach

For each snapshot: measure every qubit in a uniformly random X/Y/Z basis (rotate the statevector,
sample a Z-basis bitstring). The unbiased single-qubit estimator of a Pauli `P_q` is `3 s_q` when the
measured basis matches and 0 otherwise; a Pauli string's estimate is the product over its support;
the energy estimate is `Σ_k c_k ·` (term estimate). Average over snapshots. Reference: exact
`⟨ψ|H|ψ⟩`; variance bound: `shadow_norm = Σ_k |c_k|² 3^{w_k}`.

## 4. Public interface

```
classical_shadows.collect_classical_shadow(statevector, n_qubits, n_shots, seed=0) -> (bases, signs)
classical_shadows.shadow_energy_samples(bases, signs, pauli_op) -> np.ndarray   # per-shot; mean = <H>
classical_shadows.shadow_norm(pauli_op) -> float                                 # HKP variance bound
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_classical_shadows_spec.py` (test-first). Seeded (deterministic); exact
statevector; pyscf/qiskit, no block2.

- **G1 — unbiased on HF and FCI states.** The shadow energy mean is within `4·stderr` of the exact
  `⟨ψ|H|ψ⟩` for both the HF reference and the FCI ground statevector of H₂, at high shots.
- **G2 — `1/√shots` convergence (definition of done).** The |estimate − exact| error at a high shot
  budget is smaller than at a low budget (seeded means), consistent with `~1/√shots`.
- **G3 — shadow-norm variance bound.** The empirical single-shot variance of the energy estimator is
  `≤ shadow_norm` (the HKP bound) on H₂, and `shadow_norm` is finite.
- **G4 — the 3^{w} cost (the finding).** The shadow norm is materially inflated by high-weight Pauli
  terms: `shadow_norm` is strictly larger than a weight-capped reference `Σ_k |c_k|² 3^{min(w_k,1)}`,
  and the weight-≥3 terms carry a non-trivial share (measured ≈ 22% on H₂, max weight 4) — random-
  Pauli shadows pay `3^{weight}` per term.

> Definition of done: **G1 + G2**. If the estimator is found biased (e.g. a qubit-ordering or
> single-qubit-estimator error), it fails G1 loudly — that is the test's main job. The 3^{w} finding
> (G4) is the honest cost caveat, not a defect.

## 6. Implementation plan (test-first)

1. Write `tests/test_classical_shadows_spec.py` encoding G1–G4 (initially failing — module absent).
2. Add `classical_shadows.py` (random-Pauli snapshots; per-shot estimator; shadow norm).
3. Iterate to green via `make gates` (own process; pyscf/qiskit, no block2).

## 7. Out of scope

- Grouped / derandomized / locally-biased shadows and Clifford (global) shadows (the cost reductions).
- RDMs, fidelities, and the full many-property prediction; observable batching.
- Shot-level hardware noise on the measurements.

## 8. Caveats and risks

- **R1 — RNG flakiness.** Estimates are random. *Mitigation:* seed the RNG and gate within `4·stderr`
  (a >99.99% band), not a hand-tuned absolute.
- **R2 — measurement-simulation cost.** Each snapshot rotates and samples the full statevector.
  *Mitigation:* small systems (H₂, n=4) keep the gate fast.
- Honest limitation: exact statevector, minimal-basis tiny molecules; random-Pauli shadows pay
  `3^{weight}` — a correctness + variance study, not an efficient-measurement claim.

## 9. Deliverables

- `classical_shadows.py` — `collect_classical_shadow`, `shadow_energy_samples`, `shadow_norm`.
- `tests/test_classical_shadows_spec.py` — gates G1–G4.
- Results summary (unbiasedness + `1/√shots` + the shadow-norm bound and 3^{w} cost, with §2/§7
  caveats) in the PR description.
