# SPEC: Hamiltonian-moment energies — PDS is a variational upper bound that converges to FCI

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); `moment_expansion.py` merged. PDS(K) variational
(≥ FCI) at every K and converges (H₄ 67→10→2.0→0.43 mHa over K=1..4, LiH 0.68 mHa at K=4, H₂ exact
at K=2); PDS(1)=⟨H⟩. CMX(2) dips below FCI on H₂ (−0.27 mHa) — non-variational, the boundary.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

Estimate the ground-state energy from the Hamiltonian *moments* `μ_n = ⟨φ|H^n|φ⟩` of the
Hartree–Fock reference — no time evolution, no Hilbert-space eigenproblem. Claim: the PDS(K)
functional (smallest root of `P_K(E)` built from the moment matrix) is a **variational upper bound**
that converges to FCI as K grows, while the connected-moment expansion CMX(2) — the same moment data,
resummed differently — is **not** variational and can dip *below* FCI. The claim is false if PDS ever
drops below FCI, if PDS(K) does not reach chemical accuracy at attainable K, or if CMX is found to be
variational here.

## 2. Background and honest framing

- **Prior art / reference.** Peng & Kowalski, *Variational quantum solver employing the PDS energy
  functional*, Quantum 5, 473 (2021), arXiv:2101.08526 (PDS); Kowalski & Peng, J. Chem. Phys. 153,
  201102 (2020) (quantum CMX); Peeters–Devreese–Soldatov original. A distinct method family from the
  Krylov/subspace rungs here.
- **Ground truth.** FCI = dense diagonalization of the same qubit Hamiltonian (`ground_state_energy`).
- **What we can claim if gates pass.** PDS(K) from the HF reference is a variational upper bound that
  reaches chemical accuracy on H₂/H₄/LiH by K=4, tightening monotonically with K; and PDS's
  variational guarantee is a genuine advantage over CMX, which is shown non-variational.
- **What we cannot claim.** (a) No quantum advantage — exact statevector moments, tiny systems; a
  correctness study (on hardware the cost is measuring `⟨H^n⟩`, out of scope). (b) **Reference-quality
  limited:** convergence rate depends on the HF reference overlap; a poor reference needs higher K.
  (c) **High-K ill-conditioning:** the moment matrix M becomes singular as `|φ⟩ → exact` and at large
  K (high powers) — so K is kept low (≤ ~4–5); this is the paper's own caveat, recorded not hidden.

## 3. Approach

Compute `μ_0..μ_{2K-1}` by repeated sparse mat-vec on the HF statevector. PDS(K): `M_ij =
⟨H^{2K-i-j}⟩`, `Y_i = ⟨H^{2K-i}⟩`, solve `M X = -Y`, energy = smallest root of
`P_K(E) = E^K + Σ X_i E^{K-i}` (Peng–Kowalski Eqs. 11–13). CMX(2): connected moments
`I_1=μ_1, I_2=μ_2-μ_1², I_3=μ_3-3μ_1μ_2+2μ_1³`, `E = I_1 - I_2²/I_3`. Reference: FCI.

## 4. Public interface

```
moment_expansion.hamiltonian_moments(mh, max_order) -> (moments: np.ndarray, offset: float)
moment_expansion.pds_energy(moments, order, offset=0.0) -> float   # variational upper bound
moment_expansion.cmx2_energy(moments, offset=0.0) -> float          # non-variational contrast
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_moment_pds_spec.py` (test-first). Exact statevector; pyscf/qiskit, no block2.

- **G1 — PDS(1) = ⟨H⟩.** PDS(1) equals the HF expectation value `μ_1 + offset` (= the Rayleigh
  quotient, the RHF energy) to `< 1e-9` on H₂/H₄/LiH.
- **G2 — variational + convergence (definition of done).** For H₂, H₄, LiH: PDS(K) `≥ E_FCI - 1e-9`
  at every K ∈ {1,2,3,4} (variational), and PDS(4) is within chemical accuracy `< 1.6 mHa` of FCI
  (measured: H₂ exact at K=2, H₄ 0.43, LiH 0.68 mHa at K=4).
- **G3 — monotone tightening.** The upper bound improves with order: `PDS(K+1) ≤ PDS(K) + 1e-9` for
  K ∈ {1,2,3}.
- **G4 — CMX is not variational (the boundary).** On H₂, CMX(2) lies *below* FCI
  (`cmx2 < E_FCI - 1e-4`, measured −0.27 mHa) while PDS(2) stays `≥ E_FCI` — the same moments,
  resummed, lose the variational guarantee. (Records why PDS is the robust choice.)

> Definition of done: **G2**. If PDS does not reach chemical accuracy by K=4 on some system (poor HF
> overlap), raise K and record where the moment matrix becomes ill-conditioned — that is the finding.

## 6. Implementation plan (test-first)

1. Write `tests/test_moment_pds_spec.py` encoding G1–G4 (initially failing — `moment_expansion` absent).
2. Add `moment_expansion.py` (moments by sparse mat-vec; PDS linear solve + polynomial roots; CMX(2)).
3. Iterate to green via `make gates` (own process; pyscf/qiskit, no block2).

## 7. Out of scope

- The PDS(K)-VQS variational-circuit optimizer and energy gradients (the paper's main contribution);
  here PDS is a static functional on a fixed HF reference.
- Hardware moment measurement / Pauli grouping cost; shot noise on `⟨H^n⟩`.
- Excited states from the higher roots of `P_K(E)` (a natural follow-up).
- High-order PDS (K ≳ 5) where the moment matrix is severely ill-conditioned.

## 8. Caveats and risks

- **R1 — ill-conditioning at high K.** `det M → 0` as the reference approaches the exact state and at
  large K. *Mitigation:* gate K ≤ 4; record the conditioning limit rather than pushing K blindly.
- **R2 — reference dependence.** A poor HF overlap slows convergence. *Mitigation:* gate on
  near-equilibrium geometries with a decent HF reference; a stretched/multireference case would need
  higher K (a follow-up finding).
- Honest limitation: exact statevector moments, minimal-basis tiny molecules — a correctness study.

## 9. Deliverables

- `moment_expansion.py` — `hamiltonian_moments`, `pds_energy`, `cmx2_energy`.
- `tests/test_moment_pds_spec.py` — gates G1–G4.
- Results summary (PDS convergence + variational bound, the CMX non-variational contrast, §2/§7
  caveats) in the PR description.
