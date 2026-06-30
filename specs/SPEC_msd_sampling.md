# SPEC: Mirror subspace diagonalization lowers the quantum-Krylov sampling cost

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); `msd.py` merged. Measured: N₂ CAS(6,6)
order-8 + energy shift gives fd1=5.48 < λ=14.75, MSD median error ≈ 3.2× below KQD at 10⁵ shots;
H₂ (λ/W≈1.3) gives fd1>λ and no advantage (the boundary). Modest vs the paper's 10–10⁴× — that
needs larger λ/W than dense-diagonalizable systems reach.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

Mirror subspace diagonalization (MSD) estimates the projected Hamiltonian matrix element
`H_{0k} = ⟨φ₀|H e^{−iHkτ}|φ₀⟩` as `i·S′(kτ)` — a central finite-difference of the overlap function
`S(t) = ⟨φ₀|e^{−iHt}|φ₀⟩` at symmetrically shifted timesteps — instead of an LCU/Hadamard-test sum
over the Pauli terms of H. Its sampling variance scales with the stencil 1-norm `fd1 = ‖w‖₁/δ`
rather than the Hamiltonian 1-norm λ. Claim: with an energy-level shift (centering the reachable
spectrum so the difference need only resolve the spectral width W) and a high-order stencil
(`fd1 ∝ ‖w‖₁/bias^{1/order}`), `fd1 < λ` is achievable at chemical-accuracy bias, and MSD then
estimates the ground-state energy to the same accuracy as standard KQD with **fewer shots**. The
claim is false if the finite-difference construction does not reproduce KQD/FCI, or if MSD does not
beat KQD at matched shots in the regime where `fd1 < λ`.

## 2. Background and honest framing

- **Prior art / reference.** Mirror subspace diagonalization, arXiv:2511.20998; KQD sampling-error
  analysis (Lee-Lee-Huh, Quantum 8, 1477; Kirby, Quantum 8, 1457; Epperly-Lin-Nakatsukasa, 2022).
  Builds on the validated real-time Krylov primitive and the noise machinery here.
- **What we can claim if gates pass.** The MSD finite-difference + energy-shift construction is
  correct (reproduces KQD/FCI), and on N₂ CAS(6,6) MSD gives a real, measured sampling-cost
  reduction (≈ 3× lower median error at 10⁵ shots ⇒ several-fold fewer shots for a target accuracy).
- **What we cannot claim (stated up front).** (a) **The advantage is a `λ/W` effect and is modest at
  this scale.** Minimal-basis validatable systems have `λ/W = O(1–3)`, so the realized advantage is
  ~2–3×, *not* the paper's 10–10⁴× (which needs much larger `λ/W` than dense-diagonalizable systems
  provide). (b) On H₂ (`λ/W ≈ 1.3`) there is **no advantage** — `fd1 > λ` — the recorded boundary.
  (c) Idealized i.i.d. Hadamard-test shot noise on the matrix elements (the paper's variance model),
  not a full circuit/Trotter simulation. (d) Exact statevector overlaps; a correctness + sampling
  study, not a hardware result.

## 3. Approach

Build the exact KQD Toeplitz elements (`S_{0k}`, `H_{0k}`) and the MSD elements (`i·S′(kτ)` by a
central stencil of overlaps at `kτ ± jδ`) on the energy-shifted Hamiltonian `H − μI`, μ centering
the HF-reachable spectrum. Model per-element shot noise with the paper's variance: overlaps
`∝ 2(2−1/d)/m`, Hamiltonian elements `∝ scale²(2−1/d)/m` with `scale = λ` (KQD) or `fd1` (MSD).
Solve the Hermitian-Toeplitz GEVP with a noise-aware overlap floor; compare the **median** ground-
energy error (robust to GEVP blow-ups) over many seeds. Reference: the noiseless KQD ground value
(= FCI up to Krylov order) and `MolecularHamiltonian` dense diagonalization.

## 4. Public interface

```
msd.central_difference_weights(order)                    -> (offsets, weights)   # order in {2,4,6,8}
msd.build_msd_problem(mh, n=8, order=8, delta=None, bias_target=1.6e-3) -> MSDProblem
msd.sample_ground_energy(prob, shots, method, seed)      -> float                # method in kqd/msd
msd.rms_error(prob, shots, method, seeds)                -> float
MSDProblem: .s .h_kqd .h_msd .lam .fd1 .width .ref .offset .msd_bias
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_msd_sampling_spec.py` (test-first). Median over ≥ 120 seeds. PySCF/qiskit, no
block2. (N₂ CAS(6,6) is dense-diagonalized — this gate is on the slower side.)

- **G1 — construction correct.** On N₂ CAS(6,6) with the energy shift + order-8 stencil, the
  finite-difference bias `|MSD − KQD|_noiseless < 1.6 mHa`, and the noiseless KQD ground equals FCI
  to `< 2 mHa` (Krylov order n=8).
- **G2 — sampling advantage (definition of done).** On N₂ CAS(6,6) at 10⁵ shots, the MSD median
  ground-energy error is below KQD's, by `> 1.5×` (measured ≈ 3.2×).
- **G3 — boundary: no advantage at small λ/W.** On H₂ (`λ/W ≈ 1.3`), `fd1 > λ` and MSD's median
  error is *worse* than KQD's — the honest regime boundary (advantage requires `λ/W` large enough
  that a high-order stencil drives `fd1 < λ`).
- **G4 — mechanism.** A higher-order stencil holds the bias at a larger δ and so shrinks `fd1`:
  `fd1(order 8) < fd1(order 4)` at matched (`< 1.6 mHa`) bias on N₂, and order-8 reaches `fd1 < λ`
  (the advantage regime).

> Definition of done: **G2 + G3** together — a measured advantage *and* the boundary that says where
> it vanishes. If on a larger system `λ/W` grows and the advantage exceeds the modest value here,
> record the scaling (it is the paper's point); if a system in the `fd1 < λ` regime shows *no*
> advantage, that falsifies the mechanism and must be recorded.

## 6. Implementation plan (test-first)

1. Write `tests/test_msd_sampling_spec.py` encoding G1–G4 (initially failing — `msd` absent).
2. Add `msd.py` (energy shift, central-difference H, per-element shot-noise model), reusing
   `solve_generalized_eig`.
3. Iterate to green via `make gates` (own process; pyscf/qiskit, no block2).

## 7. Out of scope

- A full circuit-level / Trotterized Hadamard-test simulation and hardware noise (the i.i.d.
  matrix-element variance model is the paper's analysis tool; we use it directly).
- The Moment-Lanczos error-mitigation extension and resource estimation of arXiv:2511.20998.
- Demonstrating the large `λ/W` (10–10⁴×) regime — needs systems beyond dense diagonalization.

## 8. Caveats and risks

- **R1 — regularization-policy sensitivity.** The advantage magnitude depends on the GEVP overlap
  floor; RMS is outlier-sensitive when KQD blows up. *Mitigation:* a common overlap-noise floor for
  both methods (they measure S identically) and the **median** error metric; gate a conservative
  `> 1.5×`, not the headline value.
- **R2 — finite-difference bias floor.** At high shots MSD becomes bias-limited; the advantage is in
  the variance-limited regime. *Mitigation:* gate at 10⁵ shots and keep bias `< 1.6 mHa`.
- Honest limitation: idealized shot noise, minimal-basis tiny systems, modest advantage — a
  mechanism + scaling-boundary demonstration, not a hardware speedup claim.

## 9. Deliverables

- `msd.py` — `central_difference_weights`, `build_msd_problem`, `sample_ground_energy`, `rms_error`,
  `MSDProblem`.
- `tests/test_msd_sampling_spec.py` — gates G1–G4.
- Results summary (the measured advantage, the H₂ boundary, the `λ/W` scaling, with §2/§7 caveats)
  in the PR description.
