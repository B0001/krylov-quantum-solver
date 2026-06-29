# SPEC: Sample-based Krylov quantum diagonalization (SKQD) reproduces the exact-Krylov floor

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

Adding a **sample-based** Krylov path — diagonalize the Hamiltonian in the subspace of Slater
determinants *sampled* from real-time-evolved states |Ψ_k⟩ = e^(−ikΔtH)|Ψ_HF⟩, rather than from
the dense ⟨φ_i|H|φ_j⟩ matrix elements — recovers the same ground-state energy as the validated
exact-evolution `QuantumKrylovSolver`, to within chemical accuracy, while never dropping below the
FCI floor. Claim: on H₄ and N₂(CAS(6,6)), as the number of sampled determinants and the Krylov
depth grow, the SKQD energy converges to FCI from above.

## 2. Background and honest framing

SKQD (Yu et al. 2025, `arXiv:2501.09702`) and its quantum-chemistry-friendly qDRIFT variant
**SqDRIFT** (Piccinelli et al., `arXiv:2508.02578`) bridge the project's two existing near-term
pillars: real-time quantum Krylov (`quantum_krylov_solver.py`) and sample-based diagonalization
(`run_nbn_sqd_sweep.py`, `qiskit-addon-sqd`). Instead of forming Krylov matrix elements ⟨ψ_i|H|ψ_j⟩
(which off-diagonal need controlled time evolution / Hadamard tests — see `hardware_krylov.py`),
SKQD *samples* bitstrings from each |Ψ_k⟩, unions them into a determinant subspace, and diagonalizes
H classically there. The two papers prove convergence guarantees under the same assumptions as QPE,
**provided the ground state is "concentrated"** (supported on poly(n) determinants) and the
reference has non-negligible overlap with it.

- **What we can claim if the gates pass:** our SKQD implementation reproduces the SKQD method's
  central behavior on validated small systems — convergence to FCI from above as samples/depth grow —
  with the convergence checked against our own exact-Krylov solver and PySCF FCI. It connects the
  Krylov and SQD code paths under one tested interface.
- **What we cannot claim:** no novelty (this reproduces a published method), and **no quantum
  advantage at this scale** — sampling is simulated classically from the exact statevector, the
  active spaces are FCI-tractable, and "concentration" holds trivially for these small ground
  states. We are validating the algorithm, not demonstrating utility. We do not run on hardware and
  do not (in this spec) implement the qDRIFT circuit compilation — sampling is from the exact
  |Ψ_k⟩ amplitudes.

## 3. Approach

Reuse the exact real-time propagation already in `QuantumKrylovSolver._ensure_basis` to produce
|Ψ_k⟩ statevectors. For each k, sample `n_shots` computational-basis bitstrings from |⟨b|Ψ_k⟩|²
(restricted to the correct particle-number sector, as SQD does), union the determinants across all
k into a subspace D, build the Hamiltonian projected onto D, and diagonalize. Reference values: the
exact-evolution `QuantumKrylovSolver` energy at the same Δt and depth, and PySCF FCI (the variational
floor). Determinant ↔ qubit-bitstring convention follows `molecular_hamiltonian` / the SQD path.

## 4. Public interface

```
hybrid_quantum_solver.skqd.SampleKrylovSolver(molecular_hamiltonian, dt=None,
    n_shots=int, depth=int, seed=None)
    .solve() -> SKQDStep(depth, n_dets, energy)          # energy in Ha, incl. offset
    .convergence(shot_schedule) -> list[SKQDStep]
benchmark_skqd.py   -> CSV: system, depth, n_shots, n_dets, E_skqd, E_exactKrylov, E_fci, err_mHa
```

Compose, don't rewrite: propagation from `QuantumKrylovSolver`; the subspace eigensolve from the
SQD machinery already in the tree; FCI from `dmrg_reference.fci_energy`.

## 5. Acceptance criteria (validation gates)

Gates live in `tests/test_skqd_spec.py` (test-first), small/fast cases.

- **G1 — variational floor.** For H₄ (1.0 Å, 8 qubits) and N₂ CAS(6,6), at every depth and shot
  count tested, `E_skqd ≥ E_fci − 1e-6 Ha`. (A sampled-subspace Rayleigh quotient can never beat
  FCI; a violation means a determinant-mapping or projection bug.)
- **G2 — convergence to FCI.** For H₄, with a large shot budget (e.g. 5×10⁴) and depth ≥ 6,
  `|E_skqd − E_fci| < 1.6 mHa` (chemical accuracy). **This gate is the definition of done.**
- **G3 — agreement with exact Krylov.** At matched Δt/depth and high shots, SKQD agrees with the
  exact-Krylov energy to < 1.6 mHa on H₄. **FINDING (gate revised during implementation):** the
  original wording ("SKQD lands at or above the exact-Krylov energy") was wrong. With enough samples
  SKQD's determinant subspace spans *more* than the depth-dimensional span{|ψ_k⟩}, so it reaches
  essentially full FCI, which lies *below* the finite-M Krylov estimate (measured: SKQD −2.16638745
  vs exact-Krylov −2.16638436 Ha on H₄, depth 8). The real variational floor is **FCI** (G1), not
  the dense Krylov energy; the revised gate checks the FCI floor + ≤1.6 mHa agreement.
- **G4 — monotone improvement.** Increasing the shot budget at fixed depth does not raise the energy
  beyond noise: `E_skqd(more shots) ≤ E_skqd(fewer shots) + tol`, tol set from the sampling spread.
- (Stretch, not a gate) qDRIFT/Trotter-compiled sampling instead of exact-amplitude sampling, to
  connect to `trotter_krylov.py` and quantify the depth/accuracy trade-off from §III of 2508.02578.

## 6. Implementation plan (test-first)

1. Write `tests/test_skqd_spec.py` encoding G1–G4 (initially failing).
2. Add `hybrid_quantum_solver/skqd.py` composing existing propagation + subspace eigensolve + FCI.
3. Iterate to green via `make gates` (own process — shares pyscf, so keep it out of the block2 group).

## 7. Out of scope

- Hardware execution and device-noise sampling (a later spec).
- qDRIFT circuit synthesis and depth reduction (stretch above; its own spec).
- Configuration recovery / self-consistent determinant augmentation (the SQD-side improvement).
- Any system where FCI is intractable (no reference → not this spec; that is the DMRG-referenced
  TM-active-space backlog item's job).

## 8. Caveats and risks

- **R1 — concentration may fail at stretched geometries.** Strongly multireference points spread the
  ground state over many determinants, so a fixed shot budget under-samples and G2 can fail.
  *Mitigation:* gate at geometries where the state is concentrated (near equilibrium); record the
  stretched-geometry failure as a finding (it is the honest content of the concentration hypothesis).
- **R2 — sampling RNG flakiness.** Seed the RNG; set tolerances from the measured spread, not by
  hand-tuning, so the gate stays falsifiable.
- Honest limitation: classical sampling from the exact statevector is an idealization — it omits the
  Trotter/qDRIFT and hardware-noise errors that dominate a real device.

## 9. Deliverables

- `hybrid_quantum_solver/skqd.py` — `SampleKrylovSolver`.
- `benchmark_skqd.py` — convergence CSV vs exact-Krylov and FCI.
- `tests/test_skqd_spec.py` — gates G1–G4.
- Results summary (with §2/§7 caveats) in the PR description.
