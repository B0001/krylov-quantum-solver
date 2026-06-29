# SPEC: Quantum Krylov eigenstates yield molecular properties (dipoles, oscillator strengths)

**Status:** CLOSED — gates G1–G4 PASS (2026-06-29); `eigenstates` + `qksd_properties` +
`build_dipole_operators` merged. No regression (excited/Krylov/noise/GPU/reference suites green).
Builds on [`SPEC_qksd_excited.md`](SPEC_qksd_excited.md).

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

The Krylov subspace carries more than energies: each Ritz eigenstate |Ψ_m⟩ = Σ_i C[i,m]|φ_i⟩ is a
genuine Hilbert-space vector, so ⟨Ψ_m|Ô|Ψ_n⟩ gives properties. Claim: the exact-evolution
`QuantumKrylovSolver` reproduces, to chemical/property accuracy, the FCI **permanent dipole**,
**transition dipoles**, and **oscillator strengths** of small molecules — recovering bright
transitions, *and* recovering symmetry-forbidden (dark) transitions as ≈0. The claim is false if a
QKSD property deviates from the dense-diagonalization reference beyond tolerance, or if the dipole
operator is mis-constructed (a nonzero dipole for a centrosymmetric molecule).

## 2. Background and honest framing

- **Prior art / reference.** Oumarou et al., *Molecular Properties from Quantum Krylov Subspace
  Diagonalization*, arXiv:2501.05286 — derives QKSD energy derivatives and relaxed 1-/2-RDMs of
  Krylov eigenstates (incl. excited states "if a suitable reference state is used"), and notes the
  framework gives "any property... including dipole moments". The hardware-facing machinery is QSP
  eigenstate prep + RDM measurement; the **classically-checkable core** is the matrix element
  ⟨Ψ_m|Ô|Ψ_n⟩, which is what we validate. Builds directly on the excited-state rung
  ([`SPEC_qksd_excited.md`](SPEC_qksd_excited.md)).
- **Ground truth.** Dense diagonalization of the same qubit Hamiltonian and the same dipole
  operators (`MolecularHamiltonian`, O(4ⁿ)) — the FCI reference for properties, the same idiom in
  which `ground_state_energy()` is the FCI reference for energies. The dipole operator is
  additionally anchored to PySCF (the qubit ground-state z-dipole is the correlated FCI dipole,
  LiH/STO-3G ≈ −1.817 a.u. vs RHF −1.912).
- **What we can claim if gates pass.** QKSD recovers FCI dipoles / transition dipoles / oscillator
  strengths on these systems, and the dipole operator is the correct physical observable (zero for
  H₂ by symmetry, nonzero and FCI-accurate for polar HeH⁺).
- **What we cannot claim.** (a) No quantum advantage — exact statevector, tiny systems, a
  correctness target. (b) Full orbital space only here; the active-space dipole (frozen-core
  contribution folded in) is a follow-up. (c) Relaxed-RDM nuclear gradients and the QSP/RDM
  measurement cost analysis of 2501.05286 are out of scope (§7) — we validate the property values,
  not the hardware measurement scheme.

## 3. Approach

Reuse the excited-state subspace verbatim. Add `ritz_pairs` (Ritz values **and** eigenvectors in
the Krylov basis) and `QuantumKrylovSolver.eigenstates` (reconstruct |Ψ_m⟩ as statevectors). Build
the per-axis total dipole operator `nuclear·I − electronic_op` (JW-mapped) with
`build_dipole_operators`. Properties are then `⟨Ψ_m|μ̂|Ψ_n⟩`: permanent dipole = diagonal,
transition dipole = ground→n off-diagonal, oscillator strength `f_n = (2/3)(E_n−E_0)|μ_{0n}|²`.
Reference = the identical matrix elements between the dense-diagonalized exact eigenstates that are
reachable from |HF⟩ (nonzero HF overlap, as in the excited-state spec).

## 4. Public interface

```
molecular_hamiltonian.build_dipole_operators(atom, basis, charge, spin, mapper)
    -> list[SparsePauliOp]                              # [mu_x, mu_y, mu_z], total dipole per axis
quantum_krylov_solver.ritz_pairs(H, S, threshold, noise_floor)
    -> (vals: np.ndarray, C: np.ndarray, rank: int)     # Ritz values + Krylov-basis eigenvectors
QuantumKrylovSolver.eigenstates(krylov_dim, n_states=None)
    -> (energies: list[float], states: np.ndarray)      # states[m] is |Psi_m> (k x N)
qksd_properties.property_matrix(states, operator)        -> (k x k) <Psi_m|O|Psi_n>
qksd_properties.transition_dipoles(states, dipole_ops)   -> (3 x k x k)
qksd_properties.oscillator_strengths(energies, states, dipole_ops) -> (k,)
```

`ritz_spectrum` / `solve_excited` / `solve_generalized_eig` keep their signatures (now delegating
to `ritz_pairs`), so every existing test is untouched.

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_qksd_properties_spec.py` (test-first). Exact statevector, pyscf/qiskit, full
orbital space. Reference = reachable exact eigenstates of the same operators. Tol `1e-3` a.u. for
dipoles unless noted.

- **G1 — permanent dipole vs FCI.** QKSD ground-state permanent dipole matches the
  dense-diagonalization value: HeH⁺ z-dipole ≈ 1.07 a.u. (nonzero), H₂ ≡ 0 by symmetry (`< 1e-6`).
- **G2 — transition dipoles vs FCI.** On HeH⁺, the ground→excited transition-dipole magnitudes to
  the reachable states match FCI within tol, and the first transition is **bright**
  (|μ_{01}| > 0.5 a.u., measured ≈ 0.85).
- **G3 — oscillator strengths + dark recovery (definition of done).** HeH⁺ oscillator strengths
  match FCI within `1e-3`; H₂'s reachable g→g transition is recovered **dark** (`f < 1e-6`,
  |μ_{01}| < 1e-6) — the symmetry the operator must respect.
- **G4 — correctness invariants (can't-be-faked).** The property matrix is Hermitian
  (`||O − O†|| < 1e-9`) and the Krylov eigenstates are normalized (`|⟨Ψ_m|Ψ_m⟩ − 1| < 1e-9`); a
  violation is a reconstruction/normalization bug, the property analogue of the energy variational
  floor.
- (Stretch, not a gate) relaxed-RDM nuclear gradients and the QSP/RDM measurement-variance analysis
  of arXiv:2501.05286; active-space dipoles.

> Definition of done: **G3**. If a property needs deeper Krylov depth than the energy to converge
> (as excited *energies* did, SPEC_qksd_excited G3), **raise the depth and record it** — that depth
> is the finding, not a failure.

## 6. Implementation plan (test-first)

1. Write `tests/test_qksd_properties_spec.py` encoding G1–G4 (initially failing).
2. Add `ritz_pairs` + `eigenstates`, `build_dipole_operators`, and `qksd_properties`.
3. Iterate to green via `make gates` (own process; pyscf/qiskit, no block2).

## 7. Out of scope

- Relaxed 1-/2-RDM nuclear gradients (the actual focus of arXiv:2501.05286) and geometry
  optimization.
- The QSP eigenstate-preparation / coherent-RDM measurement scheme and its variance analysis
  (hardware measurement cost — we validate the property *values*, not the measurement protocol).
- Active-space dipole operators (frozen-core contribution); shot-noise/Trotter property degradation.

## 8. Caveats and risks

- **R1 — properties converge slower than energies.** Off-diagonal matrix elements may need a
  larger Krylov dimension than the ground-state energy. *Mitigation:* gate at a depth where the
  reachable subspace is saturated (as in the excited-state spec); record the depth.
- **R2 — degenerate states make individual transition dipoles gauge-dependent.** Only the summed
  strength over a degenerate manifold is well-defined. *Mitigation:* gate on non-degenerate
  low-lying states (HeH⁺) and on the symmetry-forced dark value (H₂).
- Honest limitation: exact statevector, minimal-basis tiny molecules — reproduces known dipoles, a
  correctness target, not a hardware or large-system result.

## 9. Deliverables

- `hybrid_quantum_solver/quantum_krylov_solver.py` — `ritz_pairs`, `eigenstates`.
- `hybrid_quantum_solver/molecular_hamiltonian.py` — `build_dipole_operators`.
- `hybrid_quantum_solver/qksd_properties.py` — `property_matrix`, `transition_dipoles`,
  `oscillator_strengths`.
- `tests/test_qksd_properties_spec.py` — gates G1–G4.
- Results summary (with §2/§7 caveats) in the PR description.
