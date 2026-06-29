# SPEC: Real-time quantum Krylov resolves the low-lying excited spectrum, not just the ground state

**Status:** CLOSED — gates G1–G4 PASS (2026-06-29); `solve_excited` merged into
`quantum_krylov_solver.py`. No ground-state regression (existing Krylov/noise/GPU suites green).
Measured depth ladder recorded in §5 G3.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

The validated `QuantumKrylovSolver` already builds the projected subspace matrix `H_proj`; today it
discards every Ritz value except the lowest. Claim: the *full* Ritz spectrum of that same subspace
converges, from above, to the low-lying **excited** eigenvalues of the molecular Hamiltonian that
are reachable from the Hartree–Fock reference — and the first excitation energy (a physical
observable, not just an energy) lands within chemical accuracy. The claim is false if a Ritz value
ever dips below the corresponding exact eigenvalue, or if the lowest few reachable eigenvalues are
not recovered to chemical accuracy at attainable Krylov depth.

## 2. Background and honest framing

- **Prior art / reference.** Excited-state quantum Krylov subspace diagonalization — Klymko et al.,
  *Quantum Krylov subspace algorithms for ground and excited state energy estimation*,
  arXiv:2109.06868 (PRX Quantum 3, 020323, 2022); Cortes & Gray, PRA 105, 022417 (2022). The
  ground-state primitive these build on is the one already validated here
  (`quantum_krylov_solver.py`).
- **Ground truth.** The exact spectrum is the dense eigendecomposition of the *same* qubit
  Hamiltonian (`MolecularHamiltonian`, O(4ⁿ), small systems only). Real-time Krylov from |HF⟩ spans
  span{Hᵏ|HF⟩}, whose closure is exactly the span of the eigenstates with **nonzero HF overlap**.
  So the honest reference for "what the method converges to" is *the exact eigenvalues whose
  eigenstate has |⟨v_i|HF⟩|² above a small threshold* — no PySCF spin-sector bookkeeping required,
  and reachability is defined by the same operator the solver sees.
- **What we can claim if gates pass.** The exact-evolution QKSD recovers the low-lying reachable
  excited energies and the first excitation gap to chemical accuracy, with every Ritz value
  variationally above its exact target (Cauchy interlacing).
- **What we cannot claim.** (a) No quantum advantage — this is exact statevector evolution on tiny
  systems, a correctness target, not a hardware result. (b) States with **zero** HF overlap (a
  different spatial/spin symmetry than the closed-shell reference) are unreachable from this single
  reference — recovering them needs a multi-reference Krylov space (out of scope, §7). (c) Trotter
  and shot-noise degradation of excited states is a separate, harder question (§7).

## 3. Approach

Reuse the existing subspace build verbatim (`_subspace_matrices`, the noise-aware thresholded
canonical orthogonalization in `solve_generalized_eig`). The only change: after projecting H onto
the well-conditioned subspace, return **all** sorted Ritz values, not just the minimum. The lowest
`k` of them are the excited-state estimates.

Reference (in-test, from the same `MolecularHamiltonian`): dense-diagonalize the qubit Hamiltonian,
compute |⟨v_i|HF⟩|², keep eigenvalues whose overlap exceeds `tol`, sort ascending, lift by the
energy offset. These are the *reachable* exact energies. The first excitation gap is
`reachable[1] − reachable[0]`.

## 4. Public interface

```
quantum_krylov_solver.ritz_spectrum(H, S, threshold, noise_floor)
    -> (np.ndarray ascending Ritz values, rank: int)          # solve_generalized_eig wraps this
QuantumKrylovSolver.solve_excited(krylov_dim, n_states=None)
    -> ExcitedKrylovStep(dim, rank, energies: list[float])    # energies include the offset
```

`solve_generalized_eig` keeps its current `(energy, rank)` signature (now delegating to
`ritz_spectrum`), so the ground-state path and all existing tests are untouched.

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_qksd_excited_spec.py` (test-first). Exact statevector, pyscf/qiskit only (no
block2). Reference = reachable exact spectrum of the same `MolecularHamiltonian`.

- **G1 — variational interlacing (the can't-be-faked invariant).** For H₂ and H₄ at every Krylov
  dimension, the i-th ascending Ritz value is ≥ the i-th ascending exact eigenvalue of the full
  Hamiltonian, to within `1e-9` Ha. A Ritz value below its exact target is a projection/sign bug.
- **G2 — ground-state regression.** `solve_excited(M).energies[0] == solve(M).energy` exactly (the
  excited API must not perturb the validated ground-state path), and matches FCI.
- **G3 — reachable excited energies converge (definition of done).** On H₄ (8 qubits), at attainable
  depth, the lowest `k` reachable reference energies are each matched by a Ritz value within
  `1.6 mHa`, and the number of reachable states is recovered (rank saturates at the reachable
  subspace dimension). `k ≥ 2` (ground + at least one excited). **MEASURED:** H₄'s reachable
  subspace is **12-dimensional**; the 3rd state is the binding constraint — 6.2 mHa off at M=16 but
  < 0.5 mHa by M=20 and < 0.05 mHa by M=24 as the kept rank grows (8→11). Gate runs at M=24, k=3.
  Excited states demonstrably need a deeper Krylov space than the ground state (which is converged
  by M≈8) — that depth ladder is the finding.
- **G4 — first excitation gap.** The Krylov first gap `energies[1] − energies[0]` matches the exact
  reachable gap within `1.6 mHa` on H₂ and H₄.
- (Stretch, not a gate) excited-state degradation under shot noise / Trotter; transition properties
  (oscillator strengths) à la arXiv:2501.05286.

> Definition of done: **G3**. If a reachable excited state proves unreachable at attainable depth
> (overlap too small to resolve before the subspace saturates), **lower `k` / raise the depth and
> record the overlap at which resolution failed** — that crossover is the finding.

## 6. Implementation plan (test-first)

1. Write `tests/test_qksd_excited_spec.py` encoding G1–G4 (initially failing — `solve_excited`
   does not exist yet).
2. Factor `ritz_spectrum` out of `solve_generalized_eig`; add `solve_excited` + `ExcitedKrylovStep`.
3. Iterate to green via `make gates` (own process; shares pyscf, stays out of the block2 group).

## 7. Out of scope

- Multi-reference Krylov spaces (multiple references to reach zero-HF-overlap symmetries).
- Trotter/qDRIFT-compiled and shot-noisy excited states (the ground-state noise gate exists; excited
  states under noise are harder and deserve their own spec).
- Transition matrix elements / oscillator strengths / response properties (arXiv:2501.05286).
- Any system where the dense reference is intractable (no reference → not this spec).

## 8. Caveats and risks

- **R1 — near-degenerate Ritz values mis-pair with the reference.** Matching by sorted order can
  swap near-degenerate states. *Mitigation:* match the lowest `k` as a sorted set within tolerance,
  not index-by-index identity; keep `k` to the cleanly-resolved low-lying states.
- **R2 — depth vs resolution.** Excited states need a larger Krylov dimension than the ground state;
  too small `M` under-resolves them. *Mitigation:* gate at a depth where the rank has saturated the
  reachable subspace; record the depth.
- Honest limitation: exact statevector on minimal-basis tiny molecules — a correctness target that
  reproduces known spectra, not a hardware or large-system result.

## 9. Deliverables

- `hybrid_quantum_solver/quantum_krylov_solver.py` — `ritz_spectrum`, `solve_excited`,
  `ExcitedKrylovStep`.
- `tests/test_qksd_excited_spec.py` — gates G1–G4.
- Results summary (with §2/§7 caveats) in the PR description.
