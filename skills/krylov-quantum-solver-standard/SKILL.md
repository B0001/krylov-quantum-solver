---
name: krylov-quantum-solver-standard
description: qDRIFT real-time quantum Krylov diagonalization for molecular ground states in a CAS active space (PySCF + Qiskit Nature), with mandatory self-checks against CASCI/FCI. Use when asked to estimate an active-space ground-state energy with a Krylov/qDRIFT quantum-subspace method, or to build, simulate, or review such a pipeline.
version: 1.1.0
author: Benjamin
license: MIT
---

# Krylov Quantum Solver — Standard Tier

Use this skill to estimate a molecule's active-space ground-state energy with **real-time quantum
Krylov diagonalization**, where each time evolution e^{−iΔtH} is approximated by **qDRIFT**
(Campbell, PRL 123, 070503 (2019); Krylov: Klymko et al., PRX Quantum 3, 020323 (2022)).

**Rule zero:** never report an energy without the reference check in step 6. Every number this
skill produces is compared against CASCI/FCI for the same active space. If the check can't run,
say so instead of reporting the number.

## Workflow

1. **Active space.** Run `PySCFDriver` and then `ActiveSpaceTransformer(n_el, n_orb)`, and map the
   result with `JordanWignerMapper`. Keep `sum(problem.hamiltonian.constants.values())` (nuclear
   repulsion plus frozen core) as an offset; the qubit operator does not include it.
2. **Reference state.** Use the Hartree–Fock state (`qiskit_nature...HartreeFock` circuit). Never
   use |00…0⟩: under Jordan–Wigner that is the 0-electron vacuum, and time evolution conserves
   electron number, so a Krylov space built from it can't reach the ground state. Check that the
   reference lies in the (N_α, N_β) sector.
3. **qDRIFT walk.** Set λ = Σ|c_j|, Δt = π / (spectral width), and τ = λΔt/L. Each sample applies
   **only the sampled term**: e^{−iτ·sgn(c_j)·P_j} with j drawn with probability |c_j|/λ. One
   random walk supplies a **whole** basis: |φ_k⟩ is read off after every L rotations.
4. **Subspace matrices.** S = ⟨V†V⟩ and H = ⟨V†HV⟩, averaged over R independent walks. Averaging
   whole-basis walks keeps (H, S) a Gram pair, which guarantees the variational bound. Averaging
   each (i, j) entry over separately drawn walks does not, so don't do that.
5. **Solve.** Use thresholded canonical orthogonalization: diagonalize S, drop eigenvalues below
   ~1e-6 × max, project H, then diagonalize. Never call `scipy.linalg.eigh(H, S)` directly,
   because S from a Krylov basis is near-singular and that call either fails or returns values
   below the true minimum.
6. **Validate.** Compute the exact (N_α, N_β)-sector minimum, plus PySCF `mcscf.CASCI` on the same
   orbitals as an independent check. Then assert:
   (a) qubit sector minimum == CASCI within 1e-6 Ha,
   (b) E ≥ E_CASCI (variational),
   (c) report E − E_CASCI. Claim chemical accuracy only if it is below 1.6 mHa.

Working template: `examples/basic_h2o_active.py`. It is standalone, runs its own checks, and
exits non-zero if any check fails.

## Pitfalls this skill exists to prevent

| Mistake | Symptom |
|---|---|
| Start from `psi0[0] = 1` | Energy of the 0-electron sector (−73.30 Ha vs −74.96 Ha for H₂O CAS(2,2)) |
| Sample a term, then evolve under the full H | Every Krylov vector is a copy of the same state up to a phase; S has rank 1 |
| `eigh(H_sub, S_sub)` with no threshold | `LinAlgError: not positive definite` |
| Full-space Hamiltonian compared to a CAS energy | Comparing two different problems; H₂O full space is 14 qubits and ~4 GB per dense matrix |
| Forgetting the constant offset | Energy off by roughly the nuclear repulsion |
| Quoting cost as O(λ²/ε²) | qDRIFT needs N = ⌈2λ²t²/ε⌉ samples. Cost grows with λ² and t², and only linearly with 1/ε |

## Particle-number leakage (read this before trusting the floor)

A single Pauli string doesn't conserve electron number, so each walk leaks out of the N-electron
sector (H₂O CAS(4,4): up to 15 % of the norm at L=50, 0.5 % at L=2000). Without postselection the
guaranteed floor is the minimum over **all** electron numbers (the Fock space). When that
minimum lies in another sector, the estimate can drop **below** CASCI: 10–42 mHa below in a
μN̂-shifted test. **Project each state onto the N-electron sector** (a particle-number symmetry
check; it costs discarded shots on hardware). That restores CASCI as a rigorous floor, and it also
lowered the error by about 10–25 % in our tests.

## Measured performance (exact statevector, H₂O/STO-3G, M=6, R=20, median of 3 seeds)

| Active space | Qubits | λ | L=50 | L=400 | L=2000 |
|---|---|---|---|---|---|
| CAS(2,2) | 4 | 2.44 | 0.28 mHa | 0.037 mHa | 0.009 mHa |
| CAS(4,4) | 8 | 11.2 | 6.9 mHa | 3.4 mHa | 0.89 mHa |

The L needed for chemical accuracy grows with λ. Choose L from a convergence scan; don't assume a
value.

## Honest scope

- This combines known methods correctly. It is not a new algorithm.
- Results come from exact statevector simulation of the qDRIFT channel. Shot noise on H and S
  removes the variational guarantee unless S is thresholded at the noise level.
- Dense diagonalization (used for Δt and the reference) limits this template to about 12–14
  qubits. It does **not** reach "100 qubits". Larger systems need a DMRG reference and a Lanczos
  estimate of the spectral width.
- Minimal basis and small active spaces: this is a validation-scale tool, and it shows no quantum
  advantage.

## Requirements

`pyscf>=2.4`, `qiskit>=1.0`, `qiskit-nature>=0.7`, `numpy>=1.24`, `scipy>=1.10,<1.16`. scipy
1.16+ breaks the qiskit-nature import chain.
