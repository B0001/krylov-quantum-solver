# Krylov Quantum Solver — Scientific Audit & Refactor Plan

**Status:** the pipeline does not compute the quantities it claims to compute. Every headline
number in the README is a numerical artifact of a broken physics core, not a measurement of
electron correlation. This document records the evidence (all reproduced empirically against
exact references in a `qiskit 2.4.2 / pyscf 2.13` environment) and lays out a phased plan to
make the project compute real, validated chemistry.

The good news: the *scaffolding* (PySCF integral extraction, a CLI, a sweep harness, CSV
telemetry) is reasonable. The physics engine in the middle is what needs to be rebuilt, and it
can be rebuilt on solid, well-known foundations.

---

## 1. What was verified

All claims below were checked by running the repo's own classes on H₂/STO-3G, where the exact
answer is known, and on the supplied benchmark CSVs.

| # | Claim in repo | Reality (measured) | Severity |
|---|---------------|--------------------|----------|
| 1 | "Rigorous multi-body Jordan–Wigner mapping" | For H₂/STO-3G the correct JW Hamiltonian has **15 Pauli terms, ground state −1.852388 Ha (= FCI)**. The repo produces **4 terms, lowest eigenvalue −1.547 Ha**, ~0.30 Ha (191 kcal/mol) off, and **no eigenvalue anywhere near FCI**. | Critical |
| 2 | Spin-orbital Hamiltonian on "16 qubits" | Integral indices are **spatial** orbitals but the code allocates `2×n_spatial` qubits. Half the qubits are never touched; **electron spin is dropped entirely**. "CAS(8,8) → 16 qubits" is really an ~8-mode, spinless, malformed operator. | Critical |
| 3 | "Resolves exact ground states via quantum Krylov subspace projection" | The "Krylov basis" appends **one random single-Pauli rotation with an infinitesimal angle** per step. The states are ~identical: overlap matrix **S is rank-1 (cond ≈ 6×10³⁰⁰)**, off-diagonal overlaps ≈ 0.9999. It spans no useful subspace. | Critical |
| 4 | Reference state for the molecule | The sampler starts from the **vacuum \|00…0⟩ (zero electrons)**. `⟨0\|H\|0⟩ = 0`, which is why low-M energies collapse to ~10⁻¹⁵ Ha. | Critical |
| 5 | "M≤16 lacks depth → nominal zero baseline; M=32 resolves −11.39 Ha of correlation" | The identical value **1.2906×10⁻¹⁵ Ha for M=4, 8, 16** is machine-epsilon collapse, not physics. The repo's *own* `sweep_hybrid_solver.py` calls this a **"Krylov basis collapse … solver fell through to the zero-vector null space."** The README reinterprets the same artifact as a feature. | Critical (integrity) |
| 6 | Energies are physical measurements | With injected noise the solver returns **−5, −39, −191, −813 Ha** — far below the true minimum (−1.852 Ha). **A variational/Krylov method can never go below the true ground state.** The README's −48.7 Ha is the same pathology. | Critical |
| 7 | "QCIVET cryptographic integrity verification … detects hardware noise drift" | It only checks `max\|A − Aᵀ\| > 1e-6`. An **absurd but symmetric** H (all −10⁶ Ha) → `VERIFIED_SAFE`; a **correct** H with a 10⁻⁵ asymmetry → `REJECTED`. The SHA-256 hash is computed and **never used** to verify the returned matrices. | High (misleading) |
| 8 | "Hardware noise" model | `h_matrix += np.random.normal(0, noise_variance, …)` — arbitrary Gaussian numbers added asymmetrically to the Hamiltonian. Not shot noise, not a device channel. (`noise_variance` is also passed as the **std-dev** arg, so it's a std-dev, not a variance.) | High |
| 9 | "Tensor-network local optimizations" / IBM QPU gateway | `tensor_network_contractor.py` is a zeros-filled MPS toy **never imported by the solver**. `ibm_quantum_gateway.py` runs one fixed-angle `EfficientSU2` estimate and **fabricates** the whole subspace matrix from that single number (`h_val = exp_value*0.1*(1/(i-j+1))`). | High |
| 10 | The code runs as written | The orchestrator's own `__main__` calls `execute_molecular_query(...)`, **which does not exist** (`AttributeError`). The demo path was never executed. | Medium |
| 11 | Tests validate the physics | Unit tests assert the broken mapping's **own output** (e.g. `scale = 0.5·w/16`). No test compares against FCI or a reference library, so 100% of tests pass while the physics is wrong. | High |
| 12 | Benchmark timings | Up to **68,327 s (19 h) per point** for a 16-qubit statevector job — implausible for the described computation; the O(M²) grid recomputes full statevectors. Numbers are not reproducible as reported. | Medium |

### Reproduce the audit
```bash
# in the chem env
python - <<'PY'
# (the three verify scripts used for this audit are in the response;
#  key result:)
# correct JW  : 15 terms, Emin = -1.852388 Ha  == FCI
# repo mapping :  4 terms, Emin = -1.547168 Ha  != FCI, no eig near FCI
# <0|H|0> = 0 ; S eigenvalues ~ rank-1 ; noisy E -> -812 Ha (below true min)
PY
```

---

## 2. Root cause

Four independent errors compound, and each alone is fatal:

1. **Wrong Hamiltonian.** The hand-rolled JW transform handles only two ERI index patterns
   (`p==r & q==s`, and 4-distinct-index) with ad-hoc `/16` factors and a single Z-corridor.
   The correct two-body JW expansion (8 terms per integral, two Z-strings, specific signs) and
   the spin-orbital structure are missing.
2. **Wrong reference.** Time-evolution Krylov needs a state with non-zero overlap on the ground
   state (Hartree–Fock). Starting from the particle-number-zero vacuum guarantees ~0 energy.
3. **No real subspace.** qDRIFT is a *stochastic approximation to `e^{-iHt}` that must be
   averaged over many realizations*; using one short sequence of infinitesimal single-Pauli
   rotations as successive "basis vectors" produces near-duplicate states and a singular `S`.
4. **Ill-posed solve + fake noise.** Adding asymmetric Gaussian noise to `H` and solving a
   generalized eigenproblem with a near-singular `S` is unbounded below, so the "energy" is
   numerical garbage that the README narrates as correlation capture.

Everything else (QCIVET, the MPS module, the IBM gateway, the 19-hour timings) is decorative and
should not be confused with the scientific core.

---

## 3. Refactor plan (phased)

Each phase ends with a **validation gate** — an automated test that must pass before moving on.
Do not optimize, parallelize, or re-benchmark anything until Phase 1's gate is green.

### Phase 0 — Honesty & a ground-truth harness (½ day)
- Mark the current README "Empirical Benchmarks" section as **retracted/known-artifact** until
  re-derived. Do not delete the data — keep it as a regression fixture of what *not* to produce.
- Add `tests/test_reference_energies.py` that builds the qubit Hamiltonian and asserts the
  ground state equals **FCI within 1e-6** for H₂, H₄ (chain), and LiH. **It must fail on today's
  code** — that failing test is the definition of "done" for Phase 1.
- Delete or quarantine dead/decorative code: `tensor_network_contractor.py`, the fabricated
  matrix logic in `ibm_quantum_gateway.py`, and the broken `__main__` demo block.

### Phase 1 — Correct the physics core (the critical path)  ✅ DONE
> Implemented in `hybrid_quantum_solver/molecular_hamiltonian.py`; gated by
> `tests/test_reference_energies.py` (8/8 passing). Verified totals: H₂ −1.137284,
> H₄ −2.166387, LiH −7.882324 Ha (= PySCF FCI to <1e-6); the active-space path matches
> PySCF CASCI; the Hartree-Fock reference reproduces RHF; the old compactor is pinned as
> failing. Next: Phase 2 wires this Hamiltonian + HF state into the subspace solver.

- **Replace the hand-rolled mapping with a vetted library.** Use **Qiskit Nature**
  (`PySCFDriver` → `ElectronicStructureProblem` → `JordanWignerMapper`) or **OpenFermion**
  (`InteractionOperator` → `jordan_wigner`). Both reproduce −1.852388 Ha for H₂/STO-3G out of
  the box (verified). This deletes `AdvancedStochasticCompactor` entirely.
- **Use the proper spin-orbital count** (`2 × n_active_spatial` modes, with the α/β structure the
  library handles) and the **Hartree–Fock reference state** (occupy the N lowest spin-orbitals:
  `|1…1 0…0⟩` in JW), not the vacuum.
- **Carry the energy frame.** Add `e_core` (frozen-core + nuclear repulsion) back to the
  subspace eigenvalue so reported totals are physical. (`sweep_hybrid_solver.py` already flags
  this `frame_mismatch`.)
- **Validation gate:** Phase 0's reference test passes (FCI to 1e-6) for H₂, H₄, LiH.

### Phase 2 — A real quantum subspace solver  ✅ DONE (exact-evolution core)
> Implemented in `hybrid_quantum_solver/quantum_krylov_solver.py`; gated by
> `tests/test_krylov_convergence.py`. Real-time Krylov basis |φₖ⟩ = e^{−ikΔtH}|HF⟩ with
> Δt = π/width (width via Lanczos), Hermitian H/S, thresholded canonical orthogonalisation,
> and a guaranteed variational floor. Verified: H₂ hits FCI at M=2; LiH converges to
> <0.1 mHa with rank growing 1→8; **no estimate ever drops below FCI** (vs the old code's
> −800 Ha), and energy decreases on every rank-increasing step. Remaining for hardware:
> swap exact `expm_multiply` for Trotter/qDRIFT circuits (folds into Phase 3 noise work)
> and rewire `run_hybrid_solver.py`/`execute_subspace_sweep` onto this solver.

- **Pick a published method and implement it faithfully.** Recommended: **real-time quantum
  Krylov / quantum subspace expansion** — basis `|ψ_k⟩ = e^{-i k Δt H}|ψ_HF⟩`, with `Δt` set
  from the spectral range (`Δt ≈ π/(E_max−E_min)`), `H_ij = ⟨ψ_i|H|ψ_j⟩`, `S_ij = ⟨ψ_i|ψ_j⟩`.
  (Refs: Klymko et al. 2022 *PRX Quantum* "real-time evolution quantum Krylov"; Stair et al.
  2020 multireference selected QK; Motta et al. QITE/QLanczos.)
- If you keep **qDRIFT**, use it correctly: each evolution `e^{-iHΔt}` is approximated by `L`
  randomly sampled single-Pauli rotations (angle `τ = λΔt/L`, sign of the coefficient), and the
  *expectation is averaged over many stochastic realizations*. One sequence ≠ one basis vector.
- **Solve the generalized eigenproblem with a principled regularizer and a variational floor:**
  thresholded eigendecomposition of `S` (drop singular values `< ε`), project `H` into the kept
  subspace, then take the Rayleigh-quotient minimum. Assert `E ≥ E_FCI − tol` by construction.
- **Validation gate:** on H₂ and LiH, energy **converges monotonically to FCI** as `M` grows
  (e.g. error < 1 mHa by `M ≈ 8`), and **never** drops below FCI.

### Phase 3 — Real noise & error mitigation (only if a QPU is a goal)  ✅ DONE (statistical core)
> Implemented in `hybrid_quantum_solver/noise.py` + the solver's noise-aware solve; gated by
> `tests/test_noise_resilience.py`. The fake asymmetric-Gaussian-on-H is gone; shot noise is now
> modelled as Hermitian-symmetric sampling error (~1/√shots) with a noise-aware overlap cutoff
> (Epperly et al.), so error degrades gracefully and stays bounded (tens of mHa) instead of the
> old −800 Ha. QCIVET is retired; `noise.py` provides `hermitize` and a real Aer `NoiseModel`
> builder. The hardware-facing layer is now in `hybrid_quantum_solver/trotter_krylov.py` (gated by
> `tests/test_trotter_circuit.py`): the Krylov basis is built from genuine Suzuki-Trotter
> `PauliEvolutionGate` circuits (the correct replacement for the single-Pauli qDRIFT step), and
> `estimate_energy_aer` runs expectation values through qiskit-aer exactly, under shot noise, and
> under a device `NoiseModel`. The full **on-hardware** solver is in
> `hybrid_quantum_solver/hardware_krylov.py` (gated by `tests/test_hardware_krylov.py`): every
> Hᵢⱼ/Sᵢⱼ is *measured* by an ancilla Hadamard test with controlled Trotter evolution (not a
> statevector inner product), reproducing FCI exactly in the noiseless limit and staying bounded
> (mHa-scale) under shot/device noise. **Zero-noise extrapolation** (global circuit folding +
> linear extrapolation in `noise.py`; opt-in via `HardwareKrylovSolver(zne_scale_factors=[1,3,5])`)
> cuts the device-noise error ~2.4× in the good-hardware regime (gated by
> `tests/test_hardware_krylov.py`) — with the honest caveat that linear ZNE has little to extrapolate
> once the deep controlled circuits are heavily decohered. Still remaining: readout-error mitigation
> (folding does not amplify readout) and scaling the O(M²) measurement beyond small systems.


- Replace additive-Gaussian-on-H with **physical noise**: finite-shot sampling via
  `Estimator`/`Sampler`, and device noise via a **Qiskit Aer `NoiseModel`** (depolarizing,
  thermal relaxation T₁/T₂, readout error).
- Replace "QCIVET" with **legitimate** subspace-noise handling: Hermitize `H`, `S` by averaging
  with the conjugate transpose; enforce `S ⪰ 0`; apply thresholding (quantum subspace methods
  are provably noise-resilient through this). Optionally add **zero-noise extrapolation** or
  **measurement-error mitigation**. If you want provenance/audit, hash the payload *and verify
  the hash on return* — and stop calling it noise detection.
- **Validation gate:** under a calibrated Aer noise model, the mitigated energy degrades
  *gracefully and stays above FCI*; rejection logic triggers on real shot/device noise, not on a
  1e-5 asymmetry.

### Phase 4 — Honest benchmarking & real discovery potential  ◐ STARTED
> `benchmark_krylov.py` produces an honest error-vs-FCI table (H₂ equilibrium + stretched, H₄,
> LiH; all within ~0.1 mHa and above the variational floor), replacing the retracted README
> table. The CLI/pipeline are rewired and report energies in the correct frame vs CASCI.
> The **N₂ dissociation rung** is now in `benchmark_n2.py` (CAS(6,6), 12 qubits): across the
> bond-breaking coordinate the Hartree–Fock error grows 101 → 625 mHa (single-reference
> breakdown) while quantum Krylov stays ≤ 0.65 mHa vs exact CASCI — the exact reference
> cross-validated against an independent PySCF CASCI solve (−107.62310177 Ha, agree to 1e-6).
> DMRG (block2) is the reference for active spaces beyond FCI's reach. **block2 0.5.3 is now
> installed and the DMRG path is executed and validated** (`hybrid_quantum_solver/dmrg_reference.py`,
> `reference_energy(method="auto")` routes through DMRG when available). The DMRG-backed ladder is
> in `benchmark_dmrg.py` (stretched N₂, 6-31g, growing active space): DMRG reproduces exact FCI to
> ~1e-10 Ha at CAS(6,6)→(10,10) and ~1e-8 Ha at CAS(12,12) (bond-dim-400 truncation), and at
> **CAS(14,14)** (≈1.18e7 determinants, FCI skipped) **DMRG carries the reference alone** — the
> beyond-FCI rung the earlier sandbox could not produce. Quantum Krylov tracks the reference where
> the statevector fits (≤16 qubits: 0.016/0.285 mHa at CAS(6,6)/(8,8)); larger active spaces are
> correctly out of reach of the statevector solver. **Resource accounting** is in
> `benchmark_resources.py` + `HardwareKrylovSolver.resource_report` (gated by
> `tests/test_hardware_krylov.py`): e.g. N₂ CAS(6,6) needs ~6.5k CX per Trotter step and ~240k CX
> for the deepest M=6 Hadamard-test circuit — honestly far beyond NISQ, which motivates shallower
> evolution (better Trotter/qubitization). Remaining: a real transition-metal active space vs DMRG,
> and choosing a validated discovery target (CIF-as-molecule caveat below). NOTE: block2 is an
> optional reference dependency (not in `requirements.txt`); the code degrades to exact FCI without it.


- **Benchmark ladder, each vs. a classical reference:** H₂ → LiH → H₄/H₆ chains (vs FCI) → N₂
  in a CAS (vs DMRG/CCSD(T)) → *then* a transition-metal target. Report error vs reference,
  active-space and basis convergence, `M`-convergence, with seeds and error bars.
- **Fix the cost model:** cache the `M` statevectors once (don't rebuild per `(i,j)`), evaluate
  expectations with sparse Pauli measurement, and report honest resources (qubits, depth, #
  measurements, shots) instead of wall-clock numbers that imply 19-hour 16-qubit jobs.
- **To genuinely raise discovery potential,** aim the validated method at a question it can
  actually answer and *check*: e.g. correlation energy of a carefully chosen NbN active space
  benchmarked against DMRG/AFQMC, or strongly-correlated model Hamiltonians (Hubbard) where
  exact/near-exact references exist. Claim only what a classical reference can corroborate.

---

## 4. File-by-file disposition

| File | Action |
|------|--------|
| `hybrid_quantum_solver/orchestrate_hybrid_pipeline.py` | Gut `AdvancedStochasticCompactor` (→ Qiskit Nature/OpenFermion). Rewrite `StabilizedSubspaceShifter` solve to enforce the variational floor. Delete `QCIVETGuard` or demote to an optional provenance log. Fix the broken `__main__`. |
| `hybrid_quantum_solver/quantum_sampler.py` | Rewrite: HF reference state, true `e^{-ikΔtH}` basis, cached statevectors, `Estimator`-based expectations. |
| `hybrid_quantum_solver/noise_resilient_compactor.py` | Remove; replace with Aer `NoiseModel` + mitigation. |
| `hybrid_quantum_solver/ibm_quantum_gateway.py` | Remove the fabricated-matrix logic; if kept, it must measure each real `H_ij`/`S_ij`, and update the retired `ibmq_qasm_simulator` backend reference. |
| `hybrid_quantum_solver/chemistry_gateway.py` | Mostly keep (PySCF extraction is fine); ensure it returns spin-orbital integrals and the HF occupation for the reference state. |
| `tensor_network_contractor.py` | Delete (unused, decorative) or rebuild into an actual DMRG/MPS reference if you want a classical baseline. |
| `test_chemistry_mappings.py` | Replace self-referential asserts with FCI/reference comparisons. |
| `README.md` | Rewrite claims to match validated results; remove "enterprise/oracle/QCIVET" framing that obscures what the code does. |

---

## 5. One-paragraph summary for the README (suggested honest framing)

> A hybrid quantum–classical pipeline that extracts active-space integrals with PySCF, maps them
> to qubits with a validated Jordan–Wigner transform (Qiskit Nature), and estimates ground-state
> energies via a real-time quantum Krylov subspace expansion. Validated to reproduce FCI within
> 1 mHa for H₂, LiH, and H₄, with documented convergence in subspace dimension and a calibrated
> shot/device-noise study. Transition-metal systems are a research target, benchmarked against
> DMRG where tractable — not yet a validated result.

That is a smaller claim than "bypasses the FCI wall for superconductors," but it would be **true**,
reproducible, and a real foundation to build discoveries on.
