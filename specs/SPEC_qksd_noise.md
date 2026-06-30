# SPEC: Excited-state quantum Krylov degrades gracefully under shot noise, but more fragile than the ground state

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); no new code (reuses `solve_excited` + the
existing noise machinery). Measured on H₂: excited error shrinks 0.062→0.012 Ha (gap) as shots go
4096→262144, and the excited state is ≈ 24–36× more fragile than the ground state.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

The ground-state Krylov path is already validated under a finite-sampling shot-noise model
(`test_noise_resilience.py`). Extend that to the **excited** spectrum from `solve_excited`. Claim:
under Hermitian-symmetric shot noise of scale σ ≈ 1/√shots, the excited-state energies and the first
excitation gap (a) stay bounded (no blow-up), (b) improve as the shot budget grows, and (c) are
**substantially more noise-fragile than the ground state** — a quantified, falsifiable finding, not
a free lunch. The claim is false if excited errors do not shrink with shots, if they blow up, or if
the excited state is *not* meaningfully more sensitive than the ground state.

## 2. Background and honest framing

- **Prior art / reference.** Lee, Lee & Huh, *Sampling error analysis in quantum Krylov subspace
  diagonalization*, Quantum 8, 1477 (2024); Kirby, *Analysis of quantum Krylov algorithms with
  errors*, Quantum 8, 1457 (2024) (linear error scaling, optimal thresholding); Epperly, Lin &
  Nakatsukasa, SIAM J. Matrix Anal. Appl. 43, 1263 (2022) (already cited in the solver). Builds on
  [`SPEC_qksd_excited.md`](SPEC_qksd_excited.md); reuses the existing noise machinery
  (`noise.shot_noise_sigma`, the Hermitian-symmetric perturbation + noise-aware overlap cutoff in
  `QuantumKrylovSolver`) — **no new solver code**, this is a characterization gate.
- **Ground truth.** The noiseless reachable excited spectrum (dense diagonalization, as in the
  excited-state spec). Errors are averaged over many RNG seeds (the noise is random; a single seed
  is not a gate).
- **What we can claim if gates pass.** Excited QKSD energies/gaps degrade gracefully and predictably
  with shots, and we quantify the excited-vs-ground fragility ratio (≈ 24–36× on H₂).
- **What we cannot claim.** (a) No quantum advantage — idealized i.i.d. Gaussian shot noise on the
  subspace matrix elements, tiny systems. (b) The noise model omits Trotter/qDRIFT and correlated
  device noise (separate specs). (c) Weakly-HF-overlapped excited states can fall *below* the
  noise-aware overlap floor and be lost entirely at low shots (rank collapse) — recorded as a
  finding (§5 G4), the honest failure mode of the noise-aware cutoff.

## 3. Approach

For each (shots, seed): `QuantumKrylovSolver(mh, noise_sigma=shot_noise_sigma(shots), seed=s)
.solve_excited(M, n_states=2)`. Average `|E_excited − ref|`, `|gap − ref_gap|`, and
`|E_ground − ref|` over seeds at each shot budget, and compare across budgets and against the ground
state. Reference = noiseless reachable spectrum from `MolecularHamiltonian` dense diagonalization.

## 4. Public interface

No new public API — reuses `QuantumKrylovSolver.solve_excited` (SPEC_qksd_excited),
`noise.shot_noise_sigma`, and dense-diagonalization references. (If a gate proves the noise path
needs a knob, add it then and record why.)

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_qksd_noise_spec.py` (test-first). H₂ (4 qubits), `krylov_dim = 10`, ≥ 12 seeds.
Exact statevector + shot-noise model; pyscf/qiskit, no block2.

- **G1 — noiseless excited recovery.** At σ = 0, `solve_excited(10, 2)` matches the dense reachable
  excited energy and the gap to `< 1e-6` Ha (ties to SPEC_qksd_excited; the σ→0 anchor).
- **G2 — bounded under noise.** At a modest budget (4096 shots), the mean (over seeds) excited-energy
  error is chemical-to-mHa scale (`< 0.1` Ha) and the excited energy never blows up (stays
  `> ref_excited − 0.5` Ha) — nothing like the old unbounded-noise pathology.
- **G3 — improves with shots (definition of done).** The mean first-excitation-gap error at a high
  budget (262144 shots) is clearly smaller than at a low budget (4096 shots) — graceful degradation,
  the excited-state analogue of `test_noise_resilience`.
- **G4 — excited is more fragile than ground (the finding).** At a fixed budget, the mean
  excited-energy error exceeds the mean ground-energy error by a large factor (`> 5×`; measured
  ≈ 24–36× on H₂). A *recorded* observation, not gated for flakiness: weakly-overlapped excited
  states can fall below the noise-aware overlap floor and be dropped at low shots / shallow depth
  (rank collapses to 1) — the resolution threshold.

> Definition of done: **G3 + G4**. If on some system the excited state is *not* more fragile (e.g.
> a strongly-overlapped, well-separated excited state), lower the factor and record where the
> fragility vanishes — that crossover is the finding.

## 6. Implementation plan (test-first)

1. Write `tests/test_qksd_noise_spec.py` encoding G1–G4 (initially failing only because the averaged
   thresholds are asserted; no new code expected).
2. If a gate is unreachable with the current noise path, add the minimal knob and record why.
3. Iterate to green via `make gates` (own process; pyscf/qiskit, no block2).

## 7. Out of scope

- Trotter/qDRIFT-compiled excited states and correlated device-noise models (separate specs).
- ZNE / error mitigation applied to excited states (the ground-state ZNE path exists; excited-state
  mitigation is its own study).
- Larger multireference systems where the dense reference is intractable.

## 8. Caveats and risks

- **R1 — RNG flakiness.** Single-seed errors are noisy. *Mitigation:* average over ≥ 12 seeds and use
  large margins (G3/G4 separations are ≈ 5× in the data), so the gates are not borderline.
- **R2 — rank collapse.** At low shots the excited overlap can sink below the 5σ floor and vanish.
  *Mitigation:* gate at a depth (M = 10) and budget where it is resolved in all seeds; record the
  collapse as the finding (§5 G4) rather than fighting it.
- Honest limitation: idealized i.i.d. shot noise, minimal-basis tiny molecules — a characterization,
  not a hardware result.

## 9. Deliverables

- `tests/test_qksd_noise_spec.py` — gates G1–G4 (reusing the existing solver + noise machinery).
- Results summary (the excited/ground fragility ratio and the resolution threshold, with §2/§7
  caveats) in the PR description.
