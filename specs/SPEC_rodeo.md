# SPEC: The rodeo algorithm filters the spectrum to the ground-state energy

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); `rodeo.py` merged. Ground recovery ≈ 0.13 mHa
(H₂/H₄, K=12); peak sharpens with cycles (H₄ 4.63→0.13 mHa over K=3→12); off-resonance background
decays geometrically (0.35→0.13→0.037); ground-peak height = HF overlap (H₂ 0.987, H₄ 0.936).

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

The rodeo algorithm estimates eigenvalues by a stochastic spectral filter: K cycles of controlled
time evolution with random times, conditioned on the ancilla returning 0, multiply the reference by
`∏_k cos²((H−E)t_k/2)` — a band-pass filter at the target energy E. Claim: the expected survival
probability `P̄(E)` peaks at the Hamiltonian eigenvalues with height equal to the reference overlap
`|⟨E_i|HF⟩|²`; scanning E and taking the dominant low-energy peak recovers the FCI ground-state
energy; the peak sharpens and off-resonance energies are suppressed as the cycle count K grows. The
claim is false if the peaks do not sit at the eigenvalues, if the ground energy is not recovered, or
if off-resonance suppression does not improve with K.

## 2. Background and honest framing

- **Prior art / reference.** Choi, Lee, Bonitati, Qian, Watkins & Lee, *Rodeo algorithm for quantum
  computing*, Phys. Rev. Lett. 127, 040505 (2021), arXiv:2009.04092. A spectral-filtering family,
  distinct from the Krylov / variational / moment / shadow rungs here; it reuses real-time evolution
  `e^{−iHt}` (already validated).
- **Ground truth.** FCI = dense diagonalization of the same qubit Hamiltonian (eigenvalues and HF
  overlaps).
- **What we can claim if gates pass.** The rodeo survival probability filters the spectrum: peaks at
  eigenvalues with height = reference overlap, ground energy recovered to grid resolution, and the
  filter sharpens with K.
- **What we cannot claim (stated up front).** (a) No quantum advantage — exact expected-value
  (infinite-realization) simulation; tiny systems. (b) **Reference-overlap dependence is the honest
  cost:** the ground peak's height is `|⟨HF|E_0⟩|²`, so a poor reference gives a weak ground peak;
  and the K (circuit repetitions) for a target resolution grows with the spectral range — the rodeo
  trades depth for repetitions. (c) Expected survival probability (random-time average); finite
  shot/realization noise is not modeled here.

## 3. Approach

In the eigenbasis the expected (random-time-averaged, `t_k ~ N(0,σ²)`) survival probability is
`P̄(E) = Σ_i |c_i|² [(1 + e^{−(E_i−E)²σ²/2})/2]^K`. Evaluate it from the dense eigenvalues and HF
overlaps; scan E over a low-energy grid; the dominant peak gives the ground energy (+ offset).
Reference: FCI.

## 4. Public interface

```
rodeo.rodeo_survival(eigvals, overlaps, e_target, sigma, n_cycles) -> float   # P_bar(E)
rodeo.reference_spectrum(mh, overlap_tol=1e-8) -> (eigvals, overlaps, offset)
rodeo.rodeo_ground_energy(mh, sigma=2.0, n_cycles=12, n_grid=4000, window=(-0.5,1.5)) -> float
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_rodeo_spec.py` (test-first). Exact statevector; pyscf/qiskit, no block2.

- **G1 — ground recovery (definition of done).** `rodeo_ground_energy` (K=12, σ=2) matches the FCI
  ground-state energy of H₂ and H₄ to `< 3 mHa` (grid-resolution limited; measured ≈ 0.13 mHa).
- **G2 — the filter sharpens with cycles.** On H₄ the ground-energy error at K=12 is smaller than at
  K=3 (measured 0.13 vs 4.63 mHa) — more cycles narrow the peak.
- **G3 — off-resonance suppression.** At an energy between eigenvalues, `P̄` decreases monotonically
  with K (`P̄(K=12) < P̄(K=6) < P̄(K=3)`) — the filter rejects non-eigenvalues as `(<1)^K`.
- **G4 — peak height = reference overlap (the finding).** `P̄` at the ground eigenvalue equals the
  reference overlap `|⟨HF|E_0⟩|²` to `< 0.02` at K=12 (H₂ ≈ 0.987, H₄ ≈ 0.936) — the rodeo signal is
  the reference overlap, so a poor reference gives a weak ground peak.

> Definition of done: **G1**. If the dominant low-energy peak is *not* the ground state (an excited
> state has larger HF overlap), `rodeo_ground_energy` would return that instead — record it and
> scan for the lowest peak above threshold rather than the global max.

## 6. Implementation plan (test-first)

1. Write `tests/test_rodeo_spec.py` encoding G1–G4 (initially failing — module absent).
2. Add `rodeo.py` (expected survival closed form; reference spectrum; ground-energy scan).
3. Iterate to green via `make gates` (own process; pyscf/qiskit, no block2).

## 7. Out of scope

- Finite shot/realization noise and the per-cycle quantum circuit / ancilla simulation (we use the
  expected survival probability).
- Excited-state scanning and multi-peak deconvolution (a natural follow-up — the higher peaks already
  sit at the excited eigenvalues).
- Optimizing σ / the time distribution / adaptive scanning, and the spectral-range cost analysis.

## 8. Caveats and risks

- **R1 — dominant-peak assumption.** `rodeo_ground_energy` takes the global max in the window; if an
  excited state out-weighs the ground state in the reference it would be returned. *Mitigation:* HF
  near equilibrium overlaps the ground state most (G4 checks the overlap); scan the low-energy window.
- **R2 — grid resolution.** The recovered energy is limited by the scan grid. *Mitigation:* a fine
  grid (4000 points) over a narrow window; report the residual as resolution-limited.
- Honest limitation: expected-value simulation, minimal-basis tiny molecules; reference-overlap- and
  repetition-limited.

## 9. Deliverables

- `rodeo.py` — `rodeo_survival`, `reference_spectrum`, `rodeo_ground_energy`.
- `tests/test_rodeo_spec.py` — gates G1–G4.
- Results summary (ground recovery + filter sharpening + the reference-overlap finding, with §2/§7
  caveats) in the PR description.
