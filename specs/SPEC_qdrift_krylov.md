# SPEC: qDRIFT real-time Krylov, done correctly

**Status:** IMPLEMENTED — all gates green (`tests/test_qdrift_krylov_spec.py`, ~30 s).

## 1. Goal

Real-time Krylov diagonalization whose evolutions e^{−ikΔtH} are approximated by **qDRIFT**
(random single-Pauli rotations) converges to the exact-sector FCI as the number of qDRIFT samples
per step L grows, **and is never below the variational floor** — provided each random realization
supplies a whole Krylov basis and the H, S matrices are averaged over realizations. Falsified if
error does not fall with L or any run lands below E_FCI on a system where Fock-space and N-electron
minima coincide.

## 2. Background and honest framing

Triggered by an externally drafted "SqDRIFT Krylov" agent-skill template (Nov 2026) that repeated
every failure mode `REFACTOR_PLAN.md` §2 recorded for this repo's original core: vacuum start
(0-electron sector), a "qDRIFT step" that samples a term then evolves under the full H (so all
Krylov vectors are phase copies, S rank 1), an unthresholded `eigh(H, S)` (crashes on that S), a
full-space Hamiltonian compared against a CAS(2,2) energy, and a sample-count claim O(λ²/ε²)
instead of O(λ²t²/ε). On H₂O CAS(2,2) and CAS(4,4) the template raises `LinAlgError`.

- **Can claim:** a correct qDRIFT-Krylov path with a *proven* (Gram-matrix) variational floor,
  measured convergence vs exact-sector FCI, and a measured boundary in λ.
- **Cannot claim:** novelty (qDRIFT: Campbell, PRL 123, 070503 (2019); real-time Krylov: Klymko
  et al., PRX Quantum 3, 020323 (2022); SqDRIFT: Piccinelli et al., arXiv:2508.02578 — this is a
  reproduction of known ingredients); any shot-noise or device-noise robustness (exact
  statevector of the channel); any quantum advantage (≤ 8 qubits, minimal basis).

## 3. Approach

One realization = one random walk of single-Pauli rotations e^{−iτ sgn(c_j) P_j}, j ~ |c_j|/λ,
τ = λΔt/L, reading off |φ_k⟩ every L rotations from the Hartree–Fock reference. Average
S = ⟨V†V⟩, H = ⟨V†HV⟩ over R realizations, then `solve_generalized_eig` (thresholded canonical
orthogonalization, reused). Because the average equals W†(I_R⊗H)W for the stacked W, every Ritz
value ≥ min spec(H) on the reachable space. **Reference:** exact (N_α, N_β)-sector ground energy by
dense diagonalization (`sector_ground_energy`).

## 4. Public interface

```
qdrift_krylov.qdrift_krylov(mh, krylov_dim=6, samples_per_step=400, realizations=20, dt=None,
                            threshold=1e-6, postselect=False, seed=None, psi0=None)
    -> QDriftKrylovResult(energy, rank, H, S, max_leakage, dt, tau, lam)
qdrift_krylov.qdrift_sample_count(lam, t, eps) -> int        # ceil(2 λ² t² / ε)
qdrift_krylov.check_reference(mh, psi0)                      # raises on out-of-sector reference
qdrift_krylov.sector_mask(mh) / sector_ground_energy(mh)     # reference helpers
```

## 5. Acceptance criteria (validation gates)

H₂O/STO-3G, M=6, R=20, seeds {0,1,2}, L ∈ {50, 400, 2000}; medians over seeds.

- **G1 — convergence (definition of done).** Median error strictly decreasing in L; at L=2000,
  CAS(2,2) < 3e-5 Ha (measured 8.8e-6), CAS(4,4) < 1.6 mHa (measured 0.89 mHa).
- **G2 — variational floor.** Every run, both systems, with and without postselection:
  E ≥ E_sector − 1e-9. (Holds in all 36 runs.)
- **G3 — λ boundary.** λ(CAS(4,4)) = 11.2 > 4·λ(CAS(2,2)) = 2.44·4; at L=400 CAS(2,2) is within
  chemical accuracy (3.7e-5) while CAS(4,4) is not (3.4 mHa). Samples needed grow with λ.
- **G4 — particle-number leakage (the finding).** CAS(2,2) leaks 0; CAS(4,4) leaks up to 15 % of
  the norm at L=50, 0.5 % at L=2000. With H + μN̂ (μ = −2 Ha, Fock-space minimum 2.4 Ha below the
  N-electron sector), un-postselected runs land **10–42 mHa below** sector FCI; postselected runs
  stay above (+3.2–3.4 mHa). Postselection also lowers the physical-case error (~10–25 %).
- **G5 — template defects detected.** Vacuum reference raises `ValueError`; the template's
  full-H phase-copy basis has rank-1 S and unthresholded `scipy.linalg.eigh(H, S)` raises.
- **G6 — cost law.** `qdrift_sample_count` quadruples when t doubles, doubles when ε halves.

> Revision note (v1 → v1.1, during the probe): v1 planned "E ≥ E_FCI" as the floor. Reality: a
> single Pauli string does not commute with N̂, so the unconditional floor is the **Fock-space**
> minimum. G2 holds here only because Fock and sector minima coincide for neutral H₂O; G4 builds
> the case where they don't and shows the bound is violated without postselection.

## 7. Out of scope

Shot noise / Hadamard-test estimation of the averaged matrices; Trotter–qDRIFT hybrids; partial
randomization; circuit compilation and depth costing; sample-based (SKQD/SqDRIFT) diagonalization.

## 8. Caveats and risks

- The floor proof is exact-arithmetic; shot noise on H, S breaks it (as in `SPEC_qksd_noise`).
- G1 thresholds are seed-dependent; three seeds, 20 realizations. Spread across seeds at L=2000 is
  ~±15 % of the error.
- Δt = π/spectral width uses dense diagonalization — fine at ≤ 8 qubits, not scalable.
- Postselection discards the leaked norm; on hardware that costs shots in proportion to leakage.

## 9. Deliverables

- `qdrift_krylov.py`, `tests/test_qdrift_krylov_spec.py`, this spec, BACKLOG entry.
- Corrected agent skill: `skills/krylov-quantum-solver-standard/` (SKILL.md + standalone example).
