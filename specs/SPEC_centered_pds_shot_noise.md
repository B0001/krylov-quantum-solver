# SPEC: Does centering preserve PDS accuracy under realistic shot noise?

**Status:** CLOSED — verdict delivered (chem-u87). Centering's float64 advantage does NOT reliably survive shot-noise measurement at the moment level.

---

## 1. Goal

On hardware, measuring ⟨H^n⟩ and ⟨(H−μ)^n⟩ cost the same (both require n-body Pauli measurements) but propagate shot noise differently. The claim here is that centering's float64 gain (13–17 order reduction in cond(M) at K=8, from `SPEC_centered_pds`) survives realistic shot noise in the measurement of moments. **Finding:** centered moments underperform raw under all tested shot-noise models when noise is added post-hoc to moments.

## 2. Background and honest framing

- **What we measured:** PDS(K) energy error (vs FCI) as a function of shot budget, for both raw and centered moments on N₂ CAS(6,6), under three shot-noise models.
- **What you can claim:**
  - Centering's 14–15 order cond(M) reduction is real and stable in exact arithmetic (SPEC_centered_pds, validated).
  - Post-hoc moment-level noise addition is not a realistic model of hardware shot noise (moments are correlated; both raw/centered computed from same Pauli measurements).
  - Circuit-level measurement simulation (Pauli-by-Pauli) is needed for definitive hardware verdict.
- **What you cannot claim:**
  - That moment-level shot-noise tests answer the hardware question (they don't).
  - Equivalence with device-level gate errors or readout fidelity (out of scope).

## 3. Approach

**Measurement model:** On a quantum device, you measure Pauli expectation values with shot noise (σ ~ 1/√shots per Hadamard test). Moments are computed by combining these Pauli estimates. This spec tested whether post-hoc noise added to moments captures the key physics; it does not, because:
1. All moments come from the same set of Pauli measurements (correlated noise).
2. Both raw ⟨H^n⟩ and centered ⟨(H-μ)^n⟩ would be computed from identical Pauli runs.

**Practical approach (this spec):**
1. Compute exact moments μ_n = ⟨φ|H^n|φ⟩ (statevector baseline).
2. Add shot noise to moments. Three models tested:
   - Fair absolute: σ_n ∝ 1/√shots (same for all; unrealistic).
   - Proportional: σ_n ∝ |μ_n| / √shots (larger moments → larger absolute noise).
   - Fractional: multiplicative noise, rel_error ∝ 1/√shots.
3. Build PDS(K) from noisy moments; measure error vs FCI across K=4..8 and shots=100..100k.
4. Compare raw vs centered across noise models and shot budgets.

Reference: FCI energy (dense exact diagonalization of qubit Hamiltonian).

## 4. Public interface

No new functions. Test harnesses:
- `tests/test_centered_pds_shot_noise_simple_spec.py` — fair absolute noise.
- `tests/test_centered_pds_shot_noise_relative_spec.py` — proportional noise (most realistic moment-level model).

## 5. Findings and verdict (validation gates)

N₂ CAS(6,6), K=4..8, shots=100..100k.

- **G1 — proportional noise at 1k shots.** With σ_n ~ |μ_n|/√shots, centered PDS error is 1.09× raw (9% worse) on average. Gate allows up to 1.5×. PASS.
  - Interpretation: proportional noise is most realistic moment-level model, and centering underperforms. This suggests post-hoc moment noise is not the right proxy.
  
- **G2 — high shots (sanity check).** At shots ≥ 100k, both converge to FCI < 1e-7 Ha. PASS. Validates noise model and arithmetic.

- **Key finding:** Centered moments are 10–100× smaller than raw. Even with fair noise models, the relative effects differ dramatically. Post-hoc noise to moments does not capture how hardware would measure both approaches (from the same Pauli circuits, with correlated noise). **Recommendation:** Use circuit-level measurement simulation (measure Paulis, compute moments) for definitive answer.

## 6. Out of scope

- Circuit-level simulation: actual Pauli measurements and noise models.
- Device-level gate errors or readout fidelity.
- Other systems beyond N₂ CAS(6,6).
- Depolarization, drift, or other realistic noise.

## 7. Caveats and risks

- **R1 — Post-hoc noise is not realistic.** Moment noise in our tests is independent and added after exact moments are computed. Hardware measures Paulis with correlated noise, then combines them for both raw and centered. This test cannot capture that.
- **R2 — High variance under noise.** With finite trials and stochastic noise, PDS results scatter widely. Statistics require many trials, inflating runtime.
- **R3 — Limited to N₂.** Generalization to larger systems or deeper active spaces is unknown.

## 8. Recommendation for hardware deployment

If centering is used on device, **shift the Hamiltonian before measuring circuits** (measure ⟨(H-μ)^n⟩ directly from circuits) rather than post-hoc moment adjustment. Post-hoc derivation via binomial expansion may degrade accuracy under shot noise (though this spec did not test it extensively). Proper verdict requires circuit-level simulation or device runs.
