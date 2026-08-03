# SPEC: The certified energy bracket under shot noise — a probabilistic certificate

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G2 (coverage is
N-independent — shots do not buy back the guarantee).

---

## 1. Goal

Under i.i.d. shot noise (standard errors set by the Hamiltonian 1-norms λ_H, λ_{H²}), the certified
energy bracket's guarantee **breaks** — and the break is *structural*, not a finite-sample effect:
raw coverage of E₀ is ~0.4 and **N-independent**, and the variational upper bound holds only ~0.5 of
the time. Inflating by z·standard-error restores coverage ≥ 0.9, with the inflated half-width scaling
as λ_H/√N — shots buy tightness, not coverage. Falsifiable: if raw coverage climbed toward 1 with N,
or inflation failed to reach 0.9, the claim is wrong.

## 2. Background and honest framing

The certified arc (`temple_bounds`, `certified_gaps`, `certified_dipole`, `certified_thermochem`) is
exact-statevector. This module asks the hardware question — does the certificate survive finite
sampling? — and re-grounds the whole arc at once (every rung rests on ⟨H⟩, ⟨H²⟩).

- **What we can claim:** a quantified, MC-validated account of how sampling breaks the certificate,
  the surprising N-independence of the break, and a conservative inflation that restores coverage
  with a λ/√N shot-cost law.
- **What we cannot claim:** rigor under noise. The inflated bracket is a *probabilistic* certificate
  under an idealized i.i.d. Gaussian model (real grouped-Pauli measurement differs by O(1)); it
  cannot see systematic (Trotter/basis) bias. An oracle gap feeds Temple (noisy ε₁ only worsens the
  break, so the finding is conservative). ⟨H²⟩ carries the larger 1-norm — the Temple lower bound is
  the noise-expensive side.

## 3. Approach

For the exact Ritz ground state (⟨H⟩ = ρ₀, var₀), draw N-shot estimates ρ₀ + 𝒩(0, (λ_H/√N)²) and
⟨H²⟩ + 𝒩(0, (λ_{H²}/√N)²); form the noisy bracket [τ₀, ρ₀]; Monte-Carlo the coverage of the exact
reachable E₀. Inflate: [τ₀(ρ₀−z·se, var+z·se) , ρ₀ + z·se].

**Numeric result (z=2, 4000 trials):** raw coverage ≈ 0.39 (H2) / 0.40 (H4), var-upper ≈ 0.49, both
**identical across N = 10⁴, 10⁶, 10⁸**; inflated coverage ≈ 0.98; inflated half-width 54→5.4→0.54
mHa (H2) and 205→20.5→2.05 mHa (H4) as N×100 (exact 1/√N); λ_{H²} = 62.9 ≫ λ_H = 10.3 (H4).

## 4. Public interface

Reuses `temple_bounds.mean_and_variance`, the Krylov solver.

```
certified_noise.hamiltonian_one_norms(mh) -> (lambda_H, lambda_H2)
certified_noise.reachable_E0_E1(mh) -> (E0, E1)                 # dense reference
certified_noise.certified_half_width(lambda_H, shots, z=2.0) -> float
certified_noise.shot_noise_coverage(mh, m, shots, trials=4000, z=2.0, seed=0, solver=None) -> dict
certified_noise (CLI)                                          # coverage table per system/budget
```

## 5. Acceptance criteria (validation gates)

`tests/test_certified_noise_spec.py` (test-first).

- **G1 — sampling breaks the certificate.** At converged depth, raw coverage < 0.65 and the
  variational upper bound holds < 0.65 (a coin flip), for H2 and H4.
- **G2 — shots do not buy coverage (DEFINITION OF DONE).** Raw coverage moves < 0.12 across
  N = 10⁴…10⁸, and stays broken (< 0.65) — the guarantee is structural.
- **G3 — inflation restores coverage + width law.** z·se inflation gives coverage ≥ 0.9 at every
  budget; the inflated half-width shrinks 10× per 100× shots (1/√N).
- **G4 — ⟨H²⟩ is the noise-expensive side + determinism.** λ_{H²} > λ_H for both systems; coverage
  is deterministic under a fixed seed.

## 6. Implementation plan (test-first)

1. `tests/test_certified_noise_spec.py` encoding G1–G4 (initially failing — no module).
2. `certified_noise.py` — 1-norms + MC coverage of the noisy bracket + inflation.
3. `make gates` (own process; no block2).

## 7. Out of scope

- Rigorous (non-probabilistic) bounds under noise (impossible from finite samples).
- Grouped-Pauli covariances / a full measurement-allocation model (i.i.d. λ-1-norm idealization).
- Systematic (Trotter/basis) bias (the model sees only statistical error).

## 8. Caveats and risks

- **R1 — idealized noise.** Stated up front; real measurement differs by O(1), but the qualitative
  findings (break, N-independence, inflation-restores, λ/√N) are model-robust.
- The inflation constant z is a design choice (conservative); coverage ≈ 0.98 at z=2 vs target 0.9.

## 9. Deliverables

- `certified_noise.py` — noisy-bracket coverage study + shot-cost width + CLI.
- `tests/test_certified_noise_spec.py` — gates G1–G4.
- `specs/SPEC_certified_noise.md` — this spec.
- `specs/BACKLOG.md` — entry with the finding recorded.
