# SPEC: Quantum imaginary-time evolution reaches FCI; a truncated operator domain does not

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); `qite.py` merged. Exact ITE monotone + variational,
→FCI (H₂ β=4, LiH β=15); full-domain QITE (McLachlan update) reproduces exact ITE and reaches FCI on
H₂; step error shrinks with Δτ. Finding: weight-≤2 domain stalls at Hartree–Fock (+20.5 mHa) — the
correlation lives in a weight-4 operator; only the full domain reaches FCI.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

QITE realises imaginary-time evolution `|ψ(β)⟩ = e^{−βH}|ψ⟩/‖·‖` — which projects variationally onto
the ground state — by replacing each non-unitary step `e^{−ΔτH}` with a unitary `e^{−iΔτÂ}`, Â a
Hermitian sum of Paulis over a *domain*, found from the McLachlan linear system `S a = b`. Claim:
exact imaginary-time evolution is monotone-decreasing, variational, and converges to FCI; the
**full-domain** QITE unitary update reproduces it (validating the algorithm's equations) and reaches
FCI; and a **truncated (low-weight) domain** cannot — QITE's accuracy is set by whether the domain
spans the operators that build the correlation. The claim is false if exact ITE is not variational,
if full-domain QITE misses FCI, or if the truncated-domain failure is not observed.

## 2. Background and honest framing

- **Prior art / reference.** Motta et al., *Determining eigenstates and thermal states … using
  quantum imaginary time evolution*, Nature Physics 16, 205 (2020), arXiv:1901.07653. A distinct
  family from the Krylov / variational / moment / shadow / filter rungs here.
- **Ground truth.** Exact imaginary-time evolution (`expm_multiply`) and FCI (dense diagonalization)
  of the same qubit Hamiltonian.
- **What we can claim if gates pass.** The QITE McLachlan update is implemented correctly (full
  domain reproduces exact ITE and reaches FCI), the step error vanishes as Δτ→0, and the
  operator-domain locality is quantified (a low-weight domain stalls).
- **What we cannot claim (stated up front).** (a) No quantum advantage — exact statevector, dense
  Pauli domain, tiny systems. (b) **The full domain is 4ⁿ operators**, so only H₂ (n=4) is exercised;
  the scalable version uses *local* domains, and its accuracy/cost is exactly the locality trade-off
  G4 exposes. (c) First-order (single-step) McLachlan update; higher-order/expm step-error control is
  not studied.

## 3. Approach

Exact ITE: normalized `expm_multiply(-βH)` on the HF reference → E(β). QITE step: `S_IJ =
Re⟨ψ|σ_Iσ_J|ψ⟩`, `b_I = Im⟨ψ|σ_I H|ψ⟩` (least-squares match of `-iÂ|ψ⟩` to `-(H-⟨H⟩)|ψ⟩`), solve
`S a = b`, apply `e^{-iΔτÂ}` (renormalize). Domain = all Paulis (exact) or weight-capped (local).
Reference: exact ITE and FCI.

## 4. Public interface

```
qite.pauli_operators(n_qubits, max_weight=None) -> np.ndarray      # (n_ops, 2^n, 2^n)
qite.exact_imaginary_time(mh, betas) -> list[float]                # the variational reference
qite.qite_evolve(mh, dtau, n_steps, operators) -> list[float]      # energy after each step
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_qite_spec.py` (test-first). Exact statevector; pyscf/qiskit, no block2.

- **G1 — exact ITE is variational and converges.** For H₂ and LiH, E(β) is monotone non-increasing
  and stays `≥ E_FCI - 1e-9`; it reaches FCI at a gap-set β — H₂ `< 1e-4` at β=4, LiH's small gap
  needs β=15 (`< 1e-3`). (Very large β overflows `expm_multiply`, so β is kept moderate — a numerical
  caveat, not a physics one.)
- **G2 — full-domain QITE reproduces exact ITE and reaches FCI (definition of done).** On H₂ with the
  full Pauli domain (Δτ=0.1), QITE tracks exact ITE within the step error and lands within `1e-4` of
  FCI at β=4 — the McLachlan update equations are correct.
- **G3 — step error vanishes with Δτ.** On H₂, QITE at Δτ=0.05 is closer to exact ITE (at matched β)
  than QITE at Δτ=0.2 — the O(Δτ²) single-step error shrinks.
- **G4 — domain locality (the finding).** On H₂ a weight-≤2 Pauli domain **stalls at Hartree–Fock**
  (`E_QITE(β=4) - E_FCI ≈ +20.5 mHa`, unmoved from ⟨H⟩), because the correlation lives in a weight-4
  double-excitation operator; the full domain reaches FCI. QITE's accuracy is set by the domain.

> Definition of done: **G2 + G4** — a correct full-domain update *and* the locality boundary that
> says why the scalable (local-domain) version is limited. If a larger system shows the correlation
> captured by a low-weight domain, record where the locality assumption holds.

## 6. Implementation plan (test-first)

1. Write `tests/test_qite_spec.py` encoding G1–G4 (initially failing — module absent).
2. Add `qite.py` (Pauli domain; exact ITE reference; McLachlan-update QITE step).
3. Iterate to green via `make gates` (own process; pyscf/qiskit, no block2).

## 7. Out of scope

- Local-domain QITE on larger systems, unitary-Trotter circuit synthesis, and the measurement cost of
  `⟨σ_Iσ_J⟩`/`⟨σ_I H⟩` (we use exact statevector expectation values).
- Thermal states / METTS (the finite-temperature side of arXiv:1901.07653).
- Higher-order step-error control and adaptive Δτ.

## 8. Caveats and risks

- **R1 — S is singular.** The full-domain overlap matrix is rank-deficient. *Mitigation:* solve by
  least squares with a small cutoff (`rcond`), which selects a valid descent direction.
- **R2 — 4ⁿ domain.** Only H₂ is tractable full-domain. *Mitigation:* scope full-domain gates to H₂;
  the reference ITE gates use any size.
- Honest limitation: exact statevector, minimal-basis H₂, dense Pauli domain — a correctness +
  locality study, not a scalable-implementation result.

## 9. Deliverables

- `qite.py` — `pauli_operators`, `exact_imaginary_time`, `qite_evolve`.
- `tests/test_qite_spec.py` — gates G1–G4.
- Results summary (exact-ITE convergence, full-domain QITE = exact ITE = FCI, the Δτ step error, and
  the weight-≤2 stall-at-HF locality finding, with §2/§7 caveats) in the PR description.
