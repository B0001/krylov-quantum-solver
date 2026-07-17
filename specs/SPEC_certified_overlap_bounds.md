# SPEC-21: Certified Guiding-State Overlap Bounds

**Status:** Core implemented 2026-07-17 (Davis–Kahan γ_min, Temple shared-provenance floor,
invariants I1–I4, binding 1000-trial validity gate). `krylov_refine` remains a stub —
Lanczos-chained refinement and the degenerate-ground-space block form are SPEC-21b candidates.
**Depends on:** existing certified gap-input machinery (Temple/Lehmann path); Krylov/Lanczos core
**arXiv provenance verified:** 2026-07-16

## 1. Problem

The guided local Hamiltonian literature establishes that ground-state energy estimation becomes BQP-complete when a classical guiding state u with overlap γ = |⟨u|ψ₀⟩| ≥ 1/poly(n) is supplied (Gharibian–Le Gall arXiv:2111.09079; Cade–Folkertsma–Niesen–Weggemans arXiv:2207.10097; Gharibian–Hayakawa–Le Gall–Morimae arXiv:2207.10250), with active 2025 work on physically motivated guiding states (arXiv:2509.25815) and classical/dequantized regimes (arXiv:2411.16163, arXiv:2509.25829).

Every one of these guarantees is **conditional on γ** — and in practice γ is asserted, not certified. The orthogonality catastrophe makes the assertion increasingly suspect at scale. Nobody in this pipeline certifies the condition their guarantee stands on.

Our solver already certifies energies. This spec extends the same machinery to certify **overlap**: given a guiding state and a certified gap input, output a rigorous lower bound γ_min ≤ |⟨u|ψ₀⟩|.

## 2. Mathematical core

Let H be Hermitian on the fixed-basis Fock space, u a normalized trial state.

Computable quantities (no conditions):
- Rayleigh quotient: λ_u = ⟨u|H|u⟩
- Residual norm: r = ‖(H − λ_u)u‖

Conditional bound (Davis–Kahan sin-θ, simple ground state):
Let δ be a certified lower bound on dist(λ_u, spec(H) \ {E₀}), derived from the solver's existing certified two-sided machinery (e.g., certified E₀ upper bound + certified E₁ lower bound). If r < δ:

    sin θ(u, ψ₀) ≤ r / δ
    ⇒ γ_min = sqrt(1 − r²/δ²) ≤ |⟨u|ψ₀⟩|

Complementary energy bound from the same inputs (Temple):

    E₀ ≥ λ_u − r² / (β − λ_u)    for any certified β ≤ E₁ with λ_u < β

so one (u, gap-input) pair yields BOTH a certified overlap floor and a certified lower energy bound — they share provenance and must be reported together.

Krylov refinement: running Lanczos from u produces Ritz vectors with monotonically non-increasing residuals; each iterate yields a (possibly) tighter γ_min for the *Ritz* vector, chained back to u via computable inner products. This makes the certificate improvable at the cost of matvecs, matching the solver's existing convergence loop.

## 3. Deliverable semantics

Input: (H handle, u, gap_certificate)
Output: OverlapCertificate {
  gamma_min, lambda_u, residual_norm,
  gap_certificate_id,             // provenance chain, mandatory
  conditional: true,              // always — never claim unconditional overlap
  bqp_threshold_note              // whether gamma_min ≥ 1/poly for stated n, per guided-LH framing
}

## 4. Invariants (non-bypassable, loud-failure — house style)

- I1: No γ_min is ever emitted without a valid gap_certificate. Missing/expired certificate ⇒ raise. Never a warning, never a default gap.
- I2: If r ≥ δ, the bound is vacuous: emit an explicit VACUOUS result (γ_min = 0 with reason), never a fabricated positive number.
- I3: Every OverlapCertificate embeds the id of the gap certificate it is conditional on; serialization without it must fail.
- I4: gamma_min must never exceed 1 − ε_machine; clamp-and-flag is forbidden — if numerics produce > 1, raise (it indicates an upstream error).

## 5. Implementation sketch (for scaffolding — see HANDOFF T2)

    certified_overlap/
      __init__.py
      rayleigh.py        # rayleigh_quotient(H, u)          — implement now
      residual.py        # residual_norm(H, u, lambda_u)    — implement now
      davis_kahan.py     # gamma_min(r, delta) + I2 logic   — stub
      temple.py          # shared-provenance Temple bound   — stub (wraps existing lower-bound path; do not duplicate)
      certificate.py     # OverlapCertificate + I1/I3/I4    — stub with invariant raises wired
      krylov_refine.py   # Lanczos-chained refinement       — stub

Property tests (bind future implementation):
- Validity: on random Hermitian H (dim ≤ 200) with exact eigendecomposition, gamma_min ≤ |⟨u|ψ₀⟩| for 1000 random u. Zero tolerance.
- I1: call without gap certificate ⇒ raises.
- I2: engineered r ≥ δ case ⇒ VACUOUS, not positive.

## 6. Why this is an invention and not a feature

It inverts the direction of trust in the guided-quantum-algorithms stack: instead of a quantum-advantage claim resting on an unverified overlap assumption, the classical certifier produces the assumption's proof. Any downstream QPE/guided-LH pipeline can cite an OverlapCertificate the way it would cite an error bar. To our knowledge no existing tool emits certified γ floors with provenance-chained gap inputs; the nearest literature (arXiv:2509.25815) motivates guiding states physically but does not certify them.

## 7. Open questions for review

- Should δ derivation live here or stay exclusively in the existing gap machinery (recommended: exclusively there; this module only consumes)?
- Degenerate/near-degenerate ground spaces: extend via subspace sin-θ (Davis–Kahan block form) or explicitly refuse in v1? (Recommended: refuse loudly in v1, spec the block form as SPEC-21b.)
- MPS-format guiding states (DMRG imports, per arXiv:2509.25815): matvec-free residual evaluation feasibility.
