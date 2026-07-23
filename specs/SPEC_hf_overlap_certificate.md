# SPEC: Certified HF guiding-state overlap on molecules (SPEC-21 integration)

**Status:** Specced 2026-07-17.
**Depends on:** `certified_overlap` (SPEC-21 core), `certified_gaps.gap_bracket` (self-mode E₁
floor, premise-gated M ≥ 6), `MolecularHamiltonian.hf_state`.

## Claim

`certify_overlap`, fed the repo's own premise-gated Krylov E₁ floor (`gap_bracket.eps1`, self
mode, no oracle), yields a **valid and non-vacuous** lower bound on the Hartree–Fock guiding
state's overlap with the reachable-sector ground state for small molecules at equilibrium —
i.e. the quantity the guided-LH literature *assumes* (γ ≥ 1/poly) becomes a *certified* number
from Krylov data alone.

## Cheap check that could kill it

Dense `eigh` on the qubit Hamiltonian (small systems only) gives the exact reachable-sector
ground state ψ₀ and exact overlap |⟨HF|ψ₀⟩|. If any certified γ_min exceeds the exact overlap,
the claim (and SPEC-21's machinery) dies.

## Sector honesty

The certificate is **sector-restricted**, exactly like QKSD/temple_bounds/certified_gaps:
ψ₀ is the lowest *HF-reachable* eigenstate and eps1 floors the reachable E₁. The Davis–Kahan
decomposition only ever sees reachable components because u = HF defines the sector (its
unreachable components are zero by construction), so the proof carries over unchanged with
"spectrum" read as "reachable spectrum". This is the same scope statement the rest of the
certified suite already makes.

## Gates

- G1 (validity, killable): γ_min ≤ |⟨HF|ψ₀_reachable⟩| exactly (dense reference) on H₂
  (0.74 Å), stretched H₂ (2.0 Å), and H₄ (1.0 Å), for self mode at M ∈ {6, 8, 12} and oracle
  mode. Zero tolerance — one escape kills SPEC-21's machinery.
- G2 (usefulness): at equilibrium H₂ and H₄, the self-mode certificate is **non-vacuous** with
  γ_min ≥ 1/n_qubits — clearing the guided-LH 1/poly(n) threshold the framing needs. (Observed
  M=6: H₂ 0.9936 ≫ 1/4; H₄ 0.4791 ≫ 1/8.)
- G3 (ordering): self-mode γ_min ≤ oracle-mode γ_min ≤ exact overlap (a weaker E₁ floor ⇒ a
  smaller certified δ ⇒ a more conservative certificate; violation would mean a frame bug).
- G4 (premise boundary): self mode at M < 6 raises — the gap machinery's own M ≥ 6 boundary
  (SPEC_temple_bracket) is inherited loudly, not silently ignored.
- G5 (Krylov refinement, the finding): on H₄ the self-mode floor is **non-decreasing in M**
  and converges up toward the oracle floor (observed 0.479 → 0.539 → 0.721 → 0.776 at M=6,8,12
  vs oracle 0.776) — the matvec-priced tightening SPEC-21 §2 predicts, here measured. On H₂ the
  reachable sector is 2-dimensional so self mode equals oracle at every M (σ₁ = 0 exactly); the
  gate asserts that equality too.

## Findings (observed 2026-07-17)

- The certificate is genuinely **useful at equilibrium**: H₂ floor 0.9936 (exact 0.9936, gap
  <1e-4), H₄ floor 0.72 at M=12 (exact 0.968). These are the regimes where a guided-LH pipeline
  would *assume* γ ≥ 1/poly; here it is certified from Krylov data with no oracle.
- **Stretched-bond degradation is real and bounded:** H₂ at 2.0 Å drops the floor to 0.772
  (exact 0.844) — non-vacuous but visibly looser, the orthogonality catastrophe made numeric.

## Out of scope

Degenerate reachable ground spaces (SPEC-21 refuses them; block sin-θ is SPEC-21b), shot noise
on the residual/gap inputs, systems beyond dense-`eigh` validation reach.

## Honest caveats

Exact-statevector only (the shot cost of ⟨H²⟩ is not modeled — same caveat as temple_bounds).
The self-mode premise eps1 ≤ E₁ is unverifiable in-band below M = 6 and inherited as a hard
raise, not re-derived here.
