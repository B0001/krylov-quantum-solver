# SPEC-21b: Certified subspace overlap for (near-)degenerate ground clusters

**Status:** Specced 2026-07-23. The block-form follow-up SPEC-21 flagged twice
(`SPEC_certified_overlap_bounds.md` §7; `SPEC_hf_overlap_certificate.md` out-of-scope).
**Depends on:** `certified_overlap` (SPEC-21 core: `rayleigh_quotient`, `residual_norm`,
`OverlapCertificate` invariants I1–I4).

## Problem

SPEC-21 certifies `γ_min ≤ |⟨u|ψ₀⟩|` for a **simple** ground state. When the ground level is
degenerate — or merely near-degenerate — that quantity is the wrong target: within a degenerate
eigenspace the individual eigenvector ψ₀ is basis-arbitrary, so `|⟨u|ψ₀⟩|` is not even
well-defined, and for a *near*-degenerate cluster the single-vector certificate's separation
δ = E₁ − λ_u collapses to the tiny intra-cluster spacing and the bound goes **vacuous** (r ≥ δ).

Near-degenerate ground manifolds are not a corner case here: they are the strongly-multireference
regime (stretched bonds, square H₄, the Nb₃X₈ singlet–triplet near-crossing) — exactly where a
guiding state's quality matters most.

The physically meaningful, basis-independent quantity is the overlap with the **ground
eigenspace** S = span{ψ₀,…,ψ_{d−1}} of the lowest d levels:

    ‖P_S u‖ = cos θ(u, S),   P_S = orthogonal projector onto S.

## Mathematical core (block Davis–Kahan sin-θ)

Same residual as SPEC-21, r = ‖(H − λ_u)u‖, λ_u = ⟨u|H|u⟩. Expand u in the eigenbasis; split
the sum at the cluster boundary d:

    r² = Σᵢ |cᵢ|² (λᵢ − λ_u)² ≥ Σ_{i≥d} |cᵢ|² (λᵢ − λ_u)²
       ≥ δ² Σ_{i≥d} |cᵢ|² = δ² (1 − ‖P_S u‖²) = δ² sin²θ(u, S),

for any certified δ with 0 < δ ≤ dist(λ_u, {λᵢ : i ≥ d}). Hence

    sin θ(u, S) ≤ r / δ   ⇒   γ_min = √(1 − r²/δ²) ≤ ‖P_S u‖.

When λ_u < E_d, dist(λ_u, {λᵢ : i ≥ d}) = E_d − λ_u, so any certified floor β ≤ E_d gives the
valid (conservative) separation δ = β − λ_u. **d = 1 recovers SPEC-21 exactly** (S = {ψ₀},
β floors E₁): the single-vector certificate is the special case.

The crucial difference from v1: δ now measures to the first level **above the cluster**, not the
first excited level. Bracketing the whole near-degenerate cluster (large δ) is what rescues a
certificate the single-vector form (tiny intra-cluster δ) throws away.

## Deliverable

`ClusterGapCertificate {e_above_floor, cluster_size d, certificate_id, source}` — the consumed
gap input, flooring E_d (the first eigenvalue **above** the size-d cluster). Derivation stays in
the gap machinery; this module only consumes (SPEC-21 §7 recommendation, unchanged).

`certify_subspace_overlap(H, u, cluster_gap_certificate, n_qubits=None) -> OverlapCertificate`
with `cluster_size = d` recorded on the certificate; γ_min bounds ‖P_S u‖.

## Invariants (inherited, non-bypassable, loud)

- I1: no gap certificate ⇒ raise (never warn/default).
- I2: r ≥ δ or λ_u ≥ β ⇒ explicit VACUOUS (γ_min = 0 with reason), never a fabricated positive.
- I4: γ_min > 1 − ε_machine ⇒ raise.
- Ib: cluster_size d must be an integer ≥ 1 ⇒ raise otherwise (a subspace of dimension < 1 is
  meaningless; non-integer d is a caller bug).

## Gates

- G1 (validity, killable, zero-tol): on random Hermitian H (dim ≤ 200) with exact eigen-
  decomposition and oracle β = E_d, γ_min ≤ ‖P_S u‖ for 1000 random trial vectors across
  cluster sizes d ∈ {1, 2, 3}. One escape kills the block bound.
- G2 (v1 consistency): at d = 1, `certify_subspace_overlap` reproduces SPEC-21
  `certify_overlap`'s γ_min to machine precision (same object, same number).
- G3 (**the finding** — vacuous-vs-useful contrast): construct a Hamiltonian with a 2-fold
  near-degenerate ground cluster (spacing ε) gapped by Δ ≫ ε from the rest, and a trial u
  spread across the two near-degenerate states. The v1 single-vector certificate is VACUOUS
  (δ ≈ ε ⇒ r ≥ δ) while the d = 2 block certificate is **non-vacuous** and valid
  (γ_min ≤ ‖P_S u‖). This is the whole point, asserted.
- G4 (Ib boundary): d = 0, d = −1, and non-integer d each raise.
- G5 (I1/I2): missing certificate raises; λ_u ≥ β yields VACUOUS, not an error.

## Out of scope

Molecular reachable-cluster demonstration (square-H₄ / Nb₃X₈ singlet–triplet near-crossing via
the HF-overlap path) — a clean follow-up (`SPEC_hf_overlap_subspace`) once this core lands.
Shot noise on r/β. Automatic cluster-size *detection* (d is a caller input here; picking d from
a certified spectral-cluster structure is its own hypothesis). Degenerate excited-state targets.

## Honest caveats

Exact-statevector only (⟨H²⟩ shot cost unmodeled, inherited from temple_bounds). The certificate
bounds overlap with the *chosen* size-d cluster; choosing d smaller than the true near-degeneracy
re-collapses δ and returns (correctly) a vacuous or weak bound — the certificate never lies, but
a bad d wastes it. Picking d is the caller's physics judgment, not certified here.
