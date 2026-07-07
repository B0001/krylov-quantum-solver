# SPEC: Certified error bars on a molecular property (dipole) — no FCI oracle

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G2 (a useful, closing bar).

---

## 1. Goal

The ground-state dipole moment (and any ⟨ψ₀|A|ψ₀⟩ of a bounded Hermitian A) gets a **certified
interval** [μ ± half_width] from Krylov data alone (no FCI), via the Davis–Kahan eigenvector bound
fed by the certified GAP lower bound of `certified_gaps`. Falsifiable: a single depth where the exact
FCI dipole falls outside a finite certified interval kills it.

## 2. Background and honest framing

The certified arc bounds energies (`temple_bounds`), gaps and self-verification (`certified_gaps`,
`gap_selfcheck`); properties — the observables spectroscopy reports — were still point estimates
(`qksd_properties`). This extends the arc to properties, reusing the gap certificate directly.

- **What we can claim:** a rigorous, oracle-free dipole interval, validated (against FCI, in tests
  only) to contain the exact dipole with zero escapes and to close with depth (LiH −1.818 ± 0.065
  a.u. at M=24 vs exact −1.817); and the honest structural result that a **property certificate
  inherits the gap certificate** beneath it.
- **What we cannot claim:** more than the gap certificate supports. The sinθ bound rests on the same
  premise as `certified_gaps` (M ≥ 6; oracle-free-checkable via `gap_selfcheck`). Sector-restricted
  ground state; exact statevector (⟨H²⟩/⟨A²⟩ shot cost unmodeled); full orbital space (the
  `build_dipole_operators` convention). On a system small enough to diagonalize, certification is
  pointless.

## 3. Approach

For the ground Ritz state ψ₀ (Rayleigh quotient ρ₀, residual σ₀ = ‖(H−ρ₀)ψ₀‖), the
eigenvector-perturbation theorem gives sinθ ≤ σ₀/(E₁−ρ₀), and **E₁−ρ₀ ≥ Δ_lo**, the certified gap
lower bound. So sinθ ≤ s := σ₀/Δ_lo. The **sharp** first-order property bound (fluctuation, not
operator norm):

```
| ⟨ψ₀|A|ψ₀⟩ − ⟨exact|A|exact⟩ |  ≤  2 σ_A(ψ₀) sinθ + W_A sin²θ  ≤  2 σ_A s + W_A s²  =: half_width
```

with σ_A(ψ₀) = √(⟨A²⟩−⟨A⟩²) the dipole fluctuation and W_A = λ_max−λ_min the spectral width. Using
σ_A instead of ‖A‖ is decisive (LiH: σ_A ≈ 1.1 vs ‖μ_z‖ ≈ 6.9 → ~6× tighter). **Reference (tests
only):** the exact HF-reachable ground-state dipole by dense diagonalization.

**Numeric result:** zero escapes on HeH⁺ and LiH at all depths; LiH closes to −1.818 ± 0.065 a.u. at
M=24; half_width is finite iff s < 1, so it is vacuous exactly at the depths where Δ_lo is weak
(LiH M=8–16) and sharp where Δ_lo is healthy (M ≥ 20).

## 4. Public interface

Reuses `certified_gaps.gap_bracket`, `temple_bounds._mean_and_variance`, `build_dipole_operators`.

```
certified_dipole.certified_dipole(mh, a_sparse, m, width=None, solver=None, e1=None) -> CertifiedDipole
certified_dipole.certified_dipole_ladder(mh, a_sparse, dims, solver=None) -> list[CertifiedDipole]
certified_dipole.spectral_width(a_sparse) -> float
CertifiedDipole(m, mu, half_width, sin_theta_bound, sigma_A, gap_lower, finite)
certified_dipole (CLI)                                       # per-molecule M-ladder + inside?
```

## 5. Acceptance criteria (validation gates)

`tests/test_certified_dipole_spec.py` (test-first). Reference = HF-reachable ground (correct
particle-number sector — not the global lowest eigenvector, which for a charged species differs).

- **G1 — zero escapes.** Exact dipole ∈ [μ ± half_width] at every depth; the deepest interval is
  finite (not passing only via vacuous ∞) and still contains it.
- **G2 — closes & useful (DEFINITION OF DONE).** LiH at M=24 is finite with half_width < 0.15 a.u.
  containing exact −1.817; HeH⁺ certified to < 1e-2 a.u.
- **G3 — fluctuation beats operator norm.** half_width < 0.3 × the naive 2‖A‖ bound; σ_A < 0.25‖A‖.
- **G4 — inherits the gap certificate (the boundary).** half_width finite ⟺ s < 1 ⟺ σ₀ < Δ_lo; LiH
  shows both regimes (vacuous where Δ_lo weak, sharp where healthy).

## 6. Implementation plan (test-first)

1. `tests/test_certified_dipole_spec.py` encoding G1–G4 (initially failing — no module).
2. `certified_dipole.py` composing the residual + certified gap lower bound + Davis–Kahan.
3. `make gates` (own process; no block2).

## 7. Out of scope

- Transition dipoles / excited-state properties (needs excited-eigenvector bounds).
- Shot-noise statistics on σ₀, σ_A (exact statevector here).
- Active-space dipole convention (full orbital space only, per `build_dipole_operators`).

## 8. Caveats and risks

- **R1 — inherits the gap premise.** Mitigated by pairing with `gap_selfcheck` (oracle-free
  trustworthiness of Δ_lo) and reporting vacuous intervals honestly as ∞, never as a tight bound.
- The reference must be the reachable ground (a bug caught in implementation: the global lowest
  eigenvector of a charged species is a different sector — G1 pins this).

## 9. Deliverables

- `certified_dipole.py` — certified property interval + CLI.
- `tests/test_certified_dipole_spec.py` — gates G1–G4.
- `specs/SPEC_certified_dipole.md` — this spec.
- `specs/BACKLOG.md` — entry with the finding recorded.
