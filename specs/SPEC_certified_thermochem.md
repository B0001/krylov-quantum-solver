# SPEC: Certified error bars on a relative energy (reaction / dissociation) — no FCI oracle

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G2 (a finite, closing,
two-sided interval).

---

## 1. Goal

A relative energy Δ = E(B) − E(A) (reaction energy, barrier, dissociation/stretch energy) gets a
**certified interval** from Krylov data at each geometry, no FCI, by composing the Temple/Ritz
brackets: Δ ∈ [τ_B − ρ_A, ρ_B − τ_A]. Falsifiable: a depth where the exact in-basis FCI relative
energy falls outside kills it.

## 2. Background and honest framing

The certified arc bounds absolute energies, gaps, and properties; but chemistry's currency is
*relative* energies, which reach experiment most directly. This composes two certified absolute
brackets into a certified difference.

- **What we can claim:** a rigorous, oracle-free interval on a relative energy, validated (against
  FCI, tests only) with zero escapes and closing with depth (H4 stretch to < 0.05 eV), plus the
  structural finding that the error bar **localizes at the strongly-correlated endpoint**.
- **What we cannot claim:** more than the temple bracket supports at the *harder* endpoint (its
  premise sets the two-sidedness; a vacuous lower bound there makes Δ one-sided). Sector-restricted;
  exact statevector (⟨H²⟩ shot cost unmodeled). The certificate contains the **in-basis FCI**
  relative energy — basis-set error vs experiment is separate and uncertified.

## 3. Approach

**Reference (tests only):** the in-basis FCI relative energy (`ground_state_energy` at each
geometry). With E_A ∈ [τ_A, ρ_A] and E_B ∈ [τ_B, ρ_B] from `temple_bounds.krylov_bracket`:

```
Delta_lower = tau_B - rho_A          Delta_upper = rho_B - tau_A          (Delta = E_B - E_A)
```

**Numeric result (H4 symmetric stretch 0.9→2.3 Å):** exact Δ = 8.2255 eV; zero escapes at every M;
two-sided and tight at M=20 (< 0.01 eV); at M=6 the equilibrium bracket (width 0.001 eV) is ~25×
tighter than the stretched (0.025 eV); at intermediate M the stretched Temple lower bound is vacuous
→ Δ_lower = −∞ (one-sided), while Δ_upper stays finite and valid throughout.

## 4. Public interface

Reuses `temple_bounds.krylov_bracket`, `EnergyBracket`.

```
certified_thermochem.certified_relative_energy(mh_a, mh_b, m, solver_a=None, solver_b=None,
                                               e1_a=None, e1_b=None) -> RelEnergyBracket
certified_thermochem.certified_relative_energy_ladder(mh_a, mh_b, dims) -> list[RelEnergyBracket]
RelEnergyBracket(m, delta, delta_lower, delta_upper, width, width_a, width_b)
certified_thermochem (CLI)                                     # H4 stretch M-ladder + inside?
```

## 5. Acceptance criteria (validation gates)

`tests/test_certified_thermochem_spec.py` (test-first).

- **G1 — zero escapes.** The exact FCI relative energy ∈ [Δ_lo, Δ_hi] at every depth (a vacuous side
  counts as ±∞ — still valid one-sided).
- **G2 — two-sided closes (DEFINITION OF DONE).** At M=20 the interval is finite, two-sided,
  < 0.05 eV, and contains the exact ~8.23 eV.
- **G3 — correlated endpoint dominates.** At M=6 the equilibrium bracket width is < 0.1× the
  stretched width — the uncertainty is localized at the hard geometry.
- **G4 — inherits the temple premise; upper always holds.** Some depth has Δ_lower = −∞ (stretched
  Temple vacuous), yet Δ_upper is finite and ≥ the exact Δ at every depth.

## 6. Implementation plan (test-first)

1. `tests/test_certified_thermochem_spec.py` encoding G1–G4 (initially failing — no module).
2. `certified_thermochem.py` composing two `temple_bounds` brackets into a difference.
3. `make gates` (own process; no block2).

## 7. Out of scope

- Geometry optimization / the full PES (pointwise certified energies only).
- Basis-set / relativistic error vs experiment (uncertified, separate).
- Shot-noise statistics on the brackets (exact statevector here).

## 8. Caveats and risks

- **R1 — the harder endpoint sets the premise.** Mitigated by exposing per-endpoint widths and
  reporting one-sided certificates honestly (Δ_lower = −∞), never as a finite bound.
- Reproduction-adjacent (composition of `temple_bounds`); the value is the chemical application
  (relative energies) and the endpoint-localization finding.

## 9. Deliverables

- `certified_thermochem.py` — certified relative-energy interval + CLI.
- `tests/test_certified_thermochem_spec.py` — gates G1–G4.
- `specs/SPEC_certified_thermochem.md` — this spec.
- `specs/BACKLOG.md` — entry with the finding recorded.
