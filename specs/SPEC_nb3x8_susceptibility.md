# SPEC: Finite-T magnetic susceptibility χ(T) of the Nb₃X₈ dimer clusters

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_nb3x8_susceptibility_spec.py`).

---

## 1. Goal

*A thermodynamic study, not a method rung* (precedent: `SPEC_nb3x8_gaps.md`). From the same
exactly-diagonalizable downfolded dimer (`nb3x8_gaps.py`, cRPA parameters of `arXiv:2501.10320`),
compute the finite-temperature molar magnetic susceptibility χ(T) of every Nb₃X₈ parameter set by
an exact Boltzmann trace over the half-filling (N=2) spectrum. The single falsifiable claim: the
exact-spectrum χ(T) reproduces the analytic **Bleaney–Bowers** singlet–triplet law up to a
temperature set by the **charge scale E_s** (the first ionic singlet), not by the exchange J —
and it deviates first, at the lowest reduced temperature, for the iodides, whose spin and charge
scales are least separated (E_s/J ≈ 3). False if the exact χ(T) does not match Bleaney–Bowers in
the spin regime, if the Curie–Weiss θ ≠ −J/4, or if the deviation boundary does not track E_s.

## 2. Background and honest framing

- **Reproduction of ~70-year-old magnetochemistry, applied to new ab-initio parameters.** The
  Bleaney–Bowers dimer law (Bleaney & Bowers, *Proc. R. Soc. A* **214**, 451, 1952; Kahn,
  *Molecular Magnetism*, 1993) and the Van Vleck χ = (g²μ_B²/k_BT)⟨S_z²⟩ trace are textbook. The
  deviation as T → the charge scale is the exact two-site **Hubbard-dimer** thermodynamics
  (Anderson, *Phys. Rev.* **115**, 2, 1959, J = 4t²/U; Carrascal et al., *J. Phys. Condens.
  Matter* **27**, 393001, 2015, `arXiv:1502.05038`) — a known, expected result, not a discovery.
- **Prior Nb₃X₈ susceptibility already exists**, experimentally and theoretically: Sheckelton et
  al. (*Inorg. Chem. Front.* **4**, 481, 2017 — SQUID χ(T), the ~90 K singlet transition in
  Nb₃Cl₈), Haraguchi et al. (*Inorg. Chem.* **56**, 3483, 2017), Grytsiuk/Katsnelson/van
  Loon/Rösner (`arXiv:2305.04854`, Nb₃Cl₈ Mott–Hubbard). This spec does **not** claim first
  observation of dimer magnetism in Nb₃X₈.
- **What we can claim:** the *family-wide* exact-spectrum χ(T), Curie–Weiss θ, and room-T
  effective moments from the ab-initio-downfolded (t, U₀, U_s⊥) — numbers `arXiv:2501.10320` did
  not tabulate — cross-checking the known Nb₃Cl₈ physics, plus the clean, universal
  charge-scale validity boundary T ≲ 0.40·E_s for the pure-spin picture.
- **What we cannot claim:** isolated-single-dimer χ(T) — no inter-dimer coupling, no phonons, no
  structural transition (real Nb₃Cl₈ has a first-order one at ~90 K); density–density
  interactions only; a reference table, not a solid-state prediction.

## 3. Approach

Per parameter set, build the 4-qubit JW Hamiltonian (`dimer_cluster_integrals(**p)
.to_hamiltonian()`), restrict to the N=2 sector via the JW number operator, and evaluate
χ(T) = ⟨S_z,tot²⟩_thermal / T (reduced units: k_B = 1, energies and T in meV, χ in meV⁻¹, g = 2
absorbed into the emu conversion). References: the analytic Bleaney–Bowers form
χ_BB(T) = (2/T)/(3 + e^{J/T}) with J = `odmd_spin.dimer_exchange_analytic`; the Curie–Weiss
θ = −J/4; the exact N=2 spectrum (singlet ground / triplet at J / ionic singlets at E_s).

## 4. Public interface

```
nb3x8_susceptibility.n2_spectrum(U0, t, Us) -> (energies_meV, sz2)   # N=2, rel. to ground
nb3x8_susceptibility.ionic_singlet_energy(U0, t, Us) -> float         # E_s (first level > triplet)
nb3x8_susceptibility.susceptibility(U0, t, Us, T) -> float|ndarray    # reduced χ, meV^-1
nb3x8_susceptibility.bleaney_bowers(J, T) -> float|ndarray            # analytic reduced χ
nb3x8_susceptibility.curie_weiss_theta(J) -> float                    # -J/4 (reduced, meV)
nb3x8_susceptibility.bb_deviation_temperature(U0, t, Us, tol) -> float # T where |Δχ/χ_BB| = tol
nb3x8_susceptibility.EMU_PER_REDUCED  # 0.12931 emu·meV/mol: χ[emu/mol] = EMU_PER_REDUCED·χ_reduced
```

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_susceptibility_spec.py`; all 10 sets in `NB3X8_CLUSTERS` unless noted. Reduced
units; T in meV.

- **G1 — exact spectrum + pure-singlet ground (instant).** For every set: the N=2 triplet sits at
  J (`dimer_exchange_analytic`) within 1e-6 meV; the ground state carries ⟨S_z,tot²⟩ = 0 (< 1e-9,
  a pure singlet — NOT the 0.759 *staggered* moment of `odmd_spin`); and χ·T → 1/3 as T → ∞
  (measured 0.3333 at T = 1e6 meV, < 1e-4).
- **G2 — reproduces Bleaney–Bowers in the spin regime.** For every set, at T = E_s/20 (well below
  the charge scale): `|χ_exact − χ_BB| / χ_BB < 1e-3`. The χ·T Curie constant reaches the
  coupled-pair value 0.5 only as J/T → 0 (at finite T/J, BB itself gives 2/(3+e^{J/T}) < 0.5, e.g.
  0.4935 at T/J = 19), so **the 0.5 check is revised (during implementation)** to fire only for
  members with a clean deep-spin window (J/T = 0.01 at T ≤ E_s/20 both reachable) — the strict
  "= 0.5" over-specified the finite-T value.
- **G3 — Curie–Weiss θ = −J/4.** `curie_weiss_theta(J) == −J/4`; and for the well-separated
  members (E_s/J > 10: the chlorides, bromides, fluorides, HT phases) a C/(T−θ) fit of the exact
  χ over the spin window recovers θ = −J/4 to < 5% (measured θ/J = −0.25). The iodides
  (E_s/J ≈ 3) have **no clean Curie–Weiss window** — recorded in G4, not gated here.
- **G4 — the charge-scale boundary (DEFINITION OF DONE).** The Bleaney–Bowers 5%-deviation
  temperature tracks the charge scale, not J: `T₅% / E_s ∈ [0.38, 0.46]` for all 10 sets
  (measured 0.407–0.444). And the iodides break BB at the lowest reduced temperature —
  `T₅%/J ∈ [1.2, 1.5]` for Nb₃I₈ LT-bulk (measured 1.28), versus > 3 for the bromides and > 10³
  for the fluorides (strictly ordered by E_s/J: I < Br < Cl < F).

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_susceptibility_spec.py` encoding G1–G4 (RED — module missing).
2. `nb3x8_susceptibility.py` — N=2 Boltzmann trace + the analytic references (composition; no new
   solver). Reuses `dimer_cluster_integrals`, `dimer_exchange_analytic`, the JW mapper.
3. `make gates`; record the room-T table (below) in the module `__main__` and the Rösner note.

**Room-T record (T = 26 meV ≈ 300 K, reduced χ in meV⁻¹, μ_eff² = 3Tχ):** fluorides & HT phases
near the free-pair Curie regime (χ ≈ 1.9e-2, μ_eff ≈ 1.22); iodides deep in the singlet
(χ ≈ 6e-6, μ_eff ≈ 0.02). All materials are in the good-Bleaney–Bowers regime at 300 K
(26 meV ≪ every E_s ≥ 774 meV).

## 7. Out of scope

- Inter-dimer coupling / lattice χ(T), phonons, the Nb₃Cl₈ structural transition, anisotropy/
  g-tensor; SI/experimental fitting beyond the stated emu conversion.
- Any quantum/ODMD method (exact trace is trivial at 4 qubits — no advantage; claiming one would
  be dishonest).

## 8. Caveats and risks

- **R1:** the emu conversion assumes g = 2 (spin-only); real Nb₃ trimer moments carry orbital
  contributions. Quote reduced units as primary; emu as indicative.
- Single isolated dimer — the ~90 K Nb₃Cl₈ transition and any inter-dimer physics are absent.

## 9. Deliverables

- `nb3x8_susceptibility.py` (+ `__main__` room-T + χ(T) table); `tests/test_nb3x8_susceptibility_spec.py`.
- `BACKLOG.md` entry; a χ(T) section appended to `notes/nb3x8_cluster_gaps.md` (cross-checking
  Sheckelton/Haraguchi) for the Rösner note.
