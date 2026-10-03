# SPEC: Nb3X8 magnetocaloric effect — a quantified negative

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G1 (machinery); G3 carries
the verdict.

---

## 1. Goal

`SPEC_nb3x8_thermo` §7 and `SPEC_nb3x8_metamagnetism` §7 both name magnetocaloric S(T,B) as a
follow-up. One composition closes both: feed `nb3x8_metamagnetism.field_spectrum`'s Zeeman-
augmented spectrum into `nb3x8_thermo.entropy`'s Boltzmann-trace formula. Two claims, checked
independently:

- **(i) Machinery.** The Maxwell relation (∂S/∂B)_T = (∂M/∂T)_B is an exact thermodynamic identity
  (Clairaut's theorem on the mixed partial of F(T,h) = −T ln Z(T,h)); computing its two sides from
  two independently-derived Boltzmann traces should agree to machine precision.
- **(ii) Verdict.** |ΔS_M| peaks at the critical field B_c (already gated in
  `SPEC_nb3x8_metamagnetism`, e.g. 572 T for Nb₃Cl₈) where the full spin-sector entropy becomes
  thermally accessible, but at any laboratory-achievable field the response is exponentially
  suppressed — **Nb₃X₈ dimers are quantitatively useless magnetocalorics below megagauss fields.**

Falsifiable: the machinery dies if the two derivatives disagree by more than 1e-6 relative
anywhere in the spin window swept. The verdict — the physically hoped-for outcome, so its death
would be a real finding, not a footnote — dies if any halide reaches ≥10% of R ln2 per formula
unit at B ≤ 100 T and T ≥ 2 K.

## 2. Background and honest framing

**This is a bounding result by construction** (precedent: `SPEC_senseforge` §3–4) — it reports an
upper bound on a laboratory-accessible effect, not a measured curve. Composes two already-gated
primitives with no new physics: `field_spectrum` (`SPEC_nb3x8_metamagnetism`, direct
diagonalization of H0 − h·Ŝz) and the reduced Boltzmann-trace entropy formula
(`SPEC_nb3x8_thermo`, S = ln Z + ⟨E−E₀⟩/T, k_B = 1).

- **What we can claim:** a real cross-implementation check of the Maxwell relation (a bug in
  either module's sign/scale convention would show up here, since the two sides are computed by
  genuinely different code paths — entropy's partition-function route vs magnetization's ⟨Sz⟩
  route); a quantified, unit-converted, sourced comparison against GGG showing the isolated-dimer
  magnetocaloric response is negligible under any field a pulsed magnet can reach.
- **What we cannot claim:** an experimentally observed magnetocaloric curve (same caveat as the
  parent specs); anything at or near B_c, which remains beyond even destructive pulsed-field
  records (`SPEC_nb3x8_metamagnetism` G4); a GGG comparison against a specific measured ΔS_M(T,B)
  curve — see the sourcing caveat below.
- **GGG sourcing caveat (read before trusting the ratio):** this session's web-fetch tools
  produced inconsistent, uncitable numbers for GGG's magnetocaloric response on repeated attempts
  (one fetch reported a GGG molar mass of 644.37 g/mol, ~36% off the formula-weight arithmetic
  below — a sign of fabrication, not a real quote). Rather than cite an unverifiable number, this
  spec uses GGG's **rigorous, textbook, code-derived full paramagnetic-entropy ceiling**: Gd³⁺ is
  4f⁷ with L=0, S=7/2 (Hund's-rule ground term ⁸S₇/₂, an orbital singlet isolated from the first
  excited term by several eV — precisely why GGG is chosen as a magnetic refrigerant), so its
  total magnetic entropy content is 3·R·ln(8) per mole (3 Gd per formula unit) — an upper bound on
  any real ΔS_M GGG could ever show. Comparing against this ceiling is **conservative in the
  direction of the negative verdict**: if Nb₃X₈ is tiny against GGG's theoretical maximum, it is
  tinier still against any real, smaller, achieved GGG number.
- Inherits g=2, density-density-only, isolated-single-dimer from the parent specs — not re-derived
  here (see their own §8 caveats).

## 3. Approach

**References:**
1. Machinery: the Maxwell relation itself (an exact identity of any well-defined free energy),
   checked via two independently-coded finite differences — `(∂S/∂h)_T` from this spec's
   `entropy_field`, `(∂M/∂T)_h` from the already-gated `nb3x8_metamagnetism_thermal.magnetization_thermal`.
2. Verdict: R ln2 per formula unit (the parent specs' own plateau reference, `SPEC_nb3x8_thermo`
   G2/G3) as the kill threshold; GGG's 3·R·ln(8) ceiling (derived, §2) as the physical-units yardstick.
3. Molar masses: `pyscf.data.elements.MASSES` (natural-abundance standard atomic weights), not
   hand-copied constants — code-derived, per the repo's "the code produces it" rule.

## 4. Public interface

Reuses `nb3x8_metamagnetism.field_spectrum`, `nb3x8_metamagnetism.G_MU_B`,
`nb3x8_metamagnetism_thermal.magnetization_thermal`, `nb3x8_magnetometry.MEV_PER_K`,
`odmd_spin.dimer_exchange_analytic`. New module `nb3x8_magnetocaloric.py`:

```
nb3x8_magnetocaloric.entropy_field(U0, t, Us, h, T) -> float | ndarray   # S(T,h), reduced, per dimer
nb3x8_magnetocaloric.delta_entropy_isothermal(U0, t, Us, T, h) -> float # S(T,h)-S(T,0), reduced
nb3x8_magnetocaloric.ds_dh_numeric(U0, t, Us, h, T, dh) -> float        # (dS/dh)_T, entropy route
nb3x8_magnetocaloric.dm_dT_numeric(U0, t, Us, h, T, dT) -> float        # (dM/dT)_h, magnetization route
nb3x8_magnetocaloric.delta_s_m_fraction_of_r_ln2(U0, t, Us, T_K, B_T) -> float
nb3x8_magnetocaloric.delta_s_m_j_per_kg_k(U0, t, Us, T_K, B_T, mass) -> float
nb3x8_magnetocaloric.scan_max_delta_s(U0, t, Us, b_max_tesla=100, t_min_kelvin=2, ...) -> (max, B*, T*)
nb3x8_magnetocaloric.nb3x8_formula_mass(halide_symbol) -> float          # g/mol
nb3x8_magnetocaloric.NB3X8_FORMULA_MASS_G_PER_MOL                       # dict
nb3x8_magnetocaloric.GGG_MOLAR_MASS_G_PER_MOL                           # float, g/mol
nb3x8_magnetocaloric.GGG_ENTROPY_CEILING_J_PER_KG_K                     # float
nb3x8_magnetocaloric (CLI __main__)                                     # family table + finding
```

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_magnetocaloric_spec.py` (test-first).

- **G1 — Maxwell relation (DEFINITION OF DONE, machinery).** Across a documented spin window
  (T/J ∈ {0.1, 0.15, 0.2, 0.3, 0.5, 0.8}, h/J ∈ {−0.3, 0.05, 0.3, 0.6, 0.9, 1.1, 1.4, 1.8}), central
  finite differences `ds_dh_numeric` and `dm_dT_numeric` agree to < 1e-6 relative, for Cl/Br/I.
  **Scope note:** this window avoids T/J < 0.1 and h/J very close to the exact crossing/edge
  points, where the true derivative passes through zero and floating-point/finite-difference noise
  dominates the relative comparison — including **h/J = 0 exactly**, where both sides vanish
  identically by the h → −h symmetry (a 0/0, not a disagreement) (verified empirically — worst case
  in the chosen window is ~1.9e-7; points outside it were seen to reach ~2% relative or a literal
  0/0, an ill-conditioning artifact, not a machinery bug).
- **G2 — the peak is really near B_c (supports the verdict's premise).** For each halide: (a)
  `|Delta S_M(T,B)|` is non-decreasing in B over a grid spanning [0, 100] T at representative T
  (confirms the search-grid max legitimately sits at the boundary, not an interior artifact); (b)
  `|Delta S_M|` at `B = B_c` (low T) exceeds `|Delta S_M|` at B=100T by at least 10x — a cheap,
  directly-computed check that the search window is genuinely far from the interesting physics.
- **G3 — the 10%-of-R-ln2/f.u. kill criterion (VERDICT — take seriously if it fails).** For each of
  Nb₃Cl₈/Br₈/I₈: `scan_max_delta_s` over B ∈ [0, 100] T, T ∈ [2 K, 5·(J/k_B)] gives
  `delta_s_m_fraction_of_r_ln2 < 0.10` at the found maximum.
- **G4 — sourced, unit-matched GGG comparison.** `NB3X8_FORMULA_MASS_G_PER_MOL` and
  `GGG_MOLAR_MASS_G_PER_MOL` match hand-checked reference values (guards silent drift in the
  underlying atomic-weight table); `GGG_ENTROPY_CEILING_J_PER_KG_K` matches its closed form
  `3*R*ln(8)/M_GGG`; Nb₃Cl₈'s `delta_s_m_j_per_kg_k` at its G3 maximum is < 1% of
  `GGG_ENTROPY_CEILING_J_PER_KG_K`.

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_magnetocaloric_spec.py` encoding G1–G4 (initially failing — no module).
2. `nb3x8_magnetocaloric.py` composing `field_spectrum` + the thermo entropy formula + the GGG
   reference constants.
3. `make gates` (own process; no block2 — same qiskit-nature + numpy footprint as
   `nb3x8_metamagnetism_thermal.py`).

## 7. Out of scope

- Any GGG comparison against a specific measured ΔS_M(T,B) curve (sourcing caveat, §2) — the
  ceiling used here is a theoretical maximum, not a lab-achieved number.
- T beyond ~5·(J/k_B) or fields beyond B_c (would need the charge-scale intrusion already scoped
  out by the parent specs).
- Any renormalization from the isolated cluster to the real lattice (same caveat as
  `SPEC_nb3x8_metamagnetism` §7).

## 8. Caveats and risks

- **R1 — G1's window is empirically chosen, not derived.** The Maxwell relation is an exact
  identity everywhere; the 1e-6 gate is a numerical-implementation check whose finite-difference
  comparison is ill-conditioned near points where the true derivative vanishes (documented in G1).
- **R2 — GGG ceiling, not GGG's realized ΔS_M.** See the sourcing caveat in §2; the comparison is
  conservative (generous to GGG), not an apples-to-apples experimental match.
- Inherits R1–R3 of `SPEC_nb3x8_metamagnetism` (g=2, density-density only, isolated dimer) and the
  scope limits of `SPEC_nb3x8_thermo` (no phonons/lattice, no structural transition).

## 9. Deliverables

- `nb3x8_magnetocaloric.py` — new module (S(T,B), Maxwell cross-check, GGG comparison, CLI).
- `tests/test_nb3x8_magnetocaloric_spec.py` — gates G1–G4.
- `specs/BACKLOG.md` — entry marked closed with the finding recorded.
