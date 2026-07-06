# SPEC: Nb3X8 magnetic heat capacity & entropy — the Schottky anomaly and the R ln 2 plateau

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G3 (the finding).

---

## 1. Goal

The exact-spectrum magnetic heat capacity C_m(T) and entropy S_m(T) of the downfolded Nb₃X₈ dimer
show (a) a Schottky anomaly pinned at T ≈ 0.352 J by the analytic two-level (singlet/triplet)
result — a J-scale fingerprint whose ratio to the χ(T) peak is universal — and (b) a localized-moment
entropy plateau R ln 4/dimer (R ln 2/cluster) that is clean **only** when the charge scale E_s ≫ J,
degrading monotonically with E_s/J and vanishing for the iodide. Falsifiable: a Schottky peak off
0.352 J, or a clean plateau surviving for the iodide, would break it.

## 2. Background and honest framing

Completes the thermodynamic triad — `nb3x8_susceptibility` gave χ(T); this adds C and S from the same
exact N=2 Boltzmann trace (cRPA parameters of `arXiv:2501.10320`). The two-level Schottky anomaly and
the R ln(2S+1) magnetic-entropy plateau are textbook (Gopal, *Specific Heats at Low Temperatures*,
1966); the contribution is the family-wide ab-initio numbers and the charge-scale plateau boundary.

- **What we can claim:** family-wide C_m/S_m tables, the universal Schottky/χ peak-ratio cross-tie,
  and the E_s/J-controlled plateau boundary — numbers the source paper did not report.
- **What we cannot claim:** novelty of the physics (reproduction); isolated single dimer, no
  inter-dimer coupling, no phonon/lattice heat capacity, no structural transition (real Nb₃Cl₈ has a
  first-order one at ~90 K); density-density only; Nb₃F₈ (J ~ 0.05 meV, below neglected terms) has no
  resolvable spin feature. A reference table, not a solid-state heat-capacity prediction.

## 3. Approach

**Reference = analytic two-level Schottky/entropy limits** + the exact N=2 spectrum. Reduced units
(k_B = 1): C_m = Var(E)/T², S_m = ln Z + ⟨E−E₀⟩/T. The two-level (g₀=1, g₁=3) Schottky peak
(T_pk/Δ = 0.3515) and the R ln 4 plateau are the checkable limits.

**Numeric result (LT-bulk):** Schottky peak at T_pk/J = 0.351/0.351/0.357 (Cl/Br/I), i.e. 270/486/
1018 K; C-peak/χ-peak = 0.564/0.564/0.580 (≈ 0.3515/0.625, universal); entropy plateau S/ln4 =
0.987/1.016/1.171 with flatness (min |dS/dlnT|) = 0.061/0.207/0.253 — strictly worsening as E_s/J =
16.9/8.1/3.1 shrinks.

## 4. Public interface

Reuses `n2_spectrum`, `ionic_singlet_energy`, `dimer_exchange_analytic`, `chi_max_temperature`.

```
nb3x8_thermo.heat_capacity(U0, t, Us, T) -> float | ndarray     # reduced C_m(T)
nb3x8_thermo.entropy(U0, t, Us, T) -> float | ndarray           # reduced S_m(T), per dimer
nb3x8_thermo.two_level_schottky_ratio(g_ratio=3.0) -> float     # analytic T_pk/Delta ~ 0.3515
nb3x8_thermo.schottky_peak_temperature(U0, t, Us) -> float      # meV
nb3x8_thermo.entropy_plateau(U0, t, Us) -> (S_plateau, flatness)
nb3x8_thermo (CLI)                                              # family table + finding
```

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_thermo_spec.py` (test-first).

- **G1 — Schottky peak pinned.** Analytic (1,3) peak at T/J ≈ 0.3515; each dimer's C_m peaks there
  within 3%; C ≥ 0; third-law S(T→0) = 0.
- **G2 — localized-moment plateau.** Nb₃Cl₈ (E_s/J ≈ 17) has a flat (min |dS/dlnT| < 0.10) entropy
  plateau within 3% of ln 4 (= 2·ln 2, per cluster); sector entropy saturates at ln 6.
- **G3 — plateau cleanliness is charge-scale-set (DEFINITION OF DONE, the finding).** Both the
  flatness metric and the deviation from ln 4 increase strictly Cl → Br → I as E_s/J shrinks; the
  iodide has no clean plateau (dev > 10%), the chloride does (dev < 5%). Same E_s boundary as χ(T).
- **G4 — cross-tie to χ + scope.** C-peak/χ-peak ≈ 0.562 across the family (within 5%,
  material-independent); Nb₃F₈ (J ~ 0) excluded.

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_thermo_spec.py` encoding G1–G4 (initially failing — no module).
2. `nb3x8_thermo.py` composing the validated `n2_spectrum` trace into C and S.
3. `make gates` (own process; no block2).

## 7. Out of scope

- Lattice/phonon heat capacity, the cooperative structural transition, inter-dimer coupling.
- Nb₃F₈ spin thermodynamics (J below neglected terms).
- A magnetocaloric / field-dependent S(T,B) treatment (a follow-up).

## 8. Caveats and risks

- **R1 — plateau metric.** `flatness` = min |dS/dlnT| on [2J, E_s/3] is a heuristic; the gate uses
  the strict Cl<Br<I ordering plus absolute dev bounds, robust to the exact window.
- Reproduction of known physics; the value is the ab-initio table + the E_s/J boundary, reported
  honestly, not a claim of new thermodynamics.

## 9. Deliverables

- `nb3x8_thermo.py` — C_m/S_m + Schottky/plateau tools + CLI.
- `tests/test_nb3x8_thermo_spec.py` — gates G1–G4.
- `specs/SPEC_nb3x8_thermo.md` — this spec.
- `specs/BACKLOG.md` — entry with the finding recorded.
