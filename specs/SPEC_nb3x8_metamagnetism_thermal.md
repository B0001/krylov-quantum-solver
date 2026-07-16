# SPEC: Nb3X8 metamagnetism at finite temperature — the crossover width law

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G2.

---

## 1. Goal

`SPEC_nb3x8_metamagnetism` gates the T=0 field-driven singlet→triplet ground-state crossing
(sharp step at h_c = J). This spec adds temperature: the full N=2-sector Boltzmann trace under
the same Zeeman term gives the finite-T magnetization M(h,T), which must (a) reduce to the sharp
T=0 step as T→0, (b) have its zero-field slope dM/dh|_{h=0} match the **independently-implemented**
`nb3x8_susceptibility.susceptibility(T)` (a real cross-check between two separately-built modules,
not the same code called twice), and (c) show the T=0 step smoothing into a crossover whose
10%–90% width in h is a **closed-form multiple of T** in the regime T ≲ 0.1·J: width = 2·ln(9)·T
(the width of a logistic/Fermi function's 10–90 interval) — **and** that this law measurably
breaks down by T ~ 0.2–0.3·J, gated as a finding rather than hidden by a loose tolerance (see §5
G4). Falsifiable: a slope mismatch with `susceptibility`, a width that doesn't track linearly with
T in the safe regime, or a family-dependent prefactor would kill the claim.

## 2. Background and honest framing

The Zeeman-augmented dimer's field-independent-apart-from-a-rigid-shift structure (established in
`SPEC_nb3x8_metamagnetism`) means the low-field, near-crossing physics of the N=2 sector reduces,
in the regime T ≪ E_s − J (charge scale far above the crossing — already established across the
whole family in `SPEC_nb3x8_susceptibility`/`SPEC_nb3x8_thermo`), to an effective two-level system:
the h-independent singlet and the h-dependent polarized branch. A two-level thermal population
gives a logistic magnetization curve in (J−h)/T, whose 10–90% width in h is exactly 2·ln(9)·T —
this spec checks whether the full 6-state trace (which also contains the h-independent Sz=0
triplet member and the far-away ionic singlets) actually reduces to that clean two-level number,
or whether the other levels measurably contaminate it.

- **What we can claim:** a genuinely new, temperature-resolved observable for this family; a
  real independent cross-check between two separately-built modules (this spec would fail if
  either module had a sign or scale bug the other didn't share); a closed-form, material-
  independent crossover width law, checked (not assumed) against the full multi-level trace.
- **What we cannot claim:** an experimentally observed magnetization curve (same caveat as
  `SPEC_nb3x8_metamagnetism` — this is a prediction, and the T=0 crossing itself already sits at
  Tesla-scale fields beyond routine lab reach for Cl/Br/I). The 2·ln(9)·T law is measured to hold
  to < 1e-4 relative for T ≲ 0.1·J and to have already broken (> 0.1% deviation) by T = 0.2·J,
  worsening past 1% by T = 0.3·J — **found while probing, not assumed**: an initial version of
  this gate at 1e-3 tolerance up to T = 0.2·J failed by a hair (0.10% vs the 0.10% bound), which is
  what prompted checking the trend explicitly rather than loosening the tolerance to paper over
  it. Not claimed at T comparable to or above E_s − J, where the ionic singlets would necessarily
  intrude further.

## 3. Approach

**References:**
1. T→0 limit: `nb3x8_metamagnetism.magnetization` (the already-gated sharp step).
2. Zero-field slope: `nb3x8_susceptibility.susceptibility(U0, t, Us, T)` — built independently
   (different file, different derivation route — the Van Vleck `<Sz²>/T` trace over the
   *zero-field* spectrum) from this spec's finite-difference `dM/dh` of the *field-augmented*
   trace. Agreement is a real consistency check, not a tautology.
2. Closed form: the two-level logistic width, 2·ln(9)·T (elementary; not a repo primitive).

## 4. Public interface

Reuses `nb3x8_metamagnetism.field_spectrum` (generalizes `zeeman_ground_state` to the full
spectrum), `nb3x8_susceptibility.susceptibility`. New module `nb3x8_metamagnetism_thermal.py`:

```
nb3x8_metamagnetism_thermal.magnetization_thermal(U0, t, Us, h, T) -> float   # <Sz> Boltzmann trace
nb3x8_metamagnetism_thermal.crossover_width(U0, t, Us, T, lo=1e-4, hi=None) -> float  # h(M=0.9)-h(M=0.1)
nb3x8_metamagnetism_thermal (CLI __main__)                                    # family table + finding
```

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_metamagnetism_thermal_spec.py` (test-first).

- **G1 — T→0 recovers the sharp step.** At T = 1e-4·J, `magnetization_thermal(h=0.9J) < 1e-3` and
  `magnetization_thermal(h=1.1J) > 1 - 1e-3`, matching `nb3x8_metamagnetism.magnetization`'s own
  step, for Cl/Br/I.
- **G2 — susceptibility cross-check (DEFINITION OF DONE).** Central finite difference
  `dM/dh|_{h=0}` at T = 0.1·J matches `nb3x8_susceptibility.susceptibility(T)` to
  `< 1e-5` relative, for Cl/Br/I — two independently-built modules agreeing on an unshared
  derivation route.
- **G3 — monotonicity.** `magnetization_thermal(h, T)` is non-decreasing in h across a grid
  spanning the crossing, at two different T, for Cl/Br/I.
- **G4 — the crossover width law, and its breakdown.** `crossover_width(T) / T` equals `2*ln(9)`
  (≈4.394) to `< 1e-4` relative for T ∈ {0.01, 0.05, 0.1}·J, for Cl/Br/I — the ratio is both
  T-independent (linear scaling) and material-independent in this regime, despite the full trace
  carrying 4 more levels than the effective two-level picture assumes. **The recorded boundary,
  gated not just noted:** at T = 0.2·J the deviation has already exceeded 0.1%, at T = 0.3·J it
  exceeds 1%, and the deviation strictly grows between the two — the two-level picture measurably
  fails once T stops being small compared to J.

> G2 is the definition of done: it is the one gate that could fail from a bug in *either*
> independently-written module, not just this one.

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_metamagnetism_thermal_spec.py` encoding G1–G4 (initially failing — no module).
2. `nb3x8_metamagnetism_thermal.py` composing `field_spectrum` + a Boltzmann trace (same style as
   `nb3x8_susceptibility.n2_spectrum`'s thermal average) plus a bisection-based `crossover_width`.
3. `make gates` (own process; no block2/qiskit-aer conflict).

## 7. Out of scope

- T comparable to or above E_s − J, where the ionic singlets would intrude on the two-level
  picture (would need the full charge-scale boundary machinery from `SPEC_nb3x8_susceptibility`).
- Deriving the 2·ln(9)·T law analytically from the full 6-state Hamiltonian (checked numerically
  against the elementary two-level closed form, not re-derived from the correlated dimer).
- Nb₃F₈ (excluded, same reason as `SPEC_nb3x8_metamagnetism`: J below the model's noise floor).
- Any experimental comparison (see §2).

## 8. Caveats and risks

- **R1:** the 2·ln(9) law is an emergent two-level approximation, not an exact identity of the
  6-state system — G4's tolerance (1e-3 relative) is empirically set from the tested range; it is
  not claimed to hold at arbitrarily large T/J.
- Inherits R1–R3 of `SPEC_nb3x8_metamagnetism` (g=2 assumption, density-density only, isolated
  single dimer).

## 9. Deliverables

- `nb3x8_metamagnetism.py` — extended with `field_spectrum` (full-spectrum generalization of
  `zeeman_ground_state`, which now composes it).
- `nb3x8_metamagnetism_thermal.py` — new module.
- `tests/test_nb3x8_metamagnetism_thermal_spec.py` — gates G1–G4.
- `specs/BACKLOG.md` — recorded finding.
