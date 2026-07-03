# SPEC: Optical absorption & exciton binding via ODMD — and the eigenstate-kick fix

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_odmd_optical_spec.py`).

---

## 1. Goal

The two-particle side of `SPEC_odmd_spectral.md`: kicking the ground state with a *same-sector*
operator (the dipole μ̂ for a molecule, the polarization P = n₁ − n₂ for the dimer clusters)
gives a survival amplitude whose DMD poles are the **bright excited states** and whose weights
are **|⟨E_n|Ô|ψ₀⟩|²** — the optical absorption spectrum, with selection rules emerging as
missing lines. Applied to the Nb₃X₈ dimers this yields numbers the source paper never reported:
the **optical gaps** (analytically pinnable: the dimer's only bright state is the odd-singlet at
exactly U₀) and the **exciton binding energies** Δ_charge − Δ_optical, which collapse from
≈ Us in the atomic limit (F) to 0.26·Us for the iodide — the exciton unbinds with hopping.

## 2. Background, the found defect, and honest framing

- **Found while probing (fixed and gated here):** `odmd_spectral.reference_signal` **hangs** when
  the kicked reference is an exact eigenstate — precisely what P|ψ₀⟩ *is* on the inversion-
  symmetric dimer (P is odd, ψ₀ even, and the singlet sector has exactly one odd state). The
  reachable width is then 0, the `max(width, 1e-12)` guard set τ = π·10¹², and `expm_multiply`
  span ~10¹⁷ substeps (two probe processes burned 35 CPU-min before being killed). Fix: a
  degenerate-reference short-circuit — the signal is exactly constant, returned without
  evolution. First real use of the #8 module found its edge case; the gate makes it permanent.
- Prior art: dipole-autocorrelation absorption spectra via real-time propagation are standard;
  `SPEC_qksd_properties.md` already validated transition dipoles on these systems (the bright
  HeH⁺ |μ| ≈ 0.85 is G1's cross-pin: weight 0.85² = 0.7224).
- **What we can claim if gates pass:** machine-exact absorption lines vs dense FCI; the Nb₃X₈
  optical gaps and exciton bindings, pinned both by dense ED and by the analytic odd-singlet
  energy `ω_opt = U₀ − E₀`, `E₀ = (U₀+Us)/2 − √(((U₀−Us)/2)² + 4t²)`; and two gated material
  trends (binding/Us and total oscillator weight, both monotone in the halide series).
- **What we cannot claim:** exact statevector signals; isolated-dimer optics (no band
  broadening, no non-density-density terms — the capstone caveats); density-density models have
  no true dipole, P is the standard stand-in; exciton "binding" here is the cluster's
  charge-vs-neutral gap difference, not a solid-state exciton dispersion.

## 3. Approach

`absorption_lines(mh, kick_op, reference)` reuses `reference_signal`/`lines_from_signal` verbatim
(they are operator-agnostic) with ω = pole − E₀(reference sector); the elastic line at ω = 0
carries ⟨Ô⟩² when nonzero (Rayleigh line — physical). Dimer helpers: `dimer_polarization()`
(JW-mapped n₁ − n₂), `dimer_optical_gap` (the analytic odd-singlet formula), and
`dimer_exciton_binding = nb3x8_device_gap.exact_gap − dimer_optical_gap`. References: dense
diagonalization; the analytic formula; the capstone's sector-FCI charge gaps.

## 4. Public interface

```
odmd_optical.absorption_lines(mh, kick_op, reference=None, n=24, amp_floor=1e-6)
    -> (omegas, weights)                       # relative to the reference-sector ground
odmd_optical.dimer_polarization() -> sparse    # P = n_1 - n_2, JW block ordering, 4 qubits
odmd_optical.dimer_optical_gap(U0, t, Us) -> float      # analytic: U0 - E0_even
odmd_optical.dimer_exciton_binding(U0, t, Us) -> float  # exact charge gap - optical gap
```

Changed: `odmd_spectral.reference_signal` — degenerate-reference short-circuit (returns the
constant signal; no behavior change for any non-degenerate reference — the pinned
`SPEC_odmd_spectral` gates re-run green).

## 5. Acceptance criteria (validation gates)

`tests/test_odmd_optical_spec.py`.

- **G1 — absorption lines are exact (HeH⁺, multi-line).** Every visible line of the μ_z-kicked
  exact ground state: pole and weight vs dense FCI < 1e-10 (measured ≤ 8e-14); the transition
  line's weight equals the `SPEC_qksd_properties` bright transition, 0.7224 ± 1e-3; the elastic
  line equals the permanent dipole².
- **G2 — the eigenstate kick terminates and is exact (the fix's gate).** P|ψ₀⟩ on the Nb₃I₈
  dimer: `absorption_lines` returns (pre-fix: hangs), with exactly ONE line, at pole = U₀ to
  < 1e-9 relative, weight = ‖P|ψ₀⟩‖². 
- **G3 — optical gaps and exciton binding (DEFINITION OF DONE).** All four LT-bulk materials:
  the ODMD optical gap matches the analytic `dimer_optical_gap` < 1e-6 meV; binding(F) within
  2% of Us (measured 0.986 — the atomic limit); binding/Us strictly decreasing F → Cl → Br → I
  (measured 0.986, 0.486, 0.358, 0.263) — the exciton unbinds with hopping.
- **G4 — selection rules and the brightness ladder.** Each dimer's spectrum holds exactly one
  bright line while its sector holds ≥ 4 levels (the odd-singlet selection rule as a gate); the
  total oscillator weight ‖P|ψ₀⟩‖² increases strictly F → Cl → Br → I (measured 1.1e-4, 0.22,
  0.44, 0.96 — a 4-orders-of-magnitude polarizability ladder).

## 6. Implementation plan (test-first)

1. `tests/test_odmd_optical_spec.py` encoding G1–G4 (RED — module missing; G2 additionally
   hangs against the unfixed `reference_signal`).
2. Fix `reference_signal`; add `odmd_optical.py` (composition + the analytic formula).
3. `make gates`; re-run `tests/test_odmd_spectral_spec.py` (pinned) after the fix.

## 7. Out of scope

- Shot noise / circuits on optical signals; solid-state excitons (dispersion, screening);
  non-density-density dimer terms; molecular systems beyond HeH⁺-scale validation.

## 8. Caveats and risks

- **R1:** P is a stand-in for the dipole in a model without geometry — trends are physical,
  absolute intensities are model-defined.
- The F dimer's bright line carries 1.1e-4 of weight — near-dark; any sampled measurement of it
  would face the visibility law head-on (compose with `SPEC_odmd_uq`/`SPEC_device_odmd` to cost
  it — a follow-up).

## 9. Deliverables

- `odmd_optical.py`; the `odmd_spectral.py` fix; `tests/test_odmd_optical_spec.py` (G1–G4).
- `BACKLOG.md` entry with the exciton-binding table.
