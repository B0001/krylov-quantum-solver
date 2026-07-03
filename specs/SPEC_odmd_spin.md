# SPEC: Spin spectroscopy — the interlayer exchange J of the Nb₃X₈ dimers

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_odmd_spin_spec.py`).

---

## 1. Goal

Complete the response trilogy (`SPEC_odmd_spectral` = charge, `SPEC_odmd_optical` = optical):
kick the dimer ground state with the **staggered magnetization** S₁ᶻ − S₂ᶻ — the operator that
cannot move charge — and the survival amplitude exposes exactly the state the polarization
leaves dark: the m=0 triplet, at ω = J, the singlet–triplet splitting. This yields the
**interlayer magnetic exchange constants of the Nb₃X₈ family** from the paper's own cRPA
parameters (J = 0.051 / 66.2 / 119.1 / 245.9 meV for F/Cl/Br/I), each pinned by the closed form
`J = √(((U₀−Us)/2)² + 4t²) − (U₀−Us)/2`, plus a falsifiable physics claim: the textbook
Heisenberg superexchange `4t²/(U₀−Us)` fails progressively across the family, by **46.5% for
the iodide** — Nb₃I₈'s interlayer dimer is beyond the Heisenberg regime.

## 2. Background and honest framing

- Machinery is 100% reuse: `odmd_optical.absorption_lines` with a spin kick; Sᶻ|ψ₀⟩ on the
  dimer is again an *exact eigenstate* (Sᶻ annihilates the ionic components), so G1 also
  re-exercises the `SPEC_odmd_optical` degenerate-reference fix on a second operator.
- **What we can claim if gates pass:** the J table (numbers `arXiv:2501.10320` did not report —
  its focus was charge gaps), the analytic pin, the Heisenberg-failure ladder, the exact
  selection-rule complementarity between the spin and optical channels, and the local-moment
  ladder (the Sᶻ spectral weight ‖Sᶻ|ψ₀⟩‖² falls from 1.000 to 0.759 F → I: charge fluctuations
  eat a quarter of the iodide's moment).
- **What we cannot claim:** isolated-dimer J (no inter-dimer/kagome in-plane exchange — this is
  the *interlayer* coupling of one breathing-trimer pair); density-density interactions only
  (the paper's few-meV non-density-density terms would shift J at that scale — relevant for F,
  where J = 0.051 meV is smaller than those neglected terms: quote F's J as "≈ 0, below the
  model's own accuracy"); exact statevector signals.

## 3. Approach

`spin_excitation_lines(mh, reference)` = `absorption_lines` with the JW-mapped S₁ᶻ − S₂ᶻ kick.
References: the closed form above (from the 2×2 even singlet block vs the triplet at Us);
`dimer_exchange_heisenberg = 4t²/(U₀−Us)` as the perturbative anchor; dense ED implicitly via
`SPEC_odmd_optical`'s validated extraction. ω convention as in the optical spec.

## 4. Public interface

```
odmd_spin.dimer_staggered_moment() -> sparse     # S1z - S2z, JW block ordering, 4 qubits
odmd_spin.spin_excitation_lines(mh, reference=None, n=24) -> (omegas, weights)
odmd_spin.dimer_exchange_analytic(U0, t, Us) -> float     # sqrt(((U0-Us)/2)^2+4t^2)-(U0-Us)/2
odmd_spin.dimer_exchange_heisenberg(U0, t, Us) -> float   # 4 t^2 / (U0 - Us)
```

## 5. Acceptance criteria (validation gates)

`tests/test_odmd_spin_spec.py`; all four LT-bulk materials.

- **G1 — one line, at exactly J.** Each dimer's Sᶻ-kicked spectrum has exactly ONE line, at the
  analytic J to < 1e-9 relative, with weight ‖Sᶻ|ψ₀⟩‖² (< 1e-9) — and the eigenstate
  short-circuit terminates instantly (the `SPEC_odmd_optical` fix, on a second operator).
- **G2 — channel complementarity.** For every material the spin line and the optical line never
  coincide (distance > 100 meV; measured ≥ 528), and the spin line lies BELOW the optical line —
  magnetic excitations are the low-energy physics of every member.
- **G3 — the J table and the Heisenberg failure (DEFINITION OF DONE).** J strictly increasing
  F → Cl → Br → I (0.051 → 245.9 meV); the Heisenberg estimate errs < 1% for F (perturbative
  anchor) and > 30% for I (measured 46.5%) with the error strictly increasing across the family
  — the iodide is beyond the Heisenberg regime.
- **G4 — the local-moment ladder.** ‖Sᶻ|ψ₀⟩‖² = 1 to < 1e-3 for F (pure-spin limit) and strictly
  decreasing F → I (measured 1.000, 0.944, 0.890, 0.759) — charge fluctuations reduce the local
  moment by 24% in the iodide.

## 6. Implementation plan (test-first)

1. `tests/test_odmd_spin_spec.py` encoding G1–G4 (RED — module missing).
2. `odmd_spin.py` — one operator + two closed forms; extraction is `absorption_lines` unchanged.
3. `make gates`.

## 7. Out of scope

- In-plane (kagome) exchange, inter-dimer coupling, magnon dispersion — needs larger clusters.
- S⁺/S⁻ channels (identical physics by SU(2) here), field-dependent spectra, finite temperature.
- Shot noise / circuits on spin signals (compose with the device/UQ specs later).

## 8. Caveats and risks

- **R1:** Nb₃F₈'s J (0.051 meV) sits below the model's own neglected terms — report it as
  "≈ 0", not as a four-significant-figure prediction.
- The dimer J is the interlayer coupling only; the materials' magnetism also involves in-plane
  physics this cluster cannot see.

## 9. Deliverables

- `odmd_spin.py` (+ `__main__` J table); `tests/test_odmd_spin_spec.py` (G1–G4).
- `BACKLOG.md` entry with the J / Heisenberg-error / local-moment ladders.
