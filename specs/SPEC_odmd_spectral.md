# SPEC: ODMD spectroscopy — Green's-function poles and weights from survival amplitudes

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_odmd_spectral_spec.py`).

---

## 1. Goal

Every validated number in this repo is an energy; experiments measure **spectra**. Claim: the
one-particle Green's function's Lehmann representation is exactly ODMD-shaped — the survival
amplitude of a particle-removed reference `a_i|ref⟩` is
`s_k = Σ_n |⟨E_n^{N−1}|a_i|ref⟩|² e^{−iE_n kτ}`, so the DMD **poles are the (N−1)-sector
eigenvalues (ionization lines)** and the **Vandermonde amplitudes are the spectral weights**;
likewise `a_i†` gives the electron-addition side. Together: the photoemission /
inverse-photoemission spectrum A(ω) from the same 1-D signals the ODMD stack already measures.
Falsifiable line by line against the exact Lehmann representation (dense diagonalization), and
cross-pinned to `SPEC_nb3x8_device_gap.md`: the distance between the Nb₃I₈ dimer's Hubbard bands
must equal the recorded 842.44 meV charge gap.

## 2. Background and honest framing

- Green's functions from real-time propagation of `a_i|ψ⟩` on quantum computers are established
  (e.g. Kosugi & Matsushita, PRA 101, 012330 (2020); Endo et al.); ODMD supplies the pole/weight
  extraction. Composition of pinned primitives; the falsifiable line-by-line packaging and the
  damping-immunity finding for *intensities* are what is new here.
- **What we can claim if gates pass:** machine-precision poles and degeneracy-aggregated weights
  against exact Lehmann; a material cluster's full A(ω) whose band gap reproduces the capstone
  number; and that under uniform damping **both** energies and intensities survive (damping goes
  into |λ|, amplitudes are untouched — extending `SPEC_device_odmd`'s phase-immunity to the
  whole spectrum).
- **What we cannot claim:** exact statevector signals (no circuits/shot noise here — compose
  with the device/UQ specs later); the default HF reference gives *HF-referenced* weights, which
  differ from the true Lehmann weights (exact ground-state reference) by ~2% on H₄ — measured
  and gated, not ignored; DMD merges exactly degenerate lines into one (their weights add —
  physically correct for A(ω), but individual degenerate components are not resolved); weak
  lines below the amplitude floor are invisible (the visibility law); isolated-cluster caveat
  for the Nb₃X₈ numbers as in the capstone.

## 3. Approach

`ladder_operator` (qiskit-nature `FermionicOp` → the vetted Jordan–Wigner mapper) applies
`a_i`/`a_i†` to the reference statevector; each (normalized) reference gets its own centered
frame (μ, τ = π/W over ITS reachable spectrum — evolution conserves N, so the sector is pinned);
`odmd_spectrum` extracts eigenphases + Vandermonde amplitudes; weights are rescaled by
`‖a_i|ref⟩‖²`; lines aggregate over spin-orbitals. ω conventions: removal `ω⁻ = E₀^N − pole`,
addition `ω⁺ = pole − E₀^N`; photoemission gap = `min ω⁺ − max ω⁻`. References: dense
diagonalization (eigenvalues, |⟨E_n|ψ_ref⟩|² Lehmann weights, sector grounds).

## 4. Public interface

```
odmd_spectral.SpectralLine                    # dataclass: omega, weight, orbital, kind
odmd_spectral.ladder_operator(kind, i, n_so) -> sparse matrix    # '-'/'+' JW ladder op
odmd_spectral.reference_signal(mh, psi_raw, n=24) -> (s, tau, mu, nrm2)
odmd_spectral.lines_from_signal(s, tau, mu, nrm2, amp_floor=1e-4, mod_window=0.2)
    -> (poles_electronic, weights)
odmd_spectral.greens_function_lines(mh, kind, orbitals=None, reference=None, n=24)
    -> list[SpectralLine]                     # reference None -> HF
odmd_spectral.photoemission_gap(mh, reference=None, n=24) -> float
```

## 5. Acceptance criteria (validation gates)

`tests/test_odmd_spectral_spec.py`. Exact references computed independently in the test.

- **G1 — Lehmann, line by line (DEFINITION OF DONE).** H₄ removal, every occupied spin-orbital:
  each visible exact line (degeneracy-aggregated weight > 1e-3) is matched by an ODMD pole
  within 1e-8 Ha with aggregated weight error < 1e-4 (measured: strong lines ≤ 1e-13 Ha; the
  weakest satellite, weight ≈ 0.002, needed K=32 rather than 24 to clear the gate — the
  visibility law at work; gates run at K=32).
- **G2 — reference honesty.** With the exact ground state as reference, weights match the TRUE
  Lehmann weights < 1e-4; with the HF reference the deviation on H₄ is real but bounded
  (max line error in (0.005, 0.05); measured ≈ 0.023) — the recorded approximation.
- **G3 — the material's spectrum, cross-pinned.** Nb₃I₈ dimer (exact N=2 reference): ≥ 2
  distinct lines on each side (both Hubbard bands), and
  `photoemission_gap = min ω⁺ − max ω⁻` reproduces the sector-FCI charge gap
  (`nb3x8_device_gap.exact_gap`) within 0.01 meV (measured exact: 1073.22 − 230.78 = 842.44).
- **G4 — intensities are damping-immune too.** H₄ removal reference, `s → 0.7^k s` (30% loss
  per step), wide modulus window: every pole moves < 1e-6 Ha and every weight < 1e-6 (measured
  ~1e-12) — uniform damping enters |λ| only; `SPEC_device_odmd`'s immunity extends from
  energies to the whole spectrum.

## 6. Implementation plan (test-first)

1. `tests/test_odmd_spectral_spec.py` encoding G1–G4 (RED — module missing).
2. `odmd_spectral.py` — ladder ops + per-reference frames; extraction is `odmd_spectrum`
   unchanged.
3. `make gates`.

## 7. Out of scope

- Shot noise / circuits / device noise on spectral signals (compose with `SPEC_device_odmd` and
  `SPEC_odmd_uq` later); broadening/self-energy extraction; two-particle (optical) response.
- Correlated reference preparation on circuits (exact ψ₀ used at validation scale only).
- Resolving individual degenerate components (physically unobservable in A(ω) anyway).

## 8. Caveats and risks

- **R1 — weight floor:** lines with aggregated weight below `amp_floor` (default 1e-3) are
  invisible; satellites in strongly correlated references may need deeper K (visibility law).
- HF-referenced weights are NOT true Lehmann weights (G2 bounds the difference at mild
  correlation; expect worse for strong correlation).
- The Nb₃X₈ spectrum is the isolated dimer's (no band broadening) — see the capstone caveat.

## 9. Deliverables

- `odmd_spectral.py` (+ `__main__`: the Nb₃I₈ dimer A(ω) table).
- `tests/test_odmd_spectral_spec.py` — gates G1–G4.
- `BACKLOG.md` entry.
