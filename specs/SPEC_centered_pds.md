# SPEC: PDS's high-K accuracy floor is mostly an artifact of the uncentered moment frame

**Status:** CLOSED — gates G1–G5 (+G2b) PASS (2026-10-01, chem-70c). **The claim survives, with one
amendment.** At K=8, rebuilding the same moments in `device_odmd.centered_frame` cuts cond(M) by
**13.0 (H₄), 11.3 (LiH), 17.4 (N₂ CAS(6,6)) orders**, well above the 8-order kill line. Centered PDS
matches raw PDS to ≤ 8e-12 Ha wherever raw cond(M) < 1e10. Centered PDS stays variational and
monotone through K=8 (N₂ 1.06 µHa, H₄ 0.06 µHa, LiH 0.091 mHa above FCI). **Amendment
(the boundary):** "det M → 0" is not *only* an artifact. When K exceeds the number of HF-reachable
eigenstates (H₂: 2), M is exactly singular in every frame and centering buys only 2.3 orders. That
part of `SPEC_moment_pds` R1 is fundamental, though harmless: PDS is already exact at
K = n_reachable.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

`SPEC_moment_pds` R1 treats "det M → 0 at large K" as fundamental and caps PDS gates at K ≤ 4. The
claim here is that this cap mostly comes from the frame. PDS is shift-covariant in exact arithmetic,
but raw electronic-frame moments grow like ‖H‖ⁿ, so the Hankel matrix M carries a dynamic range of
~‖H‖^(2K−2). Centering H on the HF-reachable band should remove most of it at identical PDS
values. **Dies if** centered ≠ raw by > 1e-6 Ha anywhere raw cond(M) < 1e10 (not the same
functional), **or** centering gives < 8 orders of cond(M) reduction at K=8.

## 2. Background and honest framing

- **Numerical hygiene, not physics.** Centering does not change PDS's convergence *rate*; it only
  pushes back the point where float64 gives out. Same functional as Peng & Kowalski (Quantum 5, 473,
  2021), with the reference and FCI taken from `SPEC_moment_pds`.
- **What we can claim.** On H₄, LiH and N₂ CAS(6,6), centered PDS(K) is a variational upper bound
  that tightens monotonically to K=8 and is deterministic across BLAS thread counts. Raw PDS breaks
  down at K≥6 on N₂ in a thread-dependent way. The centering μ need not be accurate: μ = ⟨H⟩ (the
  first moment, free) works.
- **What we cannot claim.** (a) No speedup in convergence. (b) Exact statevector moments; on
  hardware, ⟨(H−μ)ⁿ⟩ has the same measurement cost as ⟨Hⁿ⟩ but different shot-noise propagation
  (out of scope). (c) Krylov exhaustion (K > n_reachable) remains a hard limit (G2b). (d) LiH at K=8
  still has centered cond(M) ≈ 9e12, so the next wall is only a few K away. Centering buys K≈8, not
  arbitrary K.

## 3. Approach

Composition only: `hamiltonian_moments(centered_frame(mh)[0], 2K)` → `pds_energy`. For μ
sensitivity, the test applies the same `dataclasses.replace` shift as `centered_frame` at
μ_center + δ and at μ = ⟨H⟩. No library code changed. Reference: FCI (dense diagonalization).

## 4. Public interface

None new. Recommended usage for K > 4:
`pds_energy(*hamiltonian_moments(centered_frame(mh)[0], 2*K), K)`. At scale, shift by ⟨H⟩ (see G5).

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_centered_pds_spec.py`, with BLAS pinned to 1 thread.

- **G1 — same functional.** |PDS_cent − PDS_raw| < 1e-6 Ha at every K ≤ 8 where raw cond(M) < 1e10,
  on H₂/H₄/LiH/N₂ (≥ 2 K values compared per system). Measured max diff 8.2e-12 Ha.
- **G2 — the claim.** log₁₀ cond_raw − log₁₀ cond_cent ≥ 8 at K=8 on H₄/LiH/N₂.
- **G2b — boundary.** H₂ (2 reachable states): PDS(2) is exact, centered cond(M) > 1e14 at K=3, and
  the K=8 reduction is < 8 orders. Here det M → 0 is real.
- **G3 — the unlock (K=7–8 gates).** Centered PDS(K) ≥ FCI − 1e-9 and PDS(K+1) ≤ PDS(K) + 1e-9 for
  K=1..8; PDS(8) − FCI < 1e-6 (H₄), < 1.6e-3 (LiH), < 5e-6 Ha (N₂); and PDS(8) error < ½ PDS(4) error.
- **G4 — raw breaks (reproduction).** N₂: raw cond(M) > 1e17 at K=6,7,8 and raw departs from
  centered by > 1e-5 Ha somewhere in K=6..8.
- **G5 — μ sensitivity.** N₂: δ ∈ {±0.5, ±1, ±2} Ha leaves PDS(7), PDS(8) unchanged to < 1e-6 Ha.
  μ = ⟨H⟩ reproduces centered PDS at all K ≤ 8 (< 1e-6 Ha) and keeps a ≥ 8-order cut at K=8.

## 6. Measured data (centered_frame μ; err = PDS − FCI)

cond(M) at K=8: raw → centered → ⟨H⟩-shifted

| system | n_reach | raw | centered | ⟨H⟩-shift | orders cut (cent / ⟨H⟩) |
|---|---|---|---|---|---|
| H₂ | 2 | 3.8e19 | 2.1e17 | 2.2e18 | 2.3 / 1.2 |
| H₄ | 12 | 3.8e20 | 4.1e7 | 1.2e11 | 13.0 / 9.5 |
| LiH | 27 | 1.6e24 | 8.6e12 | 1.1e12 | 11.3 / 12.2 |
| N₂ CAS(6,6) | 18 | 2.3e26 | 8.6e8 | 1.8e12 | 17.4 / 14.1 |

N₂ error (mHa) by K, centered (identical at 1/2/4 BLAS threads to 5 decimals):
K=4 1.076, 5 0.260, 6 0.0540, 7 0.00587, 8 0.00106.

N₂ raw, K=6/7/8 (mHa), six fresh processes:

| BLAS threads | K=6 | K=7 | K=8 | failure |
|---|---|---|---|---|
| 1 | 0.0340 | 0.0464 | 0.0287 | non-monotone 6→7 |
| 2 | 0.0416 | 0.0477 | 0.0262 | non-monotone 6→7 |
| 4 | 0.0564 | 0.0378 | −337 045 | non-variational K=8 |
| default (run A) | 0.0572 | −18 754 | 0.0315 | non-variational K=7 |
| default (run B) | 0.0686 | 0.0371 | −23 893 | non-variational K=8 |
| default (run C) | 0.0581 | 0.0417 | −977.7 | non-variational K=8 |

The earlier scout's 0.044 → 0.049 (K=6→7) is the 1–2-thread path. Every raw run broke `SPEC_moment_pds`
G2 or G3 somewhere in K=6..8. Every centered run passed both.

μ error on N₂ (half-width W/2 = 2.39 Ha): δ ∈ [−5, +2] Ha changes PDS(6..8) by ≤ 0.07 µHa;
δ = +5 makes PDS(8) non-variational (−1.8 µHa); δ = +10 Ha breaks it (PDS(7) −311 mHa). The
response is asymmetric: shifting *up* past the band pushes the ground state away from 0 and
re-inflates the dynamic range. μ = ⟨H⟩ sits at δ ≈ −2.26 Ha, inside the safe window.

## 7. Out of scope

- Shot noise on ⟨(H−μ)ⁿ⟩ and whether centering changes noise amplification.
- Higher-root/excited-state use (`SPEC_pds_excited_roots`). *Lead, not gated:* on stretched H₄
  the centered frame kept every root above its reachable target at K=3..8 at 1 and 4 threads
  (min margin +0.0042 mHa at K=8), but this container's raw run did not reproduce that spec's
  −7.49 mHa violation, so "centering fixes it" is not yet shown.
- Changing `moment_expansion` to center by default (would alter a closed spec's API).

## 8. Caveats and risks

- **R1 — validation-scale μ.** `centered_frame` diagonalizes densely. G5 shows that μ = ⟨H⟩ is good
  enough, so this is not a blocker.
- **R2 — LiH's next wall.** Centered cond(M) at K=8 is already 8.6e12 for LiH. Expect float64 trouble
  again around K ≈ 10–11 for larger or denser spectra.
- **R3 — raw-frame G4 is path-dependent.** Only the path-independent signature is gated (cond > 1e17,
  deviation > 1e-5 Ha). The specific G2/G3 violation is recorded in §6, not asserted.
