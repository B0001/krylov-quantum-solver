# SPEC: The HF overlap certificate's noise sibling — coverage bifurcates, and inflation can hurt

**Status:** IMPLEMENTED, gates green (`tests/test_hf_overlap_noise_spec.py`, 10/10 pass).
Closes the backlog hypothesis *"The overlap certificate is the one rung with no noise sibling"*
(`specs/BACKLOG.md`, Certified-bounds arc) / bead `chem-mqu`.

---

## 1. Goal

`hf_overlap_certificate.certify_hf_overlap` (`specs/SPEC_hf_overlap_certificate.md`) is the one
certificate in the certified-bounds arc with no shot-noise study. Unlike the Ritz-state
certificates already studied (`certified_noise`, `gap_selfcheck_noise`,
`certified_thermochem_noise`), its guiding state u = HF is a **determinant** — a known
computational-basis state — so λ_u = ⟨HF\|H\|HF⟩ and r = ‖(H − λ_u)HF‖ are **exact, classical**
computations with **zero shot noise**. Every bit of noise in this certificate funnels through the
single gap input β (the self-mode Weinstein floor θ₁ − σ₁), and γ_min = √(1 − r²/δ²) with
δ = β − λ_u has an **unbounded derivative** as δ → r.

**Claim under test:** because of this structural difference, coverage does not degrade smoothly
across the margin (δ_exact − r) the way it does elsewhere in this repo's noise arc — it
**bifurcates**: near-perfect at wide margin, collapsing sharply as margin thins. **Kill
criterion (pre-registered in the bead):** dies if z ≤ 2 inflation restores coverage ≥ 0.9 on every
system tested — that would just be a sixth repetition of "inflation buys coverage" and not a new
finding.

## 2. Background and honest framing

- Builds on `certified_noise.md` (i.i.d. Gaussian shot-noise model, λ-1-norm standard errors,
  `certified_half_width = z·λ_H/√N`) and `gap_selfcheck_noise.md`'s post-hoc, single-sided pad
  convention (`_pad`: shift a lower bound further down by the half-width, never up).
- **What you can claim if the gates pass:** (1) the determinant-only shortcut is real —
  reproduced here as "only β is noisy"; (2) coverage sensitivity to noise is sharply
  non-uniform across the margin spectrum, spiking well before the certificate goes structurally
  vacuous; (3) **z-inflation is not uniformly good** for this certificate — it fixes the
  wide-margin failure mode (overclaim) but *worsens* the moderate/thin-margin failure mode
  (vacuousness), a genuinely new, opposite-sign result relative to every other noise spec in this
  repo, all of which found z uniformly buys coverage (at various costs).
- **What you cannot claim:** a closed-form threshold for where the sensitivity spike sits (measured
  on four systems, not derived); that inflation is "bad" in general — it only hurts the specific
  failure mode (vacuousness) that dominates at moderate/thin margin, and it remains the only lever
  against overclaim at wide margin; anything about the chained-Ritz overlap bound
  (`krylov_refine.py`, `SPEC_chained_overlap.md`) — **that state is a Krylov Ritz vector, not a
  determinant, so its own λ_v, r_v need ⟨v\|H\|v⟩ and ⟨v\|H²\|v⟩, which DO cost shots.** The
  shot-free-residual argument is specific to determinant guiding states; nothing here generalizes
  to the chained bound, and that boundary must not be discovered later as a surprise.
- **Reference:** `hf_overlap_certificate.exact_reachable_overlap` (dense `eigh`, the same reference
  `SPEC_hf_overlap_certificate`'s own G1 uses) — the ground truth every "valid"/"invalid" label
  below is checked against.

## 3. Approach

Reuse, don't rederive:
- `certified_gaps.gap_bracket`'s exact self-mode formula (β = θ₁ − σ₁) rebuilt from a **noisy**
  realization of (θ₁, ⟨H²⟩₁) at N shots — `certified_noise.hamiltonian_one_norms` supplies λ_H,
  λ_H² for the standard errors, exactly `certified_noise`'s model.
- `certified_noise.certified_half_width(λ_H, N, z)` for the inflation magnitude, applied **one-
  sided** (subtracted from β only — never symmetric), mirroring `gap_selfcheck_noise._pad`'s
  "-half_width on the lower bound" half applied to a single scalar instead of a two-sided bracket.
- `hybrid_quantum_solver.certified_overlap.{rayleigh_quotient,residual_norm}` for λ_u, r — computed
  **once, exactly**, never touched by the noise model (the claim this module tests).

For each of 4 systems × 3 shot counts × 4 inflation levels (≥ 6000 seeded Monte-Carlo trials per
cell), classify every trial into exactly one of:
- **covered** — non-vacuous AND γ_min(noisy) ≤ exact overlap: the useful, correct case.
- **vacuous** — r ≥ δ(noisy): conservative "can't certify," never wrong.
- **invalid** — non-vacuous but γ_min(noisy) > exact overlap: a silent overclaim, the failure mode
  a certificate must never produce.

Systems (self-mode, M = 8, matching `hf_overlap_certificate`'s default and clearing the M ≥ 6
premise), spanning the margin exactly as the bead specifies, wide to vacuous:

| system | margin = δ_exact − r | r/δ | exact overlap |
|---|---|---|---|
| H₂ eq (0.74 Å) | **+1.419** (wide) | 0.113 | 0.9936 |
| H₄ chain (0.9/1.8/2.7 Å) | +0.151 (moderate) | 0.649 | 0.9769 |
| H₄ chain (1.0/2.0/3.0 Å) | +0.053 (thin) | 0.842 | 0.9677 |
| square H₄ a = 1.2 Å | **−0.255** (vacuous already at zero noise) | n/a | 0.0000 |

## 4. Public interface

```
hf_overlap_noise.exact_margin(mh, m, solver=None) -> dict
    # lam, r, delta, margin, r_over_delta, exact_overlap -- zero-noise reference point (x-axis)

hf_overlap_noise.hf_overlap_noise_coverage(
    mh, m, shots, z=2.0, trials=6000, seed=0, solver=None,
) -> dict
    # coverage, frac_vacuous, frac_invalid, lam_h, lam_h2, half_width, + exact_margin's fields
```

## 5. Acceptance criteria (validation gates)

All measured at M=8, trials=6000, seed=0, shots ∈ {1e4,1e5,1e6}, z ∈ {0,1,2,3} — the exact grid
the bead specifies (`tests/test_hf_overlap_noise_spec.py`, 10/10 green).

- **G0 — the span is real.** margin(H₂) > margin(H₄ moderate) > margin(H₄ thin) > 0 >
  margin(square H₄) at zero noise. *Measured: 1.419 > 0.151 > 0.053 > 0 > −0.255.*
- **G1 — structural sanity (killable).** The one-sided pad is really one-sided: for fixed
  system/shots, `frac_vacuous` is non-decreasing in z and `frac_invalid` is non-increasing in z,
  at every cell. *Zero violations.*
- **G2 — wide margin is near-saturated by z=2.** H₂'s only failure mode is invalid (overclaim);
  z=2 fixes it: coverage ≥ 0.995 at every shot count. *Measured: 0.9988 / 1.0000 / 1.0000
  (1e4/1e5/1e6); frac_invalid ≤ 0.0012.*
- **G3 — THE KILL CRITERION, evaluated explicitly (definition of done).** z ≤ 2 does **not**
  restore coverage ≥ 0.9 everywhere: H₄ moderate and H₄ thin stay below 0.9 at *every* z ∈ {0,1,2}
  and every shot count. *Measured max over that sub-grid: 0.8308 (H₄ moderate, 1e6 shots, z=0) —
  not a near miss.* **The claim survives; this is a genuine sixth-ish finding, not a repeat.**
- **G4 — THE BIFURCATION.** At shots=1e4, the coverage swing from z=0→3 is 3-4× larger at thin
  margin than at wide or moderate margin — sensitivity to the same noise knob spikes sharply
  rather than growing smoothly across the margin axis. *Measured swings: H₂ +0.111, H₄ moderate
  −0.138, H₄ thin **−0.451** (3.9× the H₂ magnitude, 3.3× the H₄-moderate magnitude); square H₄
  stays pinned at 0 (already vacuous, no sensitivity left to have).*
- **G5 — THE SIGN REVERSAL (the novel part).** At shots=1e4, inflation *increases* coverage at
  wide margin (+0.111) but *decreases* it at moderate (−0.138) and thin (−0.451) margin — the
  opposite of `certified_noise`/`gap_selfcheck_noise`, where z uniformly buys coverage. Mechanism:
  the dominant failure mode flips from invalid (fixed by pushing β down) to vacuous (worsened by
  pushing β down further).
- **G6 — an already-vacuous system cannot be rescued by inflation.** square H₄ (vacuous at zero
  noise) has frac_vacuous > 0.999 and coverage < 0.001 at *every* shots/z cell tested — z can only
  push β further down, so it can never manufacture a certificate the noiseless one doesn't have.
  *Measured: frac_vacuous = 1.0000 at all 12 cells (0/6000 trials crossed).*
- **G7 — shots always help, unlike z.** For every system with a positive exact margin, coverage at
  z=0 is non-decreasing in shots (more shots shrink β's noise toward the true, non-negative
  margin). *Measured H₂: 0.889→0.987→1.000; H₄ moderate: 0.537→0.617→0.831; H₄ thin:
  0.502→0.527→0.587.* This is the contrast worth keeping: shots are unconditionally good here,
  z is conditionally good — the two knobs are not interchangeable.

> Definition of done: **G3** (the pre-registered kill criterion) plus **G5** (why it isn't a
> repeat — the sign reversal is new content, not just "more inflation needed" like
> `gap_selfcheck_noise`).

## 6. Implementation plan (test-first)

1. `tests/test_hf_overlap_noise_spec.py` encoding G0–G7 (was RED — module didn't exist).
2. `hf_overlap_noise.py`: `exact_margin` (zero-noise reference point) and
   `hf_overlap_noise_coverage` (Monte Carlo over noisy β only, one-sided z-pad, three-way
   covered/vacuous/invalid classification against the dense exact reference).
3. `uv run python -m pytest tests/test_hf_overlap_noise_spec.py -v` → green; `ruff check` clean.

## 7. Out of scope

- The chained-Ritz overlap bound (`krylov_refine.py`) — **explicitly not covered**; its guiding
  state v is a Krylov Ritz vector, so λ_v and r_v are themselves shot-noisy. Extending this study
  to the chained bound is a distinct, harder problem (two noisy inputs, not one) and a natural
  follow-up, not attempted here.
- A repair for the vacuous-side degradation (this spec measures the bifurcation and the sign
  reversal; it does not propose or test a fix).
- Oracle-mode β (only self-mode is studied — the production path the certificate is meant for).
- More than four systems / a swept continuum of margins (the four are chosen to span wide →
  vacuous per the bead, not to trace a continuous curve).
- Non-Gaussian or grouped-Pauli measurement noise (same idealization as the rest of the arc).

## 8. Caveats and risks

- **R1 — determinant-only scope (state this loudly, per the bead).** The zero-shot-noise
  shortcut for λ_u, r is specific to a guiding state that is a single computational-basis
  determinant. It does **not** carry over to the chained-Ritz bound
  (`SPEC_chained_overlap.md`/`krylov_refine.py`), whose Ritz vector v puts λ_v, r_v back on the
  shot budget. Any future reader tempted to reuse this module's "only β is noisy" framing for a
  chained or Ritz-based guiding state would be wrong; this caveat exists so that mistake doesn't
  get made silently.
- **R2 — the specific margin/swing numbers are measured on four systems, not derived.** The
  *direction* (bifurcation + sign reversal) is the falsifiable claim; the exact swing ratios
  (3–4×) and the z\* location are system-specific measurements, like every other z\* in this arc's
  specs (`gap_selfcheck_noise` R1 makes the same point about its own z\*).
  R3 — one-sided pad convention only. A symmetric two-sided pad on β (rather than
  "conservative-only, subtract half-width") was not tested; the one-sided choice follows directly
  from β being a lower bound (same logic as `gap_selfcheck_noise._pad` on `gap_lower`), but a
  different convention could in principle change the G5 sign-reversal's magnitude (not its
  direction, which follows from which failure mode dominates).

## 9. Deliverables

- `hf_overlap_noise.py` — `exact_margin`, `hf_overlap_noise_coverage`.
- `tests/test_hf_overlap_noise_spec.py` — gates G0–G7, all green.
- This spec, with measured numbers in §5 (not `data/`, which is gitignored).
