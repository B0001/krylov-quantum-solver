# SPEC: Does centering change how shot noise in ⟨Hⁿ⟩ amplifies into PDS(K)?

**Status:** CLOSED (chem-u87), with the G2 rule recorded as non-decisive. **Claim 1 survives:** with
shared-shot (circuit-level) noise, centering does not change noise amplification at all. Raw and
centered PDS agree wherever raw cond(M) < 1e10 (gate 1e-6 Ha; the observed differences are
1e-13 to 1e-10 Ha), and they fail together where the
noisy Hankel matrix is indefinite. **Claim 2 (the pre-registered rule):** not reproducible. The
hit-rate rule fired at 10⁴, at 10⁵, or at 10⁶ depending on the run, so it gives no verdict (G2,
revised). **What held in every run:** centering's float64 gain does survive shot noise *in kind*.
At 10⁶ snapshots, centered PDS(8) has a median error of 0.089 mHa against 2.4–3.3 mHa at K=4, and
it beats raw at K=7–8 by 2–4× in median. Raw PDS(7) fails catastrophically (> 10 mHa) in about
30–39% of trials at every budget ≥ 10⁵, even when shot noise is negligible (rescaled 10¹²). Centered
PDS(7) fails at that rate only from noise, and its failure rate falls to 0 as shots grow. **G3:**
the independent-per-moment noise model that a prior attempt used is about 340× too pessimistic at
K=8. That is why that attempt concluded "centering does not survive".

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

`SPEC_centered_pds` showed that centering cuts cond(M) by 17 orders at K=8 on N₂ CAS(6,6) and
unlocks variational, monotone PDS to K=8, with **exact statevector moments**. On hardware the
moments are estimated from shots. Two claims, both falsifiable:

1. **Centering does not change noise amplification.** When all moments come from one shared shot
   record (the circuit-level case), raw and centered moments are linear functionals of the *same*
   estimated state, and PDS is shift-covariant for *any* moment sequence. Raw and centered PDS
   are then the same function of the same data. Any difference comes from float64 alone. **Dies
   if** raw and centered differ by > 1e-6 Ha where raw cond(M) < 1e10 (G1).
2. **Whether the float64 gain survives shot noise** is decided by a pre-registered hit-rate rule
   (G2), not by eye.

A prior attempt (local only, not reproduced here) concluded that "centering does not survive shot
noise". It modelled independent Gaussian noise per moment. G3 tests whether that model is
pessimistic compared with correlated shots.

## 2. Background and honest framing

- **Noise model: circuit-level, shared shots.** One random-Pauli classical-shadow dataset of |HF⟩
  (the Huang–Kueng–Preskill estimator in `classical_shadows.py`) feeds every moment μ₁..μ₁₅, raw
  and centered. Snapshot ρ̂ₛ = ⊗_q (I + 3 s_q P_q)/2, and μₙ = Tr(Hⁿ ρ̄). Noise across powers is
  therefore correlated exactly as on a device that measures once and post-processes.
- **Symmetry projection.** H commutes with the (N_α, N_β) projector P, so ⟨P Hⁿ⟩ comes from the
  same snapshots (standard post-selection). Only the 400×400 sector block of ρ̄ is needed. That
  block is built as a hi/lo-qubit GEMM, which makes 10⁶ snapshots cheap. Normalizing by ⟨P⟩ is
  irrelevant to PDS, which is invariant to a common scale of the moments.
- **What we can claim.** Measured PDS(K) error versus shadow budget (10³–10⁶ snapshots, sampled
  directly) for raw and centered frames on N₂ CAS(6,6), K=4..8, plus a verdict from G2.
- **What we cannot claim.** (a) Random-Pauli shadows on the HF basis state only. Grouped,
  derandomized or Hadamard-test moment protocols are not modelled. (b) No gate noise or readout
  error: shot noise only. (c) Budgets above 10⁶ (the "scaled" rows) rescale a sampled deviation by
  √(S₀/S). This keeps the covariance exact but not the higher cumulants, so those rows are a lead
  and are not gated. (d) Symmetry projection keeps all noise weight inside the sector, whose
  minimum is E₀. Without projection, noise weight on lower-lying states in other sectors could drag
  PDS below E₀. For N₂ CAS(6,6) the global minimum equals the sector minimum, so this example
  cannot test that.

## 3. Approach

`centered_pds_shots.noise_study`: for each of T=16 trials (seeds 2000..2015), draw 10⁶ snapshots
and use nested prefixes for 10³, 10⁴ and 10⁵. Compute the weights d = diag(V†ρ̄_sec V) in the
sector eigenbasis, then the moments Σ(λᵢ−μ)ⁿ dᵢ with μ = 0 (raw) and with μ = the
`centered_frame` center (centered). Feed both to `moment_expansion.pds_energy`, unchanged.
**Reference:** E₀ = the lowest HF-reachable sector eigenvalue, which equals FCI. G0 checks that it
equals the global minimum.

**Independent contrast (G3), same marginals:** moment n of trial t is taken from trial
(t+n) mod T, separately in each frame. This is the "measure each power on its own shots" model.

## 4. Public interface

```
centered_pds_shots.sector_eigensystem(mh) -> (lam, V, hf_bits, lo, hi)
centered_pds_shots.sample_hf_shadow(hf_bits, n_shots, rng) -> (bases, signs)
centered_pds_shots.shadow_sector_weights(bases, signs, V, lo, hi) -> d   # unnormalized
centered_pds_shots.noise_study(mh, trials=16, max_k=8, seed=2000) -> dict
centered_pds_shots.hit_rate(err, tol=1.6e-3); survives(res, b)          # the G2 rule
uv run python centered_pds_shots.py [trials]                             -> the §6 table
```

## 5. Acceptance criteria (pre-registered before the gated run)

Gate thresholds were fixed after one 4-trial scout (seeds 1000..1003). The gated run uses fresh
seeds (2000..2015). Gates are in `tests/test_centered_pds_shots_spec.py`, with BLAS pinned to 1
thread. Chemical accuracy: tol = 1.6 mHa.

- **G0 — plumbing.** (a) Noiseless sector moments equal `hamiltonian_moments` to 1e-12 relative.
  Noiseless centered PDS(4..8) reproduces `SPEC_centered_pds` §6 (1.076, 0.260, 0.0540, 0.00587,
  0.00106 mHa) to 1%. The sector E₀ equals the global minimum. (b) On H₂, with snapshots from
  `classical_shadows.collect_classical_shadow`, Σλᵢdᵢ equals
  `shadow_energy_samples(·, P H P).mean()` to 1e-10. This is the same estimator, not a
  re-derivation. (c) `sample_hf_shadow` matches `collect_classical_shadow` on N₂ HF: Z-basis
  outcomes are identical and deterministic, and X/Y outcomes have mean within 4/√n of 0.
- **G1 — claim 1: same noise amplification.** At every sampled budget, every trial and every K with
  raw cond(M) < 1e10: |PDS_raw − PDS_cent| < 1e-6 Ha. Non-vacuous: at those K the median shot-noise
  error is ≥ 1e-4 Ha, so the noise is ≥ 100× the frame difference.
- **G2 — claim 2: the verdict rule.** Let hit_f(K) be the fraction of trials with
  |PDS_f(K) − E₀| < 1.6 mHa (a NaN or complex-only root counts as a miss). The float64 gain
  **survives at budget S** iff some K ∈ 5..8 has hit_cen(K) ≥ hit_cen(4) + 0.25 (going past the
  old K=4 cap helps) **and** hit_cen(K) ≥ hit_raw(K) + 0.25 (raw cannot match it there).
  **Hypothesis (from the scout): survives at S = 10⁶ and does not survive at S = 10³.**
- **G3 — the independent-noise model is pessimistic.** At S = 10⁶ the median centered |error| at
  K=8 is ≥ 10× larger with independent-per-moment noise than with shared-shot noise.

**G1 amended during implementation (the claim's tolerance is unchanged).** (i) Where noisy PDS has
no real root (an indefinite Hankel matrix), the gate requires *both* frames to return NaN together.
Previously NaN − NaN poisoned the comparison. (ii) The non-vacuity guard "K=1..3 always
comparable" was too strict: raw cond(M₃) ≈ 1e9 sits near the 1e10 mask, and 1/16 trials at 10⁶
crossed it. The guard is now K=1,2 always and K=3 in ≥ 90% of trials. (iii) In the noise-size
median, NaN counts as an infinite error (a miss).

**G2 revised (post hoc, after the rule proved non-reproducible; the original rule is kept above as
the record).** At S = 10⁶: centered median |err|(K=8) ≤ centered median |err|(K=4) / 10, and the
centered median is below the raw median at K=7 and at K=8. The rule's own verdict is reported in
§6 and is not gated.

## 6. Measured data

N₂ CAS(6,6), E₀ = −107.62310177 Ha, μ = −9.0032 Ha, BLAS pinned to 1 thread.

**Gated run (16 trials, seeds 2000..2015)**, median |PDS − E₀| in mHa, raw / centered, and
chemical-accuracy hit rate:

| S | K=4 | K=5 | K=6 | K=7 | K=8 | hit K=8 raw / cen |
|---|---|---|---|---|---|---|
| 10³ | 221 / 221 | 24.7 / 24.7 | 22.0 / 22.0 | 5.49 / 4.94 | 4.73 / 2.62 | 0.25 / 0.38 |
| 10⁴ | 13.7 / 13.7 | 10.5 / 10.5 | 4.50 / 4.51 | 8.21 / 9.07 | 204 / 1.04 | 0.06 / 0.69 |
| 10⁵ | 10.6 / 10.6 | 7.67 / 7.67 | 2.62 / 2.63 | 3.90 / 0.88 | 1.01 / 0.36 | 0.62 / 0.69 |
| 10⁶ | 2.39 / 2.39 | 1.63 / 1.63 | 0.51 / 0.46 | 0.44 / 0.23 | 0.29 / 0.089 | 0.75 / 0.88 |
| 10¹⁰ (scaled) | 1.08 / 1.08 | 0.26 / 0.26 | 0.051 / 0.053 | 0.043 / 0.0058 | 0.023 / 0.0023 | 0.75 / 1.00 |
| 10¹² (scaled) | 1.08 / 1.08 | 0.26 / 0.26 | 0.056 / 0.054 | 0.042 / 0.0059 | 0.028 / 0.0011 | 0.75 / 1.00 |

Independent-per-moment contrast (G3), centered K=8 median at 10⁶: **30.2 mHa** against 0.089 shared
(340×). At 10¹² it is still 0.31 mHa against 0.0011.

**64-trial rerun (seeds 3000..3063), not gated**:

| S | catastrophic (> 10 mHa) K=7 raw / cen | K=8 raw / cen | median K=8 raw / cen (mHa) | hit@0.1 mHa K=8 raw / cen |
|---|---|---|---|---|
| 10⁵ | 0.39 / 0.17 | 0.22 / 0.17 | 0.62 / 0.23 | 0.05 / 0.31 |
| 10⁶ | 0.39 / 0.17 | 0.20 / 0.22 | 0.35 / 0.088 | 0.19 / 0.55 |
| 10⁸ (scaled) | 0.31 / 0.08 | 0.17 / 0.22 | 0.065 / 0.023 | 0.72 / 0.72 |
| 10¹⁰ (scaled) | 0.38 / 0.00 | 0.06 / 0.00 | 0.023 / 0.0022 | 0.88 / 1.00 |

**The G2 rule's verdict by run** (S = 10³ / 10⁴ / 10⁵ / 10⁶): 16-trial script run F/T/F/F; the first
16-trial pytest run (same seeds) T at 10⁶; the 64-trial run F/F/T/F. Same seeds giving different
raw-frame hits is the BLAS-path dependence that `SPEC_centered_pds` R3 records for raw K ≥ 6. The
margin is about 1.5σ at 16 trials. **Why the rule was the wrong instrument:** noiseless PDS(4) is
already at 1.076 mHa < 1.6 mHa, so at low noise a chemical-accuracy hit rate cannot see a gain at
K>4. A 0.1 mHa threshold would have (see the 64-trial table), but it was not pre-registered.

**Other findings.**
- **Negative variance.** At 10³ snapshots, 2/16 trials estimate ⟨H²⟩ − ⟨H⟩² < 0. The Hankel matrix
  is then indefinite, P₂ has complex roots, and PDS returns no root, in both frames identically
  (G1 checks the NaN masks match).
- **Variationality under noise.** The noisy weights form a *signed* measure. At 10⁶, 12–27% of
  centered trials sit more than 1.6 mHa below E₀ at K=4..8 (64-trial run), and this falls to 0 by
  10¹⁰. PDS is not variational on shot-noisy moments.
- **Why shared shots are benign.** Shared-shot moments are the exact moments of a (signed) measure
  on the *exact* spectrum of H. Only the weights are noisy. PDS(K) converges to the bottom of that
  support, which is E₀. Independent per-moment noise destroys that structure.

## 7. Out of scope

- Grouped/derandomized shadows, Hadamard-test moments, gate noise and readout error.
- Systems without symmetry projection, where noise can see lower-lying states in other sectors
  (§2 d).
- Changing `moment_expansion` defaults.
