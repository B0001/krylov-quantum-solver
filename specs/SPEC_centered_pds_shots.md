# SPEC: Does centering change how shot noise in ⟨Hⁿ⟩ amplifies into PDS(K)?

**Status:** CLOSED — 2026-10-02 (chem-u87). 7/7 gates pass as recorded. **G1 and G3 SURVIVE.**
With one shared shot record, centering does not change noise amplification: raw and centered
agree to 1.5e-8 Ha. Independent per-moment noise is 341× too pessimistic.

**G2 is KILLED.** It pre-registered "the float64 gain survives at 10⁶ but not at 10³". On the
pinned path the rule says it survives at both, each time at exactly its margin.

**G2b is KILLED.** It was registered at review: "that verdict does not depend on the float64
path". The verdict flips with the raw frame's last bits at 10³, 10⁵ and 10⁶.

**Verdict:** centering's gain survives shot noise as path-independent results. Centered PDS is
ahead at K=7. That lead is real but roughly the size of the 25-point margin, so a 16-trial hit-rate
rule cannot certify it (§6). Supersedes `SPEC_centered_pds_shot_noise.md`.

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

A prior attempt, `SPEC_centered_pds_shot_noise.md` (now SUPERSEDED by this spec), concluded that
"centering does not survive shot noise". It added independent noise to each moment after the fact.
Its tests and the `hamiltonian_moments(..., centered=...)` argument they call never reached main.
G3 tests whether that kind of model is pessimistic compared with correlated shots.

## 2. Background and honest framing

- **Noise model: circuit-level, shared shots.** One random-Pauli classical-shadow dataset of |HF⟩
  (the Huang–Kueng–Preskill estimator in `classical_shadows.py`) feeds every moment μ₁..μ₁₅, raw
  and centered. Snapshot ρ̂ₛ = ⊗_q (I + 3 s_q P_q)/2, and μₙ = Tr(Hⁿ ρ̄). Noise across powers is
  therefore correlated exactly as on a device that measures once and post-processes.
- **Symmetry projection.** H commutes with the (N_α, N_β) projector P, so ⟨P Hⁿ⟩ comes from the
  same snapshots. This is a choice of estimator, not shot post-selection: X/Y-basis shots do not
  reveal N. Only the 400×400 sector block of ρ̄ is needed. That block is built from per-half
  snapshot patterns (§5, review additions), which makes 10⁶ snapshots cheap. Normalizing by ⟨P⟩
  does not affect PDS, which is invariant to a common scale of the moments.
- **What we can claim.** Measured PDS(K) error versus shadow budget (10³–10⁶ snapshots, sampled
  directly) for raw and centered frames on N₂ CAS(6,6), K=4..8, plus a verdict from G2.
- **What we cannot claim.** (a) Random-Pauli shadows on the HF basis state only. Grouped,
  derandomized or Hadamard-test moment protocols are not modelled. (b) No gate noise or readout
  error: shot noise only. (c) Budgets above 10⁶ (the "scaled" rows) rescale a sampled deviation by
  √(S₀/S). This keeps the covariance exact but not the higher cumulants, so those rows are a lead
  and are not gated. (d) Symmetry projection keeps all noise weight inside the sector, whose
  minimum is E₀. Without projection, noise weight on lower-lying states in other sectors could drag
  PDS below E₀. For N₂ CAS(6,6) the global minimum equals the sector minimum, so this example
  cannot test that. (e) The classical post-processing is exact: it uses the dense sector
  eigenbasis, at exponential cost. At scale, Tr(Hⁿ ρ̄) needs the Pauli expansion of Hⁿ. That is the
  same estimator at a different cost. (f) The raw frame at K ≥ 7 is one float64 path per platform.
  The pinned-path verdicts in §6 were measured on macOS 27 with Accelerate and may differ on
  another BLAS build, so the gates assert only the forms of the G2 and G2b kills that hold on any
  build.

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
centered_pds_shots.sector_eigensystem(mh, order_seed=None) -> (lam, V, hf_bits, lo, hi)
centered_pds_shots.sample_hf_shadow(hf_bits, n_shots, rng) -> (bases, signs)
centered_pds_shots.shadow_sector_rho(bases, signs, lo, hi) -> rho        # sector block of mean snapshot
centered_pds_shots.shadow_sector_weights(bases, signs, V, lo, hi) -> d   # unnormalized
centered_pds_shots.noise_study(mh, trials=16, max_k=8, seed=2000, paths=(None,)) -> [dict per path]
centered_pds_shots.pinned_n2_study(trials=16, seed=2000) -> [dict per PATHS]  # PYTHONHASHSEED=0, 1 thread
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

**Review additions (chem-u87 batch review, 2026-10-02), registered before the recorded run.** The
draft's G2 docstring already quoted a result ("survives at 10⁶ died, raw 12/16 vs centered 13–14/16
at K=7–8", plus a "64-trial check") that no file recorded. That run's PYTHONHASHSEED was unpinned,
so it could not be reproduced bit for bit. The recorded run below is fresh. A one-off probe
showed why this matters. It ran `sector_eigensystem` → `exact_weights` → raw `pds_energy` in
fresh interpreters. The noiseless raw PDS(7) came out −24 289 mHa at PYTHONHASHSEED=0 and 1, and
+0.0138 mHa at 2. The sorted Pauli labels of H are the same across hash seeds, but the sorted
coefficients are not. The hash seed changes their last bits upstream, and raw cond(M) at K=7 is
past 1/ε, so those bits pick the result. `centered_pds_shots.py` shows the same effect across
term orders in-repo (§6).

- **Pin.** The gated study runs in a spawned child with PYTHONHASHSEED=0 and BLAS at 1 thread
  (`pinned_n2_study`). Every number is then reproducible on a given platform. The draft pinned
  threads only.
- **G1 non-vacuity, measured as written.** The draft's code took "shot-noise error" as |PDS(K) − E₀|.
  At K ≤ 3 that is mostly PDS truncation error (PDS(1) = ⟨H⟩), so the check could not fail. It now
  measures |PDS_noisy(K) − PDS_noiseless(K)|, which is what the text above says. The threshold is
  unchanged.
- **G2b — the verdict is not a float64 lottery.** At every sampled budget, the G2 verdict is the
  same on all 9 float64 paths: H's own term order plus 8 seeded re-orderings of the Pauli sum
  (`PATHS`). Every path uses the same snapshots. If it fails, the verdict is path-dependent and
  is reported that way.
- **Same estimator, faster.** Snapshots are grouped by their (basis, sign) pattern on each half of
  the qubits, so the sector block is Bhᵀ·C·Bl with C a pattern-count matrix. A one-off probe on
  2·10⁴ N₂ snapshots matched the draft's per-snapshot GEMM to 8e-15 in the weights. The whole
  16-trial, 9-path study takes 51 s on one M3 core. G0b still checks the estimator against
  `classical_shadows` on H₂.
- The 64-trial figure is not used for any gate. Any run beyond 16 trials (seeds 2000..) is reported
  as a lead and cannot change the verdict.

## 6. Measured data

All numbers below come from `uv run python centered_pds_shots.py 16` (pinned path, 16 trials,
seeds 2000..2015, 51 s), unless a line says otherwise. Run on an Apple M3 under macOS 27 with
Accelerate BLAS, 2026-10-02. hits = trials out of 16 with |PDS(K) − E₀| < 1.6 mHa. "paths" = the
min–max over the 9 float64 paths. E₀ = −107.62310177 Ha, μ = −9.0032 Ha.

**Gate runs** (`uv run pytest -q tests/test_centered_pds_shots_spec.py`). With G2 and G2b as
registered: **5 passed, 2 failed**. G2 failed with `assert True is False`: the rule says survives
at 10³. G2b failed with `assert 2 == 1`: two different verdicts across paths at 10³. G0a, G0b,
G0c, G1 and G3 passed. G2 and G2b were then rewritten to assert the kills in a form that holds
on any BLAS build. G2: the pre-registered pair fails on at least one float64 path. G2b: the
verdict splits across paths at some budget, while centered hit counts never move. The
pinned-path verdicts are recorded below, not asserted. Result: **7 passed** (58 s).

**Noiseless (exact moments, same sector code).** Centered PDS(K) − E₀ for K=4..8 is 1.076, 0.2602,
0.05401, 0.005873 and 0.001056 mHa, which matches `SPEC_centered_pds` (G0a). Raw PDS(K=6, 7, 8) on
the 9 paths, in mHa:

| path | K=6 | K=7 | K=8 |
|---|---|---|---|
| pinned (H's own order) | 0.05048 | **−24 290** | 0.05236 |
| 1 / 2 / 4 / 5 | 0.049–0.058 | 0.035–0.050 | −0.0174 … 0.050 |
| 3 | 0.05352 | 0.0415 | **−224** |
| 6 | 0.04936 | **−42 900** | 0.01719 |
| 7 / 8 | 0.047–0.057 | 0.043–0.047 | −0.027 / −0.0033 |

The raw frame is a last-bits lottery before any shot noise is added. On these 9 paths, 3 have a
catastrophic K=7 or K=8 root, and 3 more are slightly non-variational at K=8 (−0.003 to −0.027
mHa).

**Shared-shot noise, K=4..8.** Median |error| is in mHa on the pinned path. Hits are given as
raw (range over paths) / centered. Centered hits are identical on all 9 paths everywhere.

| S | K | median raw | median cen | hits raw (paths) / cen | indep. median cen |
|---|---|---|---|---|---|
| 10³ | 4 | 220.9 | 220.9 | 0 / 0 | 499 |
| | 5 | 24.71 | 24.71 | 0 / 0 | 263 |
| | 6 | 22.0 | 21.95 | 0 / 0 | 409 |
| | 7 | 5.614 | 4.938 | 1 (1–3) / **5** | 218 |
| | 8 | 4.523 | 2.616 | 5 (3–5) / 6 | 302 |
| 10⁴ | 4 | 13.70 | 13.70 | 0 / 0 | 212 |
| | 5 | 10.46 | 10.46 | 0 / 0 | 60.2 |
| | 6 | 4.522 | 4.512 | 3 / 3 | 80.5 |
| | 7 | 6.919 | 9.071 | 3 (1–4) / 3 | 151 |
| | 8 | 8.435 | 1.042 | 2 (1–6) / **11** | 51.2 |
| 10⁵ | 4 | 10.55 | 10.55 | 2 / 2 | 79.8 |
| | 5 | 7.673 | 7.673 | 1 / 1 | 52.8 |
| | 6 | 2.620 | 2.629 | 6 / 6 | 69.1 |
| | 7 | 1.857 | 0.8824 | 8 (4–11) / 10 | 55.0 |
| | 8 | 1.319 | 0.3613 | 10 (9–12) / 11 | 59.3 |
| 10⁶ | 4 | 2.390 | 2.390 | 3 / 3 | 40.9 |
| | 5 | 1.629 | 1.629 | 8 / 8 | 325 |
| | 6 | 0.4545 | 0.4550 | 13 / 13 | 37.9 |
| | 7 | 1.083 | 0.2260 | 9 (8–13) / **13** | 39.6 |
| | 8 | 0.3281 | 0.08869 | 14 (8–14) / 14 | 30.2 |

At K ≤ 6, raw hits are the same on every path and equal centered hits. Under the independent
per-moment contrast, no frame reaches more than 1/16 hits at any sampled budget and K.

**G1 (PASS).** Where raw cond(M) < 1e10, max |PDS_raw − PDS_cen| = 1.48e-8 Ha. Two pairs (10³,
trials 4 and 8, K=2) have no real root in either frame, and no pair disagrees on whether a root
exists. The minimum median shot-noise deviation at K ≤ 3 is 8.2e-3 Ha, 82× the 1e-4 non-vacuity
floor. ⟨H⟩ itself is very noisy under random-Pauli shadows: the median |⟨H⟩ − ⟨H⟩_exact| is
699.6 / 410.5 / 248.2 / 135.0 mHa at 10³ / 10⁴ / 10⁵ / 10⁶. Two
G1 assertion fixes were made **after** the first study run (thresholds unchanged). (i) A pair with
no real root in both frames now counts as agreement; the draft's `diff < 1e-6` failed on those 2
NaN pairs. (ii) The draft's guard "every trial is compared at K=1..3" died at 10⁶ on trial 7, K=3,
where raw cond(M) ≥ 1e10 (a near-singular Hankel matrix made by the noise). The guard is now
"K=1..3 are compared at every budget".

**G2 (KILLED).** The rule's verdict on the pinned path at 10³ / 10⁴ / 10⁵ / 10⁶ is survives /
survives / not / survives. The pre-registered "does not survive at 10³" is false. Both
pre-registered budgets pass the rule at exactly its margin, at K=7: 10³ has 5 vs 1 of 16 and 10⁶
has 13 vs 9 of 16. The draft docstring's unpinned "died at 10⁶" was another draw of the same
lottery: raw K=7–8 hits at 10⁶ span 8–14 across paths.

**G2b (KILLED).** Surviving paths, out of 9: 10³ 3 (`TFTTFFFFF`), 10⁴ 9, 10⁵ 2 (`FFFFFTFFT`), 10⁶ 4
(`TTTFTFFFF`). The pre-registered pair (survives at 10⁶, not at 10³) holds on 2 of 9 paths. 10⁴ is
the only budget with a path-independent verdict (survives: centered 11/16 at K=8 against raw 1–6).

**G3 (PASS).** At 10⁶ and K=8, centered median |error| is 30.23 mHa with independent per-moment
noise and 0.08869 mHa with shared shots: 341× (gate ≥ 10×).

**Scaled rows (10⁸–10¹², lead, not gated).** The rule says "not survives" on every path. Centered
PDS(4) already hits 14–16/16, because the noiseless PDS(4) error of 1.076 mHa is below 1.6 mHa, so
"going past K=4 helps by 25 points" is impossible. This is the rule's blind spot: centering's
noiseless gain at K ≥ 6 is ≤ 0.05 mHa, and a 1.6 mHa hit rate sees it only through raw's
catastrophic roots. Raw K=7–8 still misses on every scaled budget. At 10¹² on the pinned path, raw
hits 13/16 at both K=7 and K=8 (paths 11–14 and 8–15) while centered hits 16/16. Centered K=8 at 10⁸ hits
only 10/16, below the 14/16 at sampled 10⁶. The rescaled rows do not behave monotonically, which is
one more reason they stay a lead.

**Verdict (chem-u87).** (1) Centering does not change how shot noise amplifies into PDS(K). With
one shared shot record, the two frames are the same function of the same data (G1). (2) Whether
centering's float64 gain survives shot noise cannot be settled by the pre-registered rule at 16
trials. Raw's K=7–8 hit counts swing by up to 7/16 across float64 paths (10⁵, K=7: 4–11), and the
rule's 4/16 margin sits inside that swing (G2, G2b KILLED). What survives on every path is
reproducibility: centered hit counts are identical on all 9 paths at every budget and K. On the
pinned path at 10⁶, the median error at K=7 / K=8 is 1.083 / 0.3281 mHa raw against 0.226 /
0.0887 mHa centered (post hoc). The 64-trial lead below says the gain is real at K=7 but about the
size of the rule's margin. (3) The independent per-moment model behind
`SPEC_centered_pds_shot_noise` overstates the noise by 341× (G3). Shared-shot noise is a signed
measure on the true spectrum, so high-K PDS stays anchored to E₀. At 10⁶ snapshots ⟨H⟩ is off by a
median 135 mHa, yet PDS(8) lands within chemical accuracy in 14/16 trials.

**64 trials (lead, not gated).** From `uv run python centered_pds_shots.py 64` (seeds 2000..2063,
244 s). It cannot change the verdict above. Surviving paths out of 9: 10³ 0, 10⁴ 7 (`TFTFTTTTT`),
10⁵ 0, 10⁶ 6 (`TTFTFTTFT`). The pre-registered pair now holds on 6 of 9 paths, but 10⁴ and 10⁶ are
still path-dependent. At 10⁶, K=7, centered hits 56/64 against raw 35–46/64 across paths: centered
is ahead on every path, by 10–21 hits against the rule's 16-hit margin. At 10⁶, K=8, raw draws
level: 39–54 against centered 53, and the pinned path's raw 54 is one ahead. At 10⁵ centered leads
on every path at K=7 (44 vs 31–38) and K=8 (55 vs 41–47), but by less than 16. G1 holds at 64
trials: max diff 3.47e-8 Ha, 10 no-root pairs, 0 disagreeing. G3 is 153× (12.81 vs 0.0835 mHa).

## 7. Out of scope

- Grouped/derandomized shadows, Hadamard-test moments, gate noise and readout error.
- Systems without symmetry projection, where noise can see lower-lying states in other sectors
  (§2 d).
- Changing `moment_expansion` defaults.
