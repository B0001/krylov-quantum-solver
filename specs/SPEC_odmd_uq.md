# SPEC: Coverage-gated error bars for ODMD — a union bootstrap from a single signal

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_odmd_uq_spec.py`).

---

## 1. Goal

Every ODMD result in this repo so far reports median error over ~100 noise seeds *against a
known FCI answer*. A real experiment gets **one noisy signal and no ground truth**. Claim: a
resampling scheme on that single signal — the **union** of a parametric bootstrap (refit modes →
rebuild the clean signal → re-noise at the *known* shot-noise scale σ) and a BOP-DMD-style
bagging ensemble (random Hankel-column subsets) — yields a confidence interval whose *empirical
coverage meets its nominal level*. Error bars are falsifiable in exactly one way, and that is
the gate: over 200 independent noise realizations per configuration, the nominal-90% interval
must cover the true energy at least 85% of the time, on every system and budget tested.

## 2. Background, the found failure, and honest framing

- Prior art: bagging/ensemble UQ for DMD is BOP-DMD (Sashidhar & Kutz, `arXiv:2107.10878`);
  parametric bootstrap is textbook. The composition for quantum eigenphase estimation with a
  known per-element noise scale, and the coverage-gated validation, are what is new here.
- **Found during probing (the real finding):** each component alone is *broken*, in
  complementary regimes. The parametric bootstrap is **anti-conservative by up to 18×**
  (coverage 0.05 at nominal 0.90 on N₂ CAS(6,6), 10⁵ shots): the DMD fit absorbs the realized
  noise into its modes, so re-noising the fitted signal produces an artificially stable ensemble
  that never sees the threshold rank-switching instability — and bootstrap mean-shift bias
  correction does *not* fix it (probed: 0.03). Bagging sees that instability but under-spreads
  few-mode signals (0.785 on H₂). The union covers both failure modes.
- **What we can claim if gates pass:** a single-signal 90% CI with measured coverage ≥ 0.85
  (observed 0.895–1.000) at a conservatism cost of ~2–3× the median error, across H₂/H₄/N₂ and
  10⁴–10⁵ shots.
- **What we cannot claim:** the interval quantifies **variance under the fitted model, not
  bias** — a too-shallow signal (N₂ at K=8: 8.7 mHa truncation bias) is invisible to any
  resampling of that signal, and coverage collapses to 0 (gated as the boundary). Always pair
  with a K-convergence check. Idealized i.i.d. Gaussian per-element noise with *known* σ (the
  repo's shot-noise conventions); one nominal level (α=0.1) validated; ODMD remains
  non-variational.

## 3. Approach

From one noisy signal `s` (length K, per-element noise scale σ known from the shot budget):

1. **Point estimate:** the pinned `odmd_energy(s, tau, svd_threshold=5σ)`.
2. **Parametric arm:** fit modes (`_dmd_modes`, same threshold) → Vandermonde least-squares
   amplitudes → clean reconstruction `s_hat` (s₀ forced to 1) → `n` resamples
   `s* = s_hat + CN(0, σ)` → re-estimate → percentile interval.
3. **Bagging arm:** `n` random subsets of 60% of the Hankel columns (the DMD least-squares is
   over columns, so `X'_sub = A X_sub` stays exact) → truncated-SVD DMD per bag → percentile
   interval, recentered at the point estimate.
4. **Union:** `[min(lows), max(highs)]`.

**Reference for validation:** exact centered ground energy; coverage counted over 200
independent outer noise realizations (seeded — the gates are deterministic).

## 4. Public interface

```
odmd_uq.ODMDInterval        # dataclass: estimate, lower, upper, half_width,
                            #   parametric + bagging component intervals, sigma, alpha
odmd_uq.odmd_confidence_interval(s, tau, sigma, n_resamples=200, alpha=0.1, seed=0)
    -> ODMDInterval
```

## 5. Acceptance criteria (validation gates)

`tests/test_odmd_uq_spec.py`; 200 trials × (200+200) resamples per configuration; nominal 90%.

- **G1 — coverage, no exceptions (DEFINITION OF DONE).** Union coverage ≥ 0.85 on ALL of:
  H₂ × {10⁴, 10⁵}, H₄ × {10⁴, 10⁵}, N₂ CAS(6,6) × {10⁴, 10⁵} shots at K=24 (measured
  0.910/0.910/0.895/0.975/0.965/1.000).
- **G2 — informative width.** Median half-width between 1× and 4× the median true |error| on
  every configuration (measured 1.9–2.9× — conservative, not vacuous).
- **G3 — the union is load-bearing.** In the same trials, the parametric arm alone fails hard
  on N₂ 10⁵ (coverage < 0.5; measured 0.05) and the bagging arm alone dips below gate level on
  H₂ 10⁵ (< 0.85; measured 0.785) — neither component suffices; the recorded finding.
- **G4 — bias is invisible (the boundary).** N₂ at K=8, 10⁶ shots (truncation bias ≫ noise):
  union coverage < 0.1 (measured 0.000). Error bars from resampling one signal CANNOT detect
  model misspecification — pair with a depth-convergence check before trusting any interval.

## 6. Implementation plan (test-first)

1. `tests/test_odmd_uq_spec.py` encoding G1–G4 (RED — `odmd_uq` missing).
2. `odmd_uq.py` — reuses `odmd._dmd_modes` / `odmd_energy`; the column-subset DMD is the only
   new numerics (~15 lines).
3. `make gates`.

## 7. Out of scope

- Device-noise signals (Aer damping changes the noise model; compose with `SPEC_device_odmd`
  later), excited-state/gap intervals, other nominal levels, studentized/double bootstrap.
- Bias detection (fundamentally outside single-signal resampling — see G4).

## 8. Caveats and risks

- **R1 — conservatism is the price of validity:** ~2–3× wider than the true spread. An
  anti-conservative error bar is worse than none; this trade is deliberate.
- σ must be known (shot budget). A misestimated σ silently miscalibrates the parametric arm.
- Coverage validated at K=24, three systems, Gaussian noise: extrapolate beyond with care.

## 9. Deliverables

- `odmd_uq.py`; `tests/test_odmd_uq_spec.py` — gates G1–G4.
- `BACKLOG.md` entry with the measured coverage table and the two-sided failure finding.
