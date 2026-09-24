# SPEC: Adaptive Non-Uniform Shot Allocation for ODMD — and where it actually pays

**Status:** v2.0 — IMPLEMENTED, headline hypothesis **FALSIFIED**, residual claim gated.
v1.0 ("beat the visibility law 5–10× with a decaying schedule") is preserved in §8 with the
measurement that killed it.
**Depends on:** `odmd.py`, `device_odmd.py`, `visibility_law.py`
**Gates:** `tests/test_adaptive_shots_spec.py` (G1, G2, G2b, G3, G4) · `adaptive_shots.py` ·
`data/adaptive_shots_bench.csv`

---

## 1. Goal

Ask whether a fixed measurement budget `S_total`, spent **non-uniformly** across the ODMD survival
amplitudes `s_k = ⟨φ|e^{−ikτH}|φ⟩`, resolves the ground state more cheaply than the uniform
allocation the visibility law prices (`shots* ∼ 1/(w²K)`). `visibility_law.py`'s own honest-scope
note names this as an open hypothesis — *"uniform shot allocation (adaptive schemes could beat it
— a hypothesis, not a bug)"*. This spec closes it.

**Answer:** no, not in the regime v1.0 assumed — and the reason is structural, not a tuning
failure. Adaptive allocation pays **only** when the signal is damped, and then by ~1.3× in error
(~1.8× in budget), not 5–10×.

## 2. What was assumed, and what is actually true

v1.0's core insight was: *"Early time-steps have high amplitudes and carry massive structural
information, whereas late steps are buried in noise."*

**Measured (G2):** for a closed system the survival amplitude is quasi-periodic, not decaying. On
N₂ CAS(6,6) at 1.1 Å, `|s_k| ≥ 0.85` for every `k` out to 23 — `[1.00, 0.958, 0.900, 0.893, 0.914,
0.924, 0.917, 0.919, …]`. There is no noise-dominated late tail to defund.

**The consequence is a theorem, not a measurement.** With every element carrying equal signal, the
total noise power under a fixed budget is

    Σ_k σ_k² = V · Σ_k 1/S_k ,   V = 2(2 − 1/dim),   Σ_k S_k = S_total

which is **convex** in the allocation, so uniform `S_k = S_total/K` is the constrained minimiser.
Any decaying schedule strictly increases it. Adaptive allocation *cannot* beat uniform here, and
G2 gates that as a ceiling on the achievable gain rather than trusting the argument.

**Where a dead tail does exist:** under the repo's validated global-depolarizing model
`s_k → f^k s_k` (`device_odmd.py`), which damps the amplitude while leaving every exact eigenphase
invariant. That is the regime this spec retargets, and it is the physically relevant one — it is
what a real device does to the signal.

## 3. Approach (as built)

1. **Schedules.** `S_k = C e^{−αk}` and `S_k = C (k+1)^{−γ}` for `k ≥ 1`, with `S_0 = 0` (`s_0 = 1`
   exactly — spending shots on it is pure waste, and allocating a `k=0` share would silently
   unbalance the matched-budget comparison). Negative `α` (growing schedules) is in range, so
   "uniform is optimal" is a reachable answer rather than an assumption.
2. **Per-element noise.** Element `k` gets its own Hadamard-test scale `σ_k² = V/S_k` (odmd.py
   conventions, `s_0` exact).
3. **Re-weighting — and why the exponential family is the specced one.** Heteroskedastic elements
   are whitened by `r_k ∝ 1/σ_k ∝ √S_k`. For a **geometric** schedule that weight is itself
   geometric, `r_k = c·β^k` — precisely the transformation that multiplies every DMD eigenvalue by
   `β` and leaves its **phase** untouched (the same algebra as ODMD's depolarizing immunity). So
   whitening an exponential schedule costs nothing: homoskedastic noise, **zero** energy bias, no
   new DMD code (`device_odmd_energy`'s wide modulus window already handles `|λ| < 1`).
   A polynomial schedule's weights are not geometric; the same whitening shifts the noiseless
   energy by tens of mHa. **G2b gates both halves** — bit-identical energies across `α` (residual
   `< 1e-7` Ha, SVD conditioning, not growing with `α`) versus `> 10` mHa for polynomial.

## 4. Public interface (`adaptive_shots.py`)

```python
exponential_schedule(K, total_shots, alpha) -> np.ndarray      # S_0 = 0, sum == total_shots
polynomial_schedule(K, total_shots, gamma)  -> np.ndarray
whitening_weights(schedule)                 -> np.ndarray      # r_k ∝ √S_k, r_0 extrapolated
sample_adaptive_odmd_energy(prob, schedule, seed, damping=1.0, ...) -> float
median_error(prob, schedule, seeds=200, damping=1.0, ...)      -> float
optimize_decay_factor(prob, K, total_shots, damping=1.0, ...)  -> float
```

**Two deviations from v1.0's signatures, both deliberate:**

- The samplers take an **`ODMDProblem`**, not a `MolecularHamiltonian` + `tau`. v1.0's
  `sample_adaptive_odmd_energy(mh, K, schedule, tau=0.5)` would rebuild the Hamiltonian, diagonalise
  it and re-propagate the signal on every one of the 1000 noise realisations G4 needs. Same
  convention as `odmd.sample_odmd_energy`; `tau` and the exact reference come from the problem.
- `optimize_decay_factor` is an **offline planner**, stated as such: it scores against the known
  exact `prob.ref`, so it sizes a schedule in simulation *before* an experiment. It is not a
  runtime estimator, and nothing in the gates treats it as one.

## 5. Acceptance gates (as gated, with the measured numbers)

- **G1 — Strict budget preservation.** `|Σ S_k − S_total| ≤ 1.0` and `S_0 = 0`, `S_k > 0` for
  `k ≥ 1`, across `K ∈ {8,12,16,24}` and both families. *Passes.* This is what makes G2/G4
  matched-budget comparisons.
- **G2 — THE KILL (replaces v1.0's G2).** Undamped N₂ CAS(6,6), `K=12`, `S_total=1e4`:
  `|s_k| ≥ 0.85` everywhere; uniform sits at **5.4 mHa** (v1.0 asserted `> 8`); and of ten
  decaying schedules (`α ∈ [0.05, 0.8]`, `γ ∈ [0.25, 2]`) **none** reaches the 1.0 mHa target and
  **none** beats uniform — best gain **0.996×**, worst 0.13×. *Passes as a ceiling
  (`gain < 1.05`).*
- **G2b — Geometric whitening preserves eigenphases.** Noiseless, `K ∈ {12,16,20}`: exponential
  whitening leaves the energy invariant across `α ∈ {0.1, 0.15, 0.3, 0.5}` to `< 1e-7` Ha;
  polynomial whitening shifts it by `> 10` mHa. *Passes.*
- **G3 — THE SURVIVOR (v1.0's G3, sharpened).** `α*` from a grid over `[−0.10, 0.50]`, `K=24`,
  `S_total=1e4`, damping `f ∈ {1.0, 0.9, 0.8}`: `α* = 0.00, +0.15, +0.15` — non-decreasing, `≤ 0`
  at the undamped end, `≥ 0.10` at both damped ends. Stable across 150–600 seeds. *Passes.*
- **G4 — The size of the win, replacing v1.0's 5× variance claim.** 1000 realisations, `K=24`,
  matched `S_total=1e4`: damped (`f=0.8`) median error improves **1.35×** (≈1.8× in budget at the
  `1/√S` scaling) — gated `> 1.25×`. Error **variance** improves **1.07×**, not `≥ 5×` — gated
  `< 1.5×` *as the falsification*. Undamped, the same schedule is strictly **worse** (0.72×).
  *Passes.*

## 6. The findings

1. **Uniform allocation is optimal for undamped ODMD, by convexity.** The visibility law's
   `1/(w²K)` is not beatable by re-allocating a fixed budget; it is the optimum, not an artefact of
   a uniform convention. `visibility_law.py`'s open hypothesis is closed **negative**.
2. **The whole advantage of non-uniform allocation is the damping rate.** `α*` is 0 when `f = 1`
   and positive as soon as `f < 1`, rising with `−log f`. The mechanism is defunding a tail that
   genuinely carries no signal — not "concentrating on high-information early steps", which was the
   v1.0 story and is wrong.
3. **`α*` saturates.** It rises `0 → 0.15` between `f = 1.0` and `f = 0.9`, then stays ≈0.15 down
   to `f = 0.7` (see `data/adaptive_shots_bench.csv`). Deeper damping does not keep buying steeper
   schedules: the optimum is a compromise between defunding the tail and preserving early-point
   precision, and the second term stops it. G3 gates monotonicity, **not** proportionality — the
   proportional version would fail.
4. **The win never exceeds ~1.35× in error (~1.8× in budget)**, and grows with depth (1.10× at
   `K=12`, 1.17× at `K=16`, 1.35× at `K=24`, all at `f = 0.8`): the deeper the window, the more
   dead tail there is to defund. Extrapolating that trend, not tuning `α`, is the only route to a
   larger number.
5. **The error variance is not a shot-allocation quantity.** Median and p90 errors both improve
   ~1.35×, but the variance moves 1.07× — it is dominated by rare catastrophic mode
   misidentifications, which a budget reshuffle does not touch. Any future scheme targeting v1.0's
   variance claim has to attack mode selection, not allocation.

## 7. Honest scope and caveats

- Exact-statevector `s_k` with the idealised i.i.d. Hadamard-test noise of `odmd.py`. Damping enters
  as the **global-depolarizing model**, not measured device noise; local gate noise is not a global
  channel (`device_odmd.py`), so `α*` on hardware must be measured, not assumed.
- One system (N₂ CAS(6,6)) and one budget (`S_total = 1e4`). The convexity argument behind G2 is
  system-independent; the `α*` values in G3 and the gains in G4 are not.
- The geometric rescale needs a **meaningful SVD truncation** to be safe. In the `σ → 0` limit at
  `K = 24` the noise-edge cutoff collapses and retained numerical-noise modes can win the
  eigenphase argmin (an artefact of the noiseless probe, not of the noisy path — every gated
  measurement here has a real `σ`). G2b is therefore gated at `K ≤ 20`.
- `optimize_decay_factor` needs the exact reference. It plans budgets; it does not estimate
  energies.
- ODMD remains **non-variational** (`SPEC_odmd.md` G2) — errors are two-sided, which is why every
  gate here uses `|E − E_exact|`.

## 8. v1.0, and the measurement that killed it

> **G2 (v1.0):** *"Under simulated noise with a tight budget `S_total = 10⁴` and `K=12` on N₂
> CAS(6,6), the optimal decaying schedule must achieve a median error `< 1.0 mHa`, while the
> standard uniform schedule yields an error `> 8.0 mHa`."*
> **G4 (v1.0):** *"Across 1,000 independent noise realizations, the adaptive scheduler must reduce
> the ground-state error variance by at least 5× compared to the uniform budget."*

Measured at exactly those settings: uniform **5.4 mHa**; best decaying schedule **5.39 mHa**
(γ=0.25, i.e. barely distinguishable from uniform); every steeper schedule worse, down to 42 mHa at
γ=2. Variance ratio at `K=24, f=0.8`: **1.07×**. Neither gate is reachable, at any `α` or `γ`, and
§2 says why it was never going to be. The spec was revised rather than the tolerances loosened —
`specs/README.md` step 5.
