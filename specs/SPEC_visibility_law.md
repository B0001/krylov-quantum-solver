# SPEC: The visibility law, made predictive — a calibrated shot-cost law for spectral lines

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_visibility_law_spec.py`).

---

## 1. Goal

Three specs recorded the same qualitative rule (`SPEC_odmd_excited`: excited modes need
p·√(dm) above the noise edge; `SPEC_device_odmd`: immunity ends at the shot floor;
`SPEC_odmd_spectral`: weak satellites need deeper K). Claim: that rule is a **quantitative,
transferable law**. For the physical (unnormalized) correlation signal
`C_k = ⟨ψ₀|Ô e^{−ikτH} Ô|ψ₀⟩`, a line of weight w is detectable iff its Hankel singular value
clears the noise edge, giving

    σ* = w·√(dm) / (c·(√d + √m))     ⇒     shots* = 2(2−1/dim)/σ*²  ∝  1/(w²·K).

Falsifiable three ways: the 50%-detection crossover must scale as **w⁻²** (a straight −2 line
over four orders of magnitude in w = eight orders in shots), a **single calibration** of the
prefactor on one line must predict every other line's crossover — including lines inside a
*multi-line* signal — and the crossover must move as **1/K** with signal depth.

## 2. Background and honest framing

- The ingredients are standard (random-matrix detection thresholds; `SPEC_odmd_excited`'s noise
  edge). What is new here is the *predictive packaging*: one measured constant → the shot budget
  for any target line of known weight, gated by transfer.
- **What we can claim if gates pass:** an experiment-planning tool — given a line's expected
  weight, `predicted_shots` returns the budget to resolve it (e.g. the Nb₃F₈ near-dark optical
  line costs ~2·10⁸ shots/element where Nb₃I₈'s costs ~2) — validated to factor ~1.1 across
  every line we can check, single- and multi-line.
- **What we cannot claim:** Gaussian per-element noise with known σ (the repo's conventions;
  device damping shrinks the *effective* w — compose with `SPEC_device_odmd` before budgeting
  hardware); detection ≠ accuracy (the law prices *seeing* a line, not resolving it to a target
  precision); the prefactor c is calibrated, not derived (Hankel noise is correlated — the same
  caveat as the noise edge itself); line *attribution* needs a tolerance below the line spacing,
  or strong lines masquerade as weak neighbors (found in probing, gated as the boundary).

## 3. Approach

Detection protocol per trial: add CN(0, σ) noise to C, form the Hankel pair, and declare the
line detected iff (i) the top singular value clears `c·σ(√d+√m)` (c = 1.2, the pinned edge) and
(ii) a retained DMD pole lies within `tol_frac·π/τ` of the line (default tol_frac = 0.03, below
every line spacing used here). Crossover = 50% detection rate over 60 seeded trials, found by
bisection in log-shots. References: the law itself (slope, transfer, depth scaling are all
overdetermined checks against measured crossovers); line weights from the pinned
`SPEC_odmd_optical`/`SPEC_odmd_spectral` machinery.

## 4. Public interface

```
visibility_law.detect_line(C, tau, sigma, e_line, seed, c=1.2, tol_frac=0.03) -> bool
visibility_law.detection_rate(C, tau, shots, e_line, dim, seeds=60, **kw) -> float
visibility_law.crossover_shots(C, tau, e_line, dim, seeds=60, **kw) -> float
visibility_law.predicted_shots(w, K, dim, c=1.2, calibration=1.0) -> float
```

## 5. Acceptance criteria (validation gates)

`tests/test_visibility_law_spec.py`. Signals: the four dimer optical correlators (K=16, single
line each, w = 1.1e-4 … 0.96) and the H₄ so=0 removal correlator (K=32, three lines
0.77/0.17/0.027 in ONE signal).

- **G1 — the −2 power law.** Measured crossovers span ≥ 6 orders of magnitude in shots and the
  log-log slope of shots* vs w is in [−2.05, −1.95] (measured −2.000).
- **G2 — one calibration predicts everything (DEFINITION OF DONE).** Calibrate the prefactor on
  Nb₃Br₈ alone; the predicted crossovers of Nb₃F₈/Cl₈/I₈ AND of all three H₄ multi-line
  components are each within a factor 1.5 of measured (measured 1.00/1.00/1.00 and
  0.86/1.11/1.09).
- **G3 — depth buys shots linearly (gate revised during implementation).** Reality showed two
  separable effects, so the gate was split (the SDD-honest revision): the **edge** component
  scales as 1/K — on Nb₃I₈, shots*(K=8)/shots*(K=32) ∈ [3, 5] with attribution wide open
  (measured 3.83; law: 4) — while the full protocol's *tight attribution* adds a pole-accuracy
  cost concentrated at shallow depth (> 30% overhead at K=8, measured 61%; < 20% at K=32,
  measured 4%). Budget with the full protocol at your working K, not the bare law, when K is
  small.
- **G4 — protocol soundness and the attribution boundary.** False-positive rate at w = 0 is
  ≤ 2% (measured 0.5%); and with a sloppy tolerance (tol_frac = 0.1 ≥ the 0.215 Ha line
  spacing) the H₄ middle line's apparent crossover falls > 10× below the law (measured 45×) —
  the strong neighbor masquerades. Attribution tolerance must sit below the line spacing.

## 6. Implementation plan (test-first)

1. `tests/test_visibility_law_spec.py` encoding G1–G4 (RED — module missing).
2. `visibility_law.py` — the probe-validated protocol, verbatim.
3. `make gates`.

## 7. Out of scope

- Device-damped signals (fold the measured damping into w first — `SPEC_device_odmd`);
  precision budgeting beyond detection; optimal (non-uniform) time sampling; deriving c from
  Hankel random-matrix theory.

## 8. Caveats and risks

- **R1:** the law prices detection at KNOWN line position/weight class; blind surveys need the
  false-positive side of the protocol (G4 bounds it at the default settings only).
- Budgets assume every element gets equal shots; adaptive allocation could beat the law
  (a follow-up hypothesis, not a bug).

## 9. Deliverables

- `visibility_law.py` (+ `__main__`: the shot-budget table for the dimer optical lines).
- `tests/test_visibility_law_spec.py` — gates G1–G4.
- `BACKLOG.md` entry.
