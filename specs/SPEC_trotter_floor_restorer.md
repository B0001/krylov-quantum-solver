# SPEC: Signal-domain Richardson for Trotter ODMD — the leak that isn't there

**Status:** v2.0 — IMPLEMENTED. **G1 FALSIFIED** (there is no variational-floor leak); **G3 is out
of scope** (no Lehmann bracket exists on this path); **G4 reversed** into a regime map; G2
re-gated. v1.0 is preserved verbatim in §8 with what resolved each claim.
**Depends on:** `trotter_odmd.py`, `odmd.py`
**Gates:** `tests/test_trotter_restorer_spec.py` (G1, G1b, G2, G3, G4, G5 — 6 tests, 24 s) ·
`trotter_restorer.py` · `data/trotter_restorer_bench.csv`

---

## 1. Goal

v1.0's goal was to plug a **"Trotter variational-floor leak"** — the claim that product-formula
discretization shifts the effective spectrum so a Trotterized Krylov/ODMD solve returns
E < E_FCI, destroying certified bounds — using Richardson extrapolation applied to the complex
survival signal *before* diagonalization.

**Answer: the leak does not exist, and the premise inverts.** Across 2 systems × reps {1,2,4} × 6
window lengths, the order-2 Suzuki Trotter ODMD energy is **above E_FCI in 36/36 configurations**
(minimum excess +0.26 mHa). Meanwhile the *extrapolation* — the proposed cure — lands **below**
E_FCI in 12 of 24 configurations here (22 of 24 over the 4-system benchmark). Trotterization does
not break the floor; mitigation is the only step in this pipeline that does.

## 2. What v1.0 assumed, and what is actually true

**The bias has the wrong sign for the story.** v1.0 asserts the ground state of
H_eff = H + Δt²·H_err "can easily drop below" E_FCI. Measured, the order-2 Suzuki bias here is
consistently **positive**: +19.6 / +4.2 / +1.0 mHa at reps 1/2/4 on stretched H₂, +10.8 / +2.6 /
+0.6 on the H₄ chain. Nothing drops. (This is about the *Trotter* bias specifically — ODMD remains
non-variational in general, per `SPEC_odmd.md` G2, and that is a separate and already-gated fact.)

**v1.0's extrapolation formula is not time-aligned.** It pairs `s_k(Δt)` with `s_k(2Δt)`. Those sit
at physical times *kΔt* and *2kΔt* — different points on the signal — so combining them at equal
index *k* subtracts unlike quantities and the O(Δt²) term cannot cancel. The correct refinement is
**`reps` at fixed τ**: every sample stays at *kτ* while dt_eff = τ/reps halves. This is what
`trotter_odmd.build_trotter_odmd_problem(mh, n, reps=r)` already provides, and G5 gates that the
three rep counts share a time axis.

**Signal-domain extrapolation has a validity window, energy-domain does not.** The Trotter bias
enters the signal as `exp(−i k τ δE)`, so a *linear* combination of two signals cancels it only
while that accumulated phase is small. Past the window the extrapolant degenerates back onto the
un-extrapolated fine signal. The eigenphase itself obeys the δE ∼ dt² law exactly, so the scalar
(energy-domain) Richardson already shipped in `trotter_odmd.richardson_energy` is **flat in K**.

## 3. Approach (as built)

1. **Time-aligned signals.** `trotter_problems(mh, n, reps=(1,2,4))` — cached, fixed τ.
2. **Signal-domain Richardson.** `signal_richardson(s_coarse, s_fine) = (w·s_fine − s_coarse)/(w−1)`,
   `w = ratio**order`, the direct analogue of the scalar formula.
3. **A diagnostic, not a promise.** `phase_budget(K, τ, bias) = K·τ·|bias|` reports the accumulated
   phase that decides whether step 2 can work at all.
4. **Three modes compared head to head** — `raw`, `signal_extrap`, `scalar_extrap`.

## 4. Public interface (`trotter_restorer.py`)

```python
trotter_problems(mh, n=32, reps=(1,2,4))                  -> dict[int, TrotterODMDProblem]
signal_richardson(s_coarse, s_fine, step_ratio=2, order=2) -> np.ndarray
phase_budget(K, tau, bias)                                 -> float   # radians
restored_energy(problems, K, mode="signal_extrap", pair=(2,4)) -> float
```

**Deviations from v1.0, all forced by §2:**

- No `generate_trotter_signal(mh, K, dt)`. Refinement is by **`reps` at fixed τ**, not by `dt`;
  a `dt` parameter would re-introduce the time-misalignment the spec was built on.
- No `certified_restored_energy(...) -> LehmannBracket`. That return type cannot be produced on
  this path — see G3. Shipping a `LehmannBracket` here would be a certificate in name only, and
  G3 gates its **absence**.
- `restored_energy` takes prebuilt `problems` rather than an `mh`, so the three rep counts are
  built once and shared; the gates evaluate 36 configurations off one build.

## 5. Acceptance gates (as gated, with the measured numbers)

- **G1 — THE KILL (replaces v1.0's G1).** Raw Trotter ODMD is above E_FCI in **36/36**
  configurations; minimum excess **+0.26 mHa**. v1.0 required a drop of **at least 1.0 mHa below**
  E_FCI. Gated as a floor on the excess. *Passes.*
- **G1b — THE INVERSION (NEW).** The extrapolated energies land below E_FCI in **12/24** here and
  **22/24** on the 4-system benchmark, by ~0.01 mHa. Richardson's residual is two-sided; the
  extrapolated energy is **not a bound**, and the gate exists so nobody later treats it as one.
  *Passes.*
- **G2 — What signal extrapolation actually buys (re-gated).** H₄ at K = 16: raw **+0.646 mHa** →
  **+0.007 mHa**, an **87×** cut, real and worth having. At K = 32 the extrapolant returns the raw
  value to within **0.003 mHa** — the window has closed. Both halves gated. *Passes.*
- **G3 — OUT OF SCOPE, not merely unmet (replaces v1.0's G3).** Lehmann needs ⟨u|H²|u⟩ for a
  physical Ritz **state**. ODMD's whole premise is that it works from the scalar overlap series
  s_k alone: it returns eigen*phases*, and no state is ever formed, so neither ⟨H⟩ nor ⟨H²⟩ of a
  Ritz vector exists to bound. v1.0's "100% sound oracle-enclosed brackets" is not a hard target,
  it is undefined on this path. The gate pins that the module exposes **no** bracket function.
  Certifying a Trotter solve needs the QKSD/state path (`excited_bounds.py`) — a different spec.
  *Passes.*
- **G4 — REVERSED into a regime map (replaces v1.0's G4).** Neither domain dominates, and a first
  draft of this gate asserting "scalar is never worse" was itself falsified by the data:
  **inside** the phase window signal wins (H₂ at K = 6: **0.0045 vs 0.0476 mHa**, 10×);
  **outside** it scalar wins by much more (K = 32: **~95×** on H₄, **~21×** on H₂), and scalar is
  flat in K to < 2e-5 mHa. The controlling quantity is the accumulated phase, not the
  extrapolation domain. *Passes.*
- **G5 — The refinement parameter is `reps`, not `dt` (NEW).** τ is identical across reps 1/2/4,
  so the signals share a time axis and the aligned pairing beats raw; `phase_budget` grows
  monotonically with the window. *Passes.*

## 6. Findings

1. **No variational-floor leak exists for order-2 Suzuki here.** The bias is positive in every
   configuration tested. The spec's reason to exist is not reproducible.
2. **The premise inverts: extrapolation is what breaks the floor.** Raw 0/24 below FCI; extrapolant
   22/24. Any pipeline wanting a variational guarantee should prefer the *un*-mitigated energy and
   pay the bias, or bound it properly — not extrapolate and call the result certified.
3. **Signal-domain and energy-domain Richardson split by accumulated phase.** K·τ·δE is the
   controlling quantity. Signal wins at small K, degenerates at large K; scalar is flat. The
   practical recommendation is the method already shipped.
4. **v1.0's extrapolation formula was time-misaligned** — `s_k(Δt)` and `s_k(2Δt)` are different
   physical times. Refining `reps` at fixed τ is the fix, and it costs nothing new.
5. **"Pre-diagonalization preserves the state" does not hold for ODMD at all.** There is no state
   to preserve: the method consumes a scalar series. The argument would apply to a QKSD/state
   pipeline, which is where a future spec should take it.

## 7. Honest scope and caveats

- **Statevector-simulated Trotter circuits** (no device noise, no native-gate transpilation),
  inherited from `trotter_odmd.py`. On hardware the bias competes with shot and device noise.
- **Order-2 Suzuki only, reps ∈ {1,2,4}, two systems, K ≤ 32.** The *sign* of the Trotter bias is
  what kills G1, and sign is a property of the formula and system — a different order, ordering,
  or a strongly non-commuting Hamiltonian could plausibly produce a negative bias. **This spec
  does not prove no leak can ever exist; it proves the specced one does not occur here.** That is
  the honest boundary, and it is the natural follow-up hypothesis.
- **The phase window is calibrated, not derived.** K·τ·δE ≈ 0.15 rad is where H₂ degenerates;
  H₄ degenerates by 0.076 rad, so the threshold is system-dependent and the gates assert the
  qualitative regimes rather than a universal constant.
- **No certification anywhere in this module** — see G3.
- ODMD remains **non-variational** in general (`SPEC_odmd.md` G2); G1 is a statement about the
  Trotter bias, not a claim that ODMD respects a floor.

## 8. v1.0, and what resolved each claim

> **G1 (v1.0):** *"…the unmitigated Trotter-Krylov ground-state energy must drop below the true
> physical ground state E_FCI by at least 1.0 mHa, successfully exposing and verifying the
> variational-floor leak."*
> **G2 (v1.0):** *"…pulling the estimated energy back above the true physical ground state
> (E ≥ E_FCI) and reducing the residual discretization error to < 0.1 mHa."*
> **G3 (v1.0):** *"…must yield 100% sound, oracle-enclosed brackets for both E₀ and E₁ across all
> step-sizes."*
> **G4 (v1.0):** *"The certified bracket width … must be at least 2× narrower than any attempt to
> construct error bounds on post-solve scalar_extrap metrics, proving the mathematical superiority
> of state-level extrapolation."*

Resolved: G1's leak never occurs (36/36 above FCI, min +0.26 mHa). G2's "back above E_FCI" is
therefore vacuous — and the extrapolant it names is the thing that goes *below*; its < 0.1 mHa
accuracy target is met at moderate K (0.007 mHa on H₄ at K = 16) and missed once the phase window
closes. G3 is undefined on a signal-only method. G4 is backwards inside the phase window and right
outside it for the wrong reason, so it becomes a regime map; "any attempt" was also unfalsifiable
as written. The spec was revised rather than the tolerances loosened — `specs/README.md` step 5.
