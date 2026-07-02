# SPEC: Device-noise ODMD — eigenphases are depolarizing-immune (until the noise edge)

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_device_odmd_spec.py`).

---

## 1. Goal

Take the ODMD rung from idealized Gaussian shot noise to **device noise**. The physics: a global
depolarizing channel damps the survival signal geometrically, `s_k -> f^k s_k`, which multiplies
every DMD eigenvalue by `f` but **leaves its phase untouched** — so the ODMD energy is *exactly*
invariant under uniform damping, while KQD's generalized eigenproblem on the same damped data has
no such protection. Realistic gate-level noise is *not* a global channel, so the falsifiable
questions are: (i) exact immunity at the channel level (machine precision, or the claim dies);
(ii) how much phase bias *local* depolarizing gate noise induces through real Hadamard-test
circuits (measured on qiskit-aer); (iii) where immunity ends — when damping pushes the signal
under the shot-noise floor, the `SPEC_odmd_excited.md` noise-edge law says visibility dies.

## 2. Background and honest framing

- Phase-based estimators' robustness to global depolarizing noise is known in the
  QCELS/ODMD-adjacent literature (damping rescales amplitudes, not eigenphases — cf.
  `arXiv:2306.01858`'s noise analysis); the quantified *local-vs-global* gap on real
  Hadamard-test circuits, against this repo's validated KQD arm, is what is new here.
- **What we can claim if gates pass:** sub-0.1 mHa eigenphases through ≥ 50% device-noise
  amplitude loss on genuinely-Trotterized, transpiled, ancilla-controlled circuits (the
  `hardware_krylov` stack, post-`SPEC_trotter_odmd`-fix), with KQD failing catastrophically on
  identical data; and a mapped boundary where the immunity ends.
- **What we cannot claim:** Aer depolarizing + readout error is still a *simulated* device (no
  coherent errors, crosstalk, drift — real hardware phase noise is adversarial in ways
  depolarizing is not); H₂-scale circuits; the energies are eigenphases of the *Trotterized*
  unitary (the reps=1 Trotter bias of +19.27 mHa is already measured and removable per
  `SPEC_trotter_odmd` — G4 references the circuit eigenphase to isolate the *noise-induced*
  bias); ODMD remains non-variational.

## 3. Approach

1. **Channel level (exact):** damp the pinned H₄ signal `s_k -> f^k s_k` (and the KQD Toeplitz
   H-row identically); compare ODMD vs KQD, noiseless and under 10⁵-shot noise.
2. **Damping-robust estimator (composition, no new DMD code):** `odmd_spectrum` with the
   noise-edge cutoff, a **wide modulus window** (the `odmd_energy` unit-modulus filter
   `||λ|−1| < 0.2` misidentifies damped signal modes at `f < 0.8` — under noise, spurious
   near-unimodular modes pass instead), and an amplitude floor from the Vandermonde refit.
3. **Aer end-to-end:** `HardwareKrylovSolver` (new public `measure_signal`) measures
   `s_k = S_0k` by ancilla Hadamard tests on transpiled controlled-Trotter circuits in the
   centered frame, under `build_depolarizing_noise_model`; reference = exact ground eigenphase
   of the *same* step circuit (via `trotter_odmd.build_trotter_odmd_problem`).

## 4. Public interface

```
device_odmd.centered_frame(mh) -> (mh_centered, tau, mu)   # dataclasses.replace, total E invariant
device_odmd.device_odmd_energy(s, tau, sigma, amp_floor=0.05, c=1.2) -> float
    # noise-edge cutoff + mod_window=2.0 + amplitude floor (damping-robust composition)
device_odmd.measure_survival_signal(mh, n, shots=None, noise_model=None, seed=None,
                                    trotter_order=2, trotter_reps=1) -> np.ndarray
HardwareKrylovSolver.measure_signal(n) -> np.ndarray        # first overlap row S_0k (new, public)
```

## 5. Acceptance criteria (validation gates)

`tests/test_device_odmd_spec.py`. Channel gates on H₄ (n=24); Aer gates on H₂ (K=8, 32768
shots/observable, median over 5 seeds).

- **G1 — exact channel immunity.** At f ∈ {0.9, 0.7}: `|E_odmd − E₀| < 1e-6 Ha` (measured
  ~1e-8), while KQD on identically damped rows errs > 0.5 mHa at f=0.9 and > 1 mHa at f=0.7
  (measured 1.02 / 6.39 — damping is *not* harmless for GEVP methods).
- **G2 — damped + shot noise, matched data (DEFINITION OF DONE).** 10⁵ shots, 100 seeds: median
  ODMD error < 1 mHa at f=0.9 and < 5 mHa at f=0.7 (measured 0.43 / 2.64); KQD-to-ODMD median
  ratio > 100 at f=0.9 and > 50 at f=0.7 (measured ~2600× / ~150×).
- **G3 — the modulus window is the mechanism.** At f=0.9 under noise, the unit-modulus
  `odmd_energy` default is ≥ 2× worse than the wide-window device estimator (measured 3.75×);
  at f=1 both agree with the exact energy (< 1e-8 Ha).
- **G4 — Aer end-to-end phase survival and its boundary.** (a) Zero noise: median
  `|E − e_circuit| < 0.01 mHa` (measured 0.001 — the full Hadamard-test stack is faithful);
  (b) cx error 3e-4: amplitude loss ≥ 50% (`|s₇|/|s₇^exact| < 0.5`, measured 0.30) yet median
  phase error < 0.5 mHa (measured 0.052); (c) cx error 1e-3: damping to < 5% of the signal
  pushes it under the shot floor and the error exceeds 1 mHa (measured 5.2) — immunity ends at
  the noise edge, as the visibility law predicts.

## 6. Implementation plan (test-first)

1. `tests/test_device_odmd_spec.py` encoding G1–G4 (RED — `device_odmd` missing).
2. `device_odmd.py` + the 3-line `HardwareKrylovSolver.measure_signal`; everything else is
   composition of pinned primitives (`odmd_spectrum`, `noise_edge`, `build_trotter_odmd_problem`,
   `build_depolarizing_noise_model`).
3. `make gates`; re-run the non-spec suite (hardware path touched).

## 7. Out of scope

- Real-hardware runs; coherent/biased noise models; ZNE on the signal (machinery exists —
  a follow-up can gate ODMD+ZNE).
- Richardson Trotter-bias removal under device noise (compose with `SPEC_trotter_odmd` later).
- Excited states under device noise (compose with `SPEC_odmd_excited` later).

## 8. Caveats and risks

- **R1 — the immunity is exact only for *global* depolarizing:** local gate noise induces a real
  (small) phase bias — G4(b) measures it rather than assuming zero; a different noise model can
  enlarge it.
- The Aer arm inherits the Trotter bias by construction; always quote device-ODMD energies
  against the circuit eigenphase or after Richardson removal.
- Amplitude floors trade robustness for reach: `amp_floor=0.05` hides states with
  `p_n < 0.05` (consistent with the visibility law; dark states stay dark).

## 9. Deliverables

- `device_odmd.py`; `hardware_krylov.py` + `measure_signal`.
- `tests/test_device_odmd_spec.py` — gates G1–G4.
- `BACKLOG.md` entry with the measured immunity/boundary numbers.
