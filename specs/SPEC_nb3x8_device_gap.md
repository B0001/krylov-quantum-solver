# SPEC: Nb₃X₈ cluster charge gaps through the simulated-hardware pipeline (capstone study)

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_nb3x8_device_gap_spec.py`).

---

## 1. Goal

*A scientific study / capstone composition, not a new method rung* (precedent:
`SPEC_nb3x8_gaps.md`). Claim: the repo's two validated threads — the Nb₃X₈ downfolded dimer
clusters (`nb3x8_gaps.py`, exact ED gaps from `arXiv:2501.10320`'s cRPA parameters) and the
device-validated ODMD stack (`SPEC_odmd` → `SPEC_trotter_odmd` → `SPEC_device_odmd`) — compose
into an end-to-end pipeline that measures a **material cluster's charge gap
Δ = E(N+1) + E(N−1) − 2E(N) through simulated quantum hardware**: three ground-state ODMD runs
(one per particle sector, each depolarizing-immune), genuinely-Trotterized Hadamard-test
circuits, an Aer device noise model, and Richardson removal of the Trotter bias. Falsifiable
against exact sector FCI at every stage, including the recorded 842.44 meV Nb₃I₈ LT-bulk gap.

## 2. Background and honest framing

- **What we can claim if gates pass:** the Nb₃I₈ dimer gap from a noisy simulated device to
  ~10 meV (1.2%), with every error source separately measured — Trotter bias (large, and it does
  **not** cancel between particle sectors), device-noise floor, and their crossover governing
  when extrapolation pays.
- **What we cannot claim:** this is the **isolated-cluster** gap — per `SPEC_nb3x8_gaps.md`'s
  corrected conclusion it is an *upper bound* on the solid gap (band broadening moves the real
  material to ~600–650 meV, vindicating the paper's Hubbard-I); density-density interactions
  only; a simulated device (no coherent errors/crosstalk/drift); charged-sector references are
  trivial for a dimer (the N=1/N=3 sectors are nearly free) — this is a pipeline demonstration
  at validation scale, not a hard correlated-electron benchmark.

## 3. Approach

Per material (`NB3X8_LT_BULK` cRPA parameters, meV units end to end): build the three sector
models (N = 1, 2, 3 electrons on the generalized Hubbard dimer; real-time evolution conserves
particle number, so each HF reference pins its sector); run ODMD per sector at three fidelity
levels — exact statevector signal (`odmd.build_odmd_problem`), exact eigenphases of the
Trotterized step circuit (`trotter_odmd.build_trotter_odmd_problem`), and Aer-measured
Hadamard-test signals (`device_odmd.measure_survival_signal`); form the gap; Richardson-remove
the Trotter bias across `reps` (`trotter_odmd.richardson_energy` — linear, so it commutes with
the gap combination). References: `fixed_filling_energy` (sector FCI, pinned by
`SPEC_nb3x8_hubbard.md`) and the recorded exact gaps of `SPEC_nb3x8_gaps.md`.

## 4. Public interface

```
nb3x8_device_gap.sector_models(U0, t, Us) -> dict[int, ModelIntegrals]   # N = 1, 2, 3
nb3x8_device_gap.exact_gap(U0, t, Us) -> float                            # sector-FCI reference
nb3x8_device_gap.statevector_gap(U0, t, Us, n=16) -> float
nb3x8_device_gap.circuit_gap(U0, t, Us, reps, n=16) -> float              # exact Trotter eigenphases
nb3x8_device_gap.device_gap(U0, t, Us, shots, noise_model, seed,
                            trotter_reps=1, n=8) -> float
nb3x8_device_gap.device_gap_richardson(U0, t, Us, shots, noise_model,
                                       seed, n=8) -> float                # reps (1,2) pair
```

Pure composition — every primitive is already spec-pinned.

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_device_gap_spec.py`. All energies in meV. Aer runs: 32768 shots/observable,
5 seeds, medians (Aer is seeded → deterministic gates).

- **G1 — statevector pipeline is exact.** `|statevector_gap − exact_gap| < 0.01 meV` for all
  four LT-bulk materials (measured ≤ 4e-10), and the Nb₃I₈ `exact_gap` reproduces the recorded
  842.44 meV (< 0.5 meV) — the cross-spec pin to `SPEC_nb3x8_gaps.md`.
- **G2 — sector Trotter biases do NOT cancel in the gap (the recorded finding).** Nb₃I₈:
  `|circuit_gap(reps=1) − exact| > 50 meV` (measured −100.9, i.e. 12% of the gap!) with the
  order-2 ratio `bias(1)/bias(2)` ∈ [3.3, 5.5] (measured 4.65); contrast Nb₃F₈:
  `|bias(reps=1)| < 1 meV` (measured 0.10 — tiny hopping → near-commuting Hamiltonian).
- **G3 — Richardson fixes the circuit-exact gap.** Reps-(2,4): residual < 1 meV (measured 0.27)
  and > 5× below the reps=4 bias (measured 19×); reps-(1,2): residual < 10 meV (measured 4.7).
- **G4 — the device measurement and the bias-vs-noise crossover (DEFINITION OF DONE).** At
  cx = 1e-4: median `|device_gap_richardson − exact| < 15 meV` (measured 10.4 — 1.2% of the
  gap) and > 5× below the raw reps=1 device gap error (measured 99.7 → 9.6×). At cx = 3e-4 the
  noise floor exceeds the reps=2 bias and Richardson stops paying
  (`median rich ≥ median plain(reps=2)`, measured 16.3 vs 9.6) — `SPEC_trotter_odmd` R1
  ("extrapolate only when bias > noise") demonstrated on a material.

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_device_gap_spec.py` encoding G1–G4 (RED — module missing).
2. `nb3x8_device_gap.py` — composition only; no new numerics.
3. `make gates`.

## 7. Out of scope

- The solid-state (broadened) gap — the coordination/TDL machinery of `SPEC_nb3x8_gaps.md` is
  the classical answer; a circuit path to it would need larger clusters.
- Larger Nb₃X₈ clusters (4+ sites) on the device path; excited states per sector.
- Real-hardware execution; non-density-density interaction terms (few-meV, per the paper).

## 8. Caveats and risks

- **R1 — quote the right number:** the device pipeline reproduces the *isolated-cluster* 842 meV
  gap; the material's gap is ~600–650 meV after broadening. Never conflate them.
- The gap combination (weights 1, −2, 1) amplifies uncorrelated sector noise by √6 ≈ 2.4× —
  budget accordingly.
- Richardson on device data inherits noise amplification (~1.9× variance); G4's crossover marks
  where that trade flips.

## 9. Deliverables

- `nb3x8_device_gap.py` (+ `__main__` table over the LT-bulk family).
- `tests/test_nb3x8_device_gap_spec.py` — gates G1–G4.
- `BACKLOG.md` entry with the measured ladder (statevector → circuit → device).
