# SPEC: Circuit-real ODMD — Trotter eigenphases, the δt² law, and Richardson bias removal

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_trotter_odmd_spec.py`).

---

## 1. Goal

Make the ODMD rung circuit-real, with the Trotter bias *measured and then removed* rather than
assumed away. Three claims: (i) a fixed-step Suzuki-Trotter circuit is **exact** evolution of an
effective Hamiltonian, so ODMD on the circuit signal returns the ground eigenphase of the circuit
unitary U_trot to machine precision — DMD stacks **no additional approximation** on top of
Trotterization; (ii) that eigenphase's bias vs FCI follows the second-order law ∝ δt²_eff
(ratio ≈ 4 per reps doubling); (iii) two-point Richardson extrapolation across reps removes it to
< 0.1 mHa. Each is numerically falsifiable against exact references (dense diagonalization of
U_trot; FCI).

## 2. Background, the found defect, and honest framing

- **Found while probing this spec (the real headline):** `Operator()` and `Statevector.evolve()`
  evaluate an *opaque* `PauliEvolutionGate` through its **exact matrix**, silently ignoring the
  attached `SuzukiTrotter` synthesis (deviation from exact: 3.5e-16; the *decomposed* circuit:
  1.8e-2). Consequently `TrotterKrylovSolver` — which evolves by `Statevector.evolve(step)` — has
  been performing **exact evolution all along**: its `trotter_order`/`trotter_reps` knobs did
  nothing, and `test_trotter_circuit.py`'s "reaches FCI within Trotter error" passed vacuously.
  This is precisely the "validation that validates nothing" failure mode this repo was rebuilt to
  eliminate (cf. the qDRIFT/QCIVET post-mortems in `REFACTOR_PLAN.md`). **Fix:** materialize the
  synthesized definition in `build_trotter_step`, so every consumer (Statevector, Operator, Aer,
  the ancilla-controlled hardware path) sees the same genuinely Trotterized unitary. G1 is the
  regression gate — it fails on the pre-fix code by ~15 orders of magnitude.
- Prior art: Trotter effective-Hamiltonian / eigenphase-error theory (Childs et al.,
  `arXiv:1912.08854`); polynomial extrapolation of Trotterized phase estimation is established
  practice (e.g. Rendon, Watkins & Wiebe, `arXiv:2212.14144`). **Reproduction of known theory**;
  the composition with ODMD and the falsifiable packaging are what is new here.
- **What we can claim if gates pass:** the ODMD pipeline runs on true circuit dynamics; the only
  systematic error is the effective-Hamiltonian shift, whose size is measured, whose scaling law
  is verified, and which extrapolation removes below 0.1 mHa — with the noisy-budget win
  quantified.
- **What we cannot claim:** statevector-simulated circuits (no device noise/transpilation to a
  native gate set); Richardson doubles the deepest circuit (reps 2→4 costs 2× depth) — depth, not
  shots, is the price; minimal-basis systems only; the Krylov *subspace* path partially
  self-corrects Trotter error (H is measured exactly there), so the raw-eigenphase bias gated
  here is the *worst case*, not what `TrotterKrylovSolver` reports.

## 3. Approach

Same centered frame as `SPEC_odmd.md` (shift by μ **before** synthesis — the identity term only
adds a global phase and does not change the splitting error). Signal
`s_k = <phi0| U_trot^k |phi0>` by repeated circuit application (`Statevector.evolve` on the
*materialized* step); ODMD via the pinned `odmd.odmd_energy`. References: ground eigenphase of
`U_trot` among HF-reachable eigenvectors (dense `eig` of the circuit's `Operator` — exact);
FCI for the bias; Richardson `(4 E_fine − E_coarse)/3` for reps ratio 2, order 2.

## 4. Public interface

```
trotter_odmd.TrotterODMDProblem            # ODMDProblem subclass (odmd.sample_odmd_energy
                                           #   applies unchanged) + reps, order, e_circuit,
                                           #   unitary_deviation, depth
trotter_odmd.build_trotter_odmd_problem(mh, n=24, reps=1, order=2) -> TrotterODMDProblem
trotter_odmd.richardson_energy(e_coarse, e_fine, step_ratio=2.0, order=2) -> float
```

Changed (the fix): `trotter_krylov.build_trotter_step` now returns the **materialized**
Suzuki-Trotter circuit (same signature/semantics on hardware paths, which always decomposed).

## 5. Acceptance criteria (validation gates)

`tests/test_trotter_odmd_spec.py`; systems H2 and H4 chain; K=24; reps ∈ {1, 2, 4}, order 2.

- **G1 — the step is genuinely Trotterized (regression gate for the fix).**
  `||U_trot − U_exact||₂ > 0.05` (H2) / `> 0.02` (H4) at reps=1, and the deviation shrinks with
  ratio ∈ [3.5, 5.5] per reps doubling (measured 4.97/4.20 and 4.32/4.07). Pre-fix code measures
  ~1e-16 and fails loudly.
- **G2 — DMD adds nothing.** `|E_odmd − ground eigenphase of U_trot| < 1e-9 Ha` at every
  system × reps (measured ≤ 5e-11).
- **G3 — the δt² eigenphase law.** Bias(reps) = E_U − FCI has ratio ∈ [3.3, 5.5] per reps
  doubling on both systems (measured H2: 4.65, 4.14; H4: 4.12, 4.03; biases 19.3 → 1.0 mHa H2,
  12.0 → 0.72 mHa H4).
- **G4 — Richardson removes the bias (DEFINITION OF DONE).** Noiseless reps-(2,4) extrapolation:
  residual < 0.1 mHa on both systems and ≥ 10× below the reps=4 bias (measured 0.046 mHa / 22×
  on H2, 0.0065 mHa / 112× on H4). Under shot noise (H4, 10⁶ shots/element, 60 seeds,
  independent draws per reps): median |Richardson − FCI| < ⅓ × median |plain reps=2 − FCI|
  (measured ~0.38 vs ~3.3 mHa, ~8.6×).

## 6. Implementation plan (test-first)

1. `tests/test_trotter_odmd_spec.py` encoding G1–G4 (RED: module missing; G1 additionally RED
   against the unfixed `build_trotter_step`).
2. Fix `build_trotter_step` (materialize `gate.definition`); update the `trotter_krylov`
   docstrings that described the vacuous behavior as real.
3. `trotter_odmd.py` — minimum code; the noisy path reuses `odmd.sample_odmd_energy` verbatim
   via subclassing.
4. Re-run `tests/test_trotter_circuit.py` (its FCI-convergence tolerance must now survive REAL
   Trotter error) and the pinned ODMD gates; `make gates`.

## 7. Out of scope

- Device noise / transpiled-depth studies on the Trotterized signal (the machinery is now real;
  a follow-up can gate it).
- Higher-order Suzuki formulas, adaptive step selection, Chebyshev/multi-point extrapolation.
- Excited-state Trotter-ODMD (compose with `SPEC_odmd_excited.md` later).

## 8. Caveats and risks

- **R1 — extrapolation trusts the asymptotic law:** at large δt (H2 reps 1+2) higher-order terms
  leave 0.9 mHa residual — extrapolate from the *fine* pair, and only when bias > noise.
- The 4× depth of reps=4 circuits is the real cost on hardware; shots are not.
- The fix changes `TrotterKrylovSolver`'s numerical results (they now contain genuine Trotter
  error — previously they silently matched the exact solver).

## 9. Deliverables

- `trotter_krylov.py` — the materialization fix + honest docstrings.
- `trotter_odmd.py` — problem builder, `richardson_energy`, `__main__` demo.
- `tests/test_trotter_odmd_spec.py` — gates G1–G4.
- `BACKLOG.md` entry recording the found defect and the measured numbers.
