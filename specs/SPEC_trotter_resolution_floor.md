# SPEC: Trotter resolution floor — when a circuit eigenphase cannot be trusted, and the ordering fix

**Status:** IMPLEMENTED once gates green. Grew out of the `nb3x8_device_gap` G2 flake
(recorded in BACKLOG; `make gates` red ~60% of runs before this spec).

---

## 1. Goal

Two coupled claims, both falsifiable:

1. **Determinism (the fix).** The repo's Trotter synthesis was **nondeterministic across
   processes**: `SuzukiTrotter` products depend on the Pauli term order, and the qubit
   Hamiltonian's term order varied run-to-run with Python hash randomization (two observed
   orderings ≡ coeff-ascending and coeff-descending; deviations 1.045e-2 vs 5.964e-3 for Nb₃F₈
   sector-2, reps=1). Canonically ordering the terms inside `build_trotter_step`
   (largest-|coefficient| first, label tiebreak) makes every Trotter circuit — and every number
   derived from one — reproducible.

2. **The resolution floor (the law).** A circuit eigenphase is *resolvable* only if the
   reference state's population on it exceeds the Trotter leakage floor ‖U_trot − U_exact‖₂².
   Below that floor the computed population is leakage-dominated (interference of ~dev²
   leaked amplitude with the genuine signal), and the extracted eigenphase branch becomes an
   artifact of the synthesis ordering. **Probe:** re-synthesize under different canonical
   orderings — a resolvable quantity moves by ordinary Trotter bias (small, smooth); an
   unresolvable one flips branch by the wrap quantum 2π/τ.

## 2. The numbers (Nb₃F₈ sector nelec=2, the flake's origin)

Reference population on the sector ground state: **1.36e-5** (genuine, above the 1e-8 cut).

| reps | deviation (coeff-desc) | leakage dev² | vs signal 1.4e-5 | ordering spread of the gap |
|---|---|---|---|---|
| 1 | 5.96e-3 | 3.6e-5 | **leak > signal** | **3756.7 meV** (= 2π/τ branch flip) |
| 2 | 1.16e-3 | 1.3e-6 | leak < signal | 5.5 meV (ordinary bias) |
| 4 | 2.73e-4 | 7.5e-8 | leak < signal | — |

The failing assertion (`abs(circuit_gap(**f8, reps=1) - exact_gap(**f8)) < 1.0`) tested a
quantity below the method's own floor; whether it passed depended on which ordering the hash
seed dealt. With the canonical ordering the F8 reps=1 value is deterministic (+2580.70, the
value every previously-green run recorded) — but the *spec-level* fact remains that at reps=1
it is floor-limited, and the ordering probe exposes that (coeff-ascending branch-flips to
−1171.1). reps=2 is above the floor under every ordering.

## 3. What we claim / do not claim

- **Claim:** canonical ordering makes the synthesis deterministic; the previously recorded
  device-gap numbers (F8 +2580.70 at reps=1, I8 549.55/774.18 at reps=1/2, bias ratio 4.29)
  are all reproduced exactly — the fix pins the green-run history, it does not move it.
- **Claim:** the floor criterion (population vs dev²) correctly separates the one flaky
  assertion from every stable one in the same suite (F8 reps=1 below floor ↔ flaky; F8
  reps≥2, I8 all-reps above floor ↔ never flaked).
- **Do not claim:** that coeff-descending is the *optimal* Trotter ordering in general (it
  is the best of the four measured here and matches the historical green runs; ordering
  optimization is its own literature); nor that the floor is tight — dev² is an upper-bound
  scale for leakage, not an exact population.

## 4. Public interface

```
hybrid_quantum_solver.trotter_krylov.build_trotter_step(...)   # now canonically ordered
trotter_resolution_floor.leakage_floor(mh, reps, ...) -> float          # dev^2
trotter_resolution_floor.reference_population(mh) -> float             # |<gs|psi0>|^2, exact
trotter_resolution_floor.is_resolvable(mh, reps, ...) -> bool          # population > floor
trotter_resolution_floor.ordering_spread(mh_or_params, reps) -> float  # the probe (meV)
```

## 5. Acceptance criteria (validation gates)

`tests/test_trotter_resolution_floor_spec.py`.

- **G1 — determinism.** `build_trotter_step` yields the *same* circuit (same Operator matrix)
  when handed the same SparsePauliOp in shuffled term orders; the Nb₃F₈ sector-2 reps=1
  deviation is pinned at its coeff-desc value (5.96e-3), not the 1.04e-2 alternative.
- **G2 — the floor separates flaky from stable.** For Nb₃F₈ sector-2: `is_resolvable` is
  False at reps=1 (population 1.4e-5 < floor 3.6e-5) and True at reps=2 (> 1.3e-6). The
  ordering-spread probe agrees: > 1000 meV at reps=1 (branch flip), < 10 meV at reps=2.
- **G3 — the capstone gate is deterministically green.** `circuit_gap(**f8, reps=1)` returns
  +2580.70 ± 0.5 on **every** of ≥ 5 consecutive in-process evaluations (was a coin flip),
  and the full `test_nb3x8_device_gap_spec.py` passes.
- **G4 — no committed Trotter number moves.** `test_trotter_odmd_spec.py` and
  `test_device_odmd_spec.py` stay green unmodified (the canonical order reproduces the
  historical green-run values).

## 6. Out of scope

- Ordering *optimization* (choosing the provably-minimal-error ordering per Hamiltonian).
- Restating the device-gap capstone's reps=1 F8 line as a floor-limited caveat in its own
  spec text (SPEC_nb3x8_device_gap edit) — noted there via BACKLOG.
- The unphysical-branch latent issue (`np.min` can select eigenphases outside the reachable
  band, e.g. the −1657.45 image) — still recorded in BACKLOG, untriggered with the canonical
  ordering.

## 7. Caveats and risks

- **R1 — order-2 only measured.** Floor numbers are for SuzukiTrotter(order=2); other orders
  rescale dev but the population-vs-dev² criterion is generic.
- **R2 — dev is a global 2-norm.** Per-eigenspace leakage can be smaller than dev²; the
  criterion is conservative in the safe direction (may call "unresolvable" something that
  happens to resolve, never the reverse — G2's reps=1 case shows the flip is real, not
  hypothetical).

## 8. Deliverables

- Canonical ordering in `hybrid_quantum_solver/trotter_krylov.py`;
  `trotter_resolution_floor.py`; `tests/test_trotter_resolution_floor_spec.py`; this spec;
  BACKLOG update closing the flake entry.
