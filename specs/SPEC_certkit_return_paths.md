# SPEC: Every ground-state energy the solver returns can carry a checked certificate

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_certkit_return_paths_spec.py`).
Bead `portfolio-cir`. Making it non-bypassable is `portfolio-0pc`, not this spec.

---

## 1. Goal

Claim: every energy-returning method of the Krylov solver can emit a certkit v1
`lambda_min_enclosure` certificate for the ground-state energy it returns, checked out of
process by the pinned checker, and **the certified ceiling is the returned number itself**. So
VERIFIED means the checker has proven that energy an upper bound on the ground state, and an
energy below the ground state cannot verify. Falsifiable: one returned energy without a
certificate, one VERIFIED certificate whose ceiling is not the returned energy or whose
enclosure misses the exact `lambda_min`, or one sub-ground-state energy that verifies, fails
the gate.

## 2. Background and honest framing

- `SPEC_certkit_regression_gate.md` certified four fixed cases from `certkit_bridge.run_case`
  and listed this bead as its coverage limit: no solver return path produced a certificate.
- **The finding that shaped the design.** `run_case` claims `[floor, <x|H|x>]`, with the
  ceiling set by the witness. That is a true enclosure for any unit vector, so it says nothing
  about the energy the solver *returned*. With shot noise the two differ: in G3's sample a
  noisy H2 energy 7.3 mHa **below** the exact ground state got a VERIFIED certificate under
  that claim (checked by mutating the producer). Anchoring the ceiling to the returned
  energy turns the same run into an ABSTAIN.
- **What we can claim if the gates pass:** on every return path, with `certify=<dir>`, the
  returned ground-state energy comes with a certificate that the pinned checker re-derives
  from the files alone, out of process. A VERIFIED one proves `lambda_min <= E` (electronic
  frame) and gives a certified floor. A noisy energy below `lambda_min` gets ABSTAIN.
- **What we cannot claim:**
  1. **Opt-in, so bypassable.** Without `certify=`, energies come back uncertified as before.
     Making certificate-or-abstention unavoidable is `portfolio-0pc`.
  2. **Global, not sector, minimum.** The claim is about the whole qubit matrix. Where the
     global minimum lies in another particle sector, a VERIFIED certificate is true but
     loose: the Hubbard dimer's returned energy is 0.17 Ha above the global `lambda_min`.
  3. **The floor is only as good as the route.** `temple_inertia` needs certkit to
     materialise the matrix (`n <= DENSE_LIMIT` = 256, i.e. up to 8 qubits) and a finite
     Temple bound. Otherwise `gershgorin_rayleigh` certifies, and its floor is loose
     (2.7 Ha below on N2 CAS(6,6)). The ceiling is tight on every route.
  4. **Ground state only.** `solve_excited` and `eigenstates` certify `energies[0]`; certkit
     v1 has no claim kind for excited levels.
  5. **The hardware witness is rebuilt classically.** `HardwareKrylovSolver` measures only
     the subspace matrices; its certificate carries `sum_k c_k U^k|HF>`, rebuilt as a
     statevector. That works because the device is simulated. On real hardware it would not
     scale.
  6. Same checker caveat as the sibling gate: out of process from the pinned release, but in
     the producer's own venv, not a dedicated one.

## 3. Approach

`certkit_bridge.certify_ground_state(mh, psi, energy, out)`:

1. Witness `x`: the lower-Rayleigh-quotient one of `Re(psi)` and `Im(psi)`, normalised. `H` is
   real symmetric, so `<psi|H|psi>` is their norm-weighted mean and `x` is never above
   `psi`'s own quotient.
2. Claim `[floor - pad, energy + pad]`, with `pad_claim` padding as in `run_case`.
3. Routes, tightest first: `temple_inertia` (full-space `beta` from `full_space_beta`,
   discharged by the checker's inertia count) when `n <= DENSE_LIMIT` and `<x|H|x> < beta`;
   then `gershgorin_rayleigh`. The first VERIFIED wins; otherwise the last abstention is
   returned. An energy below the floor makes an inverted enclosure, which the checker refuses.
4. Each certificate is checked by `certkit check` in a subprocess, as `run_case` does.

**Reference:** exact `lambda_min` by dense diagonalisation, `mh.ground_state_energy() -
mh.energy_offset` (cross-checked against PySCF FCI in `test_reference_energies.py`), and the
N2 value pinned and re-derived in the sibling gate.

## 4. Public interface

```
certkit_bridge.certify_ground_state(mh, psi, energy, out) -> Verdict
certkit_bridge.real_witness(H, psi)                       -> real unit vector
hybrid_quantum_solver.quantum_krylov_solver.certify_state(mh, psi, energy, out) -> Verdict

QuantumKrylovSolver.solve(m, certify=None)               -> KrylovStep(.certificate)
QuantumKrylovSolver.convergence(m, certify=None)         -> [KrylovStep], files in <dir>/M<m>/
QuantumKrylovSolver.solve_excited(m, n, certify=None)    -> ExcitedKrylovStep(.certificate)
QuantumKrylovSolver.eigenstates(m, n, certify=None)      -> (energies, states[, certificate])
TrotterKrylovSolver.solve / .convergence(..., certify=None)
HardwareKrylovSolver.solve / .convergence(..., certify=None)
SampleKrylovSolver.solve / .convergence(..., certify=None)     files in <dir>/shots<s>/
run_geometry / run_from_integrals(..., certify=None)     -> PipelineResult.certificates
```

`certify` is a directory; it holds `operator.json` plus one certificate file per route tried.
`certificate` is `None` when it is not given. Energies are unchanged either way (G1).

## 5. Acceptance criteria (validation gates)

- **G1 — every path, every energy, VERIFIED, by the pinned route, energy unchanged.** All
  13 paths on H2 (plus `run_from_integrals` on the Hubbard dimer) return a certificate for
  every energy. The checker's exit code and its own output line agree. The route is
  `temple_inertia` on H2 and `gershgorin_rayleigh` on the dimer. The energies match an
  uncertified run to 1e-12 Ha: not bitwise, because the default `dt` comes from ARPACK with
  a random start vector, so two solver instances already differ by ~1e-15.
- **G2 — DEFINITION OF DONE.** For every certificate, `lo <= lambda_min <= hi`, and
  `0 <= hi - E <= 1e-7` Ha: the certified ceiling is the returned energy plus float padding.
- **G3 — it can fail.** 12 seeded noisy H2 solves (`noise_sigma=3e-2`, M=6). Every energy
  below `lambda_min` ABSTAINs, by verdict and not by crash. Every VERIFIED one contains
  `lambda_min`. At least one seed dips below, so the gate is not vacuous; today six do.
- **G4 — beyond the dense limit.** N2 CAS(6,6) (4096 states): VERIFIED via
  `gershgorin_rayleigh`, exact energy inside, ceiling equal to the returned energy.

Checked by mutation: claiming the witness's quotient as the ceiling (the `run_case`
convention) turns G2 and G3 red.

## 6. Out of scope

- Certificate-or-abstention as a non-bypassable invariant (`portfolio-0pc`), and the CI gate
  over every regression test's certificate (`portfolio-zwd`).
- Top-level study scripts (`krylov_subspace_solver.py`, `qdrift_krylov.py`, ODMD and so on).
  They are harnesses, not the solver, and ODMD energies have no state to serve as a witness.
- A shared inertia certificate per operator (`temple_ref`). Each energy is checked on its
  own, so a certified `convergence` on 8 qubits pays one O(n³) interval factorisation per
  step.

## 7. Caveats and risks

- **Cost.** One checker process per certified energy, off by default. That is about 0.1 s on
  H2, seconds on 8 qubits (Temple), and about 10 s on N2 (Gershgorin, matrix-free).
- `emit` prints the checker's line for every certificate, as `run_case` does, so certified
  library calls print one line per energy.
