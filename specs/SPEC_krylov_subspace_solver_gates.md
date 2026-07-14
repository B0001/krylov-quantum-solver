# SPEC: krylov_subspace_solver — the documented bug fix gets a regression test, and the cross-check finally happens

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`krylov_subspace_solver.py` is a SECOND, independent real-time Krylov implementation — FCI-direct
(pyscf `contract_2e`), not qubit-mapped like `hybrid_quantum_solver.quantum_krylov_solver` — whose
own docstring says its purpose is to "cross-check SQD" and be "a second, independent ground-state
estimator." Neither of its two documented claims has ever been gated: (1) variational convergence to
CASCI (checked only by an informal `assert` in `__main__`), and (2) "the collapse fix" — a named,
documented historical bug (an absolute overlap-eigenvalue cutoff that "dropped almost every vector
and nulled the eigenproblem" on a poorly-scaled S) fixed by switching to a relative cutoff, with NO
regression test guarding it. This spec gates both, and finally runs the cross-check the module was
built for: does it agree with the OTHER validated Krylov implementation in this repo? False if the
relative cutoff isn't actually scale-invariant, if an independently-reimplemented absolute cutoff
doesn't reproduce the documented collapse, or if the two independent implementations disagree.

## 2. Background and honest framing

- `krylov_subspace_solver.py` already reuses validated primitives (pyscf `fci.direct_spin1` for the
  dense active-space Hamiltonian) — no new physics, only falsifiers around claims the module already
  makes but never checks beyond a single-seed `__main__` printout.
- **What you can claim if the gates pass:** the module stays variational (never below CASCI) across
  closed- and open-shell systems (H2, H4, O2 triplet); the documented "collapse fix" is proven to
  matter, not just asserted — the relative cutoff keeps an IDENTICAL number of basis vectors under a
  benign rescaling of the reference state (7/12, every tested scale), while an independently
  reimplemented absolute cutoff at the same threshold collapses to 0/12 kept vectors at scale <=
  1e-6, reproducing the exact failure the docstring describes; and the module agrees with
  `hybrid_quantum_solver.QuantumKrylovSolver` (an entirely independent, qubit-mapped statevector
  implementation) to sub-microhartree precision on the same active space — the cross-check this
  module was built to provide, run for the first time.
- **What you cannot claim:** that the specific condition numbers/kept-vector counts generalize
  beyond the tested systems; that this validates the module for active spaces beyond the "a few
  thousand determinants" the module's own docstring already scopes `_dense_active_H` to.
- **Reference:** CASCI (dense diagonalization, the module's own comparison target) for G1/G4; the
  independently-reimplemented absolute-cutoff logic (constructed in the test, not touching library
  code) for G2; `hybrid_quantum_solver.QuantumKrylovSolver`'s converged energy for G3.

## 3. Approach

Reuse `krylov_ground_state`, `krylov_convergence_sweep`, `_dense_active_H` unmodified. G1: run the
convergence sweep on three systems (H2 CAS(2,2), H4 CAS(4,4), O2 CAS(4,4) triplet — the same three
`__main__` already uses) and check variational + convergence. G2: independently reconstruct the
overlap matrix S (mirroring `krylov_ground_state`'s exact internal formula) at a fixed
`dt`/`krylov_dim`, sweep an artificial rescaling of the reference state `phi0` (a benign
normalization change that must not alter the physics), and compare kept-vector counts under the
RELATIVE threshold (the actual fix) vs. an independently-reimplemented ABSOLUTE threshold at the
same numeric cutoff. G3: build the same active-space Hamiltonian through BOTH this module's
FCI-direct path and `hybrid_quantum_solver.build_molecular_hamiltonian` + `QuantumKrylovSolver`'s
qubit-mapped path, compare converged energies. G4: sweep `dt` and confirm the docstring's "condition
numbers of 1e6+ are normal and harmless" claim against measured condition numbers and errors.

## 4. Public interface

No new library code — this spec adds only test-file assertions and one independent reconstruction
helper (the "old" absolute-cutoff logic, deliberately NOT touching `krylov_subspace_solver.py`)
around its existing public functions, reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — variational floor + convergence, pinning the existing informal assertion.** On H2
  CAS(2,2), H4 CAS(4,4), O2 CAS(4,4) triplet (open-shell): every step of `krylov_convergence_sweep`
  is variational (`energy >= CASCI - 1e-6`), and the deepest tested dimension (m=12) is within 5 mHa
  of CASCI. *Measured final deltas: H2 0.0000, H4 0.0047, O2-triplet 0.0000 mHa.*
- **G2 — THE FINDING (definition of done): the collapse fix is proven, not asserted.** Under an
  artificial rescaling of the reference state by `{1, 1e-3, 1e-6, 1e-9}` (dt=0.5, krylov_dim=12,
  H4 CAS(4,4)): the RELATIVE cutoff keeps the SAME number of vectors at every scale (scale-invariant,
  as required for a well-posed generalized eigenproblem); an independently-reimplemented ABSOLUTE
  cutoff at the identical numeric threshold collapses to 0 kept vectors at scale `<= 1e-6`.
  *Measured: relative kept=7/12 at every scale; absolute kept 8, 4, 0, 0 as scale drops
  1, 1e-3, 1e-6, 1e-9 — reproducing the docstring's "nulled the eigenproblem" exactly.*
- **G3 — the cross-check the module was built for, run for the first time.** On H4 CAS(4,4), this
  module's converged energy (FCI-direct, krylov_dim=16) agrees with
  `hybrid_quantum_solver.QuantumKrylovSolver`'s converged energy (qubit-mapped statevector,
  independent code path) to `< 1e-3` mHa. *Measured: 0.00037 mHa difference (both essentially exact
  vs CASCI).*
- **G4 — the docstring's condition-number claim, measured not asserted.** Sweeping `dt` on H4
  CAS(4,4): condition numbers reach `>= 1e6` at multiple tested `dt` values while the energy stays
  within 3 mHa of CASCI at every one — "harmless" is a checked claim, not a description.
  *Measured: cond 2.6e6 (dt=0.1, err 0.45 mHa) through 9.5e6 (dt=2.0, err ~0 mHa), all `< 3` mHa.*

> Definition of done: **G2**. A documented historical bug with no regression test is exactly the gap
> this repo's culture exists to close — G1/G3/G4 build the surrounding confidence, G2 is the fix.

## 6. Implementation plan (test-first)

1. Write `tests/test_krylov_subspace_solver_gates_spec.py` encoding G1-G4 (RED in the sense these
   checks are new, even though `krylov_subspace_solver.py`'s functions are not).
2. No changes to `krylov_subspace_solver.py` — G2's absolute-cutoff comparison is deliberately
   reimplemented in the test file, not imported from the module (it is the FORMER, buggy behavior).
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- Fixing/changing `krylov_subspace_solver.py` itself — it already has the fix; this spec only gates
  it.
- Active spaces beyond the module's own documented `_dense_active_H` scope (a few thousand
  determinants).
- A device-noise or shot-sampled version of the cross-check (both compared implementations are
  exact-statevector here).

## 8. Caveats and risks

- **R1 — G2's rescaling construction is a deliberately artificial stress test**, not something
  `krylov_ground_state` would naturally produce internally (its `phi0` is always exactly a unit
  computational-basis vector). It demonstrates the MECHANISM the fix guards against on this system's
  actual S matrix, not that the natural parameter ranges tested in G1/G4 would themselves collapse
  without the fix (they were measured to differ by at most 1 kept vector at G4's tested dt values —
  see probe data in the PR).
- Honest limitation: three systems, one active-space size regime (CAS(4,4) and smaller).

## 9. Deliverables

- `tests/test_krylov_subspace_solver_gates_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
