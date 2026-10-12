# SPEC: fix the SCF-convergence residue at its source — opt-in symmetric/tight SCF in the builder

**Status:** IMPLEMENTED. Gates in `tests/test_symmetric_scf_spec.py`. Results: §9.

---

## 1. Goal

`SPEC_eigenstate_reachability` §0 found that `QuantumKrylovSolver` started from `|HF⟩` genuinely
converges to the SCF-convergence-residue level on symmetry-forbidden geometries (square H₄ and
relatives) — so vetoing that level downstream, in the reference, was wrong; the residue is real
data for the Hamiltonian the loose SCF actually built. That spec's fix-forward closed with: "Where
the residue really belongs — a tight or symmetric SCF in the builder — is a follow-up." This is
that follow-up. **Claim:** adding an opt-in `symmetric_scf=True` flag to
`build_molecular_hamiltonian` that floors `conv_tol` to the already-validated
`reachability.TIGHT_SCF_CONV_TOL` (1e-13) whenever PySCF detects nontrivial point-group symmetry
removes the residue at its source (so `QuantumKrylovSolver` from `|HF⟩` lands on the true physical
ground state, not the artifact) while being a provable no-op — bit-identical build, to the limits
of independent-SCF-call floating noise — on every existing call site and every currently-pinned
gate, because the flag defaults to `False`.

## 2. Background and honest framing

- `conv_tol` is pinned at its loose default (1e-9) by `SPEC_scf_conv_tol` (chem-ayr) specifically
  because flipping it everywhere would move `SPEC_eigenstate_reachability`'s witnesses and every
  other gate that currently records the residue as present. This spec does not revisit that
  verdict: the default stays 1e-9, unconditionally, for every caller that does not pass
  `symmetric_scf=True`.
- **What you can claim** if the gates pass: callers who know their geometry is point-group
  symmetric and want the physical (non-residue) level can opt in with one keyword, get a build
  that is bit-identical to today's manual `conv_tol=TIGHT_SCF_CONV_TOL` workaround, and get it for
  free (a no-op) on asymmetric geometries so it is safe to pass unconditionally in a caller that
  doesn't know in advance whether its geometry is symmetric.
- **What you cannot claim:** this does not make the *default* behavior safer, does not implement a
  true symmetry-constrained SCF (it is a tolerance floor on PySCF's ordinary RHF, which can still
  spontaneously symmetry-break — see R1), and does not resolve the pre-existing, already-tracked
  Linux run-to-run flakiness in `test_scf_conv_tol_spec.py` / `test_eigenstate_reachability_spec.py`
  (chem-dy4, chem-ffe; §8 below has a platform finding relevant to those beads, logged there, not
  fixed here).

## 3. Approach

`build_molecular_hamiltonian` already accepts `conv_tol`. Add a point-group check
(`pyscf.gto.M(..., symmetry=True).groupname`, geometry-only, no SCF) and, when `symmetric_scf=True`
and the detected group is not `"C1"`, take `effective_conv_tol = min(conv_tol,
TIGHT_SCF_CONV_TOL)` before constructing `PySCFDriver`. The reference for "the residue is gone" is
the repo's own existing diagnostic machinery: `reachability._dense_hf_projection` (HF population
spectrum), `reachability.reachable_eigenpairs`/`reachable_mask` (does the residue get flagged?),
and an independent PySCF symmetry-adapted FCI (`test_eigenstate_reachability_spec._fci`) for the
"did QKSD converge to the *physical* level" claim (G7 below).

## 4. Public interface

```
hybrid_quantum_solver.molecular_hamiltonian.build_molecular_hamiltonian(
    ..., conv_tol: float = 1e-9, symmetric_scf: bool = False) -> MolecularHamiltonian
hybrid_quantum_solver.molecular_hamiltonian._point_group(atom, basis, charge, spin) -> str
```

No other module changes. `MolecularHamiltonian.build_args["conv_tol"]` now records the *effective*
tolerance (what was actually passed to `PySCFDriver`), not the caller's nominal input — load-bearing
for `reachability._driver_orbitals`, which rebuilds from `build_args` and must reproduce the same
MOs.

## 5. Acceptance criteria (validation gates)

- **G1 — defaults pinned.** `conv_tol` default is still `1e-9`; `symmetric_scf` default is `False`
  (signature introspection, no build).
- **G2 — point-group detection matches PySCF.** `_point_group` agrees with `pyscf.gto.M(...,
  symmetry=True).groupname` on a battery spanning `Dooh`/`D2h`/`C1`, with no SCF run (cheap,
  deterministic).
- **G3 — no-op on C1.** On an asymmetric 3-D H₄ geometry, `symmetric_scf=True` leaves
  `build_args["conv_tol"]` *exactly* equal to the input, and `|ΔE0| = 0`, `|ΔE_HF| = 0` against
  `symmetric_scf=False`.
- **G4 — harmless on the clean set.** On every `SPEC_scf_conv_tol.CLEAN` geometry (H₂, H₄ chain,
  square H₄ at sides where the residue is already absent at default `conv_tol`), `symmetric_scf=True`
  agrees with `symmetric_scf=False` within that spec's own G1 thresholds: `|ΔE0| < 1e-9` Ha,
  `|ΔE_HF| < 1e-8` Ha.
- **G5 — equivalent to the manual tight-tolerance workaround.** On symmetric square H₄ across a
  side-length sweep (1.05–1.40 Å, including both `SPEC_scf_conv_tol`'s witnesses), a
  `symmetric_scf=True` build starting from a deliberately loose `conv_tol=1e-6` is bit-identical to
  an explicit `conv_tol=TIGHT_SCF_CONV_TOL` build (`effective_conv_tol` floors to the same literal
  1e-13 either way) — `ΔE0 = ΔE_HF = Δ`overlap` = 0`.
- **G6 — residue collapse at a Linux-native witness (DEFINITION OF DONE, correctness half).**
  Thread-pinned (`OMP_NUM_THREADS=1` etc., the `test_hamiltonian_determinism.py`/chem-8tl
  convention — multithreaded BLAS reduction order is a separate, already-documented source of
  run-to-run noise in this exact SCF), at square H₄ side 1.2 Å, default `conv_tol`: the lowest-level
  HF population is **repeatably** `5.08e-9` (`> REACHABLE_TOL_CERTIFIED`, i.e. admitted and flagged
  by `reachable_mask` with a `RuntimeWarning`) without the flag, and **repeatably** `< 1e-20` (no
  warning) with `symmetric_scf=True`, across 5 independent builds each.
- **G7 — QKSD lands on the physical level, not the artifact (DEFINITION OF DONE, the bead's actual
  claim).** At the same witness, `QuantumKrylovSolver(mh).convergence(...)` from `symmetric_scf=False`
  reproduces `SPEC_eigenstate_reachability` G9 (converges toward the symmetry-forbidden B1g FCI
  level, `> 0.1` Ha below the physical Ag level); from `symmetric_scf=True`, every `QuantumKrylovSolver`
  energy sampled matches the Ag FCI energy to `< 1e-8` Ha and `reachable_eigenpairs(mh)[0][0]`
  equals it too.

G6 and G7 together are "done": G6 shows the residue is gone from the diagnostic's point of view:
G7 shows that absence is not cosmetic — it is the difference between the solver returning the right
answer and the wrong one.

## 6. Implementation plan (test-first)

Code landed before this spec file in the same working tree (opt-in parameter + `_point_group`
helper in `hybrid_quantum_solver/molecular_hamiltonian.py`); gates in
`tests/test_symmetric_scf_spec.py` encode G1–G7 above and were run to green against it
(`uv run pytest tests/test_symmetric_scf_spec.py -v`, §9).

## 7. Out of scope

- Flipping the global `conv_tol` default, or making `symmetric_scf` default `True` — explicitly
  rejected; see §2 and chem-ayr's prior verdict.
- A true symmetry-constrained SCF (PySCF's `symmetry=True` RHF with irrep occupation control) —
  this uses an ordinary unconstrained RHF with a tighter stopping tolerance, which fixes the
  *residue* (a convergence artifact) but not spontaneous symmetry breaking (a different, genuine
  physical phenomenon — §8 R1).
- Fixing the pre-existing Linux flakiness in `test_scf_conv_tol_spec.py` /
  `test_eigenstate_reachability_spec.py` (chem-dy4, chem-ffe) — logged as a finding in §8, not
  fixed here; this spec's own gates sidestep it by thread-pinning (G6) or by using geometries/paths
  unaffected by it (G1–G5, G7).

## 8. Caveats and risks

- **R1 — spontaneous symmetry breaking is not the residue and this does not touch it.** RHF can
  land on a symmetry-broken solution at some geometries regardless of `conv_tol` (a real physical
  degeneracy-splitting effect, `SPEC_symmetry_reachability`). Tightening `conv_tol` does not change
  *which* RHF solution is found, only how converged it is; `symmetric_scf=True` on a geometry where
  RHF spontaneously breaks symmetry will faithfully converge the broken solution to 1e-13, not
  "fix" it into the symmetric one.
- **R2 — platform finding (not fixed here, logged for chem-dy4/chem-ffe).** Measured on this Linux
  container: the RHF spontaneous-symmetry-breaking pattern at square H₄ is essentially *inverted*
  from the macOS numbers recorded in `SPEC_scf_conv_tol`/`SPEC_eigenstate_reachability` — e.g. sides
  1.05, 1.10, 1.19, 1.35, 1.40 broke symmetry on this machine (via `reachability.scf_symmetry_status`),
  while 1.00 and 1.20 stayed symmetric, roughly the opposite split from the macOS witnesses those
  specs were measured on. That is why this spec's G6/G7 witness is side 1.2, not the macOS witnesses
  1.19/1.35 — those do not reliably reproduce the residue on this platform (confirmed separately:
  `conv_tol=1e-13` builds at 1.19/1.35 on this machine kept a physical O(0.4) population, i.e. they
  are not even landing near the same solution family as the macOS run). Recorded on `chem-ffe`
  (notes) rather than a new bead, since that bead's stated purpose is exactly a Linux verification
  pass over these witnesses.
- **R3 — multithreaded BLAS nondeterminism.** Without thread-pinning, repeated default-`conv_tol`
  builds at the same geometry on this machine landed on visibly different SCF trajectories run to
  run (`chem-8tl`'s already-documented phenomenon, previously measured via Hamiltonian hash digest;
  here independently reproduced via HF population magnitude swinging ~1e-51 to ~1e-8 across runs at
  side 1.2). Not a correctness bug in this change; G6 sidesteps it by pinning
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`, the established
  convention from `test_hamiltonian_determinism.py`.

## 9. Deliverables

- `hybrid_quantum_solver/molecular_hamiltonian.py` — `_point_group` helper, `symmetric_scf` kwarg.
- `tests/test_symmetric_scf_spec.py` — G1–G7 gates.
- Results: `uv run pytest tests/test_symmetric_scf_spec.py -v` — see `sandbox-handoffs/chem-n3w.md`
  for the verbatim pass/fail line and the full regression sweep over every spec file this change
  could plausibly move.
