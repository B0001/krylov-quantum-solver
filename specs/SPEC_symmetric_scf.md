# SPEC: fix the SCF-convergence residue at its source — opt-in symmetric/tight SCF in the builder

**Status:** IMPLEMENTED. Gates in `tests/test_symmetric_scf_spec.py`. Results: §9. G8 (chem-9ud)
extends the flag to `build_dipole_operators`.

**Revised 2026-11-05 (chem-8ls, `SPEC_scf_stability_rule.md`):** §10: no lowest-level residue is left for the flag; G6 and G7 re-pinned and no longer Linux-only.

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
hybrid_quantum_solver.molecular_hamiltonian.build_dipole_operators(
    ..., conv_tol: float = 1e-9, symmetric_scf: bool = False) -> list[SparsePauliOp]   # chem-9ud
```

No other module changes. Since chem-9ud both builders run one private helper,
`_driver_problem` (the conv_tol floor plus the `pyscf.lib.num_threads(1)` pin), so a dipole built
with the Hamiltonian's arguments sits on the Hamiltonian's SCF orbitals. `MolecularHamiltonian.build_args["conv_tol"]` now records the *effective*
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

- **G8 — the dipole runs the Hamiltonian's SCF (chem-9ud).** On square H₄ at every G5 side,
  `build_molecular_hamiltonian` and `build_dipole_operators`, both with `symmetric_scf=True`, built
  at ambient `pyscf.lib.num_threads` 2 and then 8: all four `PySCFDriver` runs (spied in place)
  use the floored conv_tol 1e-13, run on one thread, and give bit-identical `mo_coeff`
  (`tobytes()` equality). G1 also pins the dipole's defaults to the Hamiltonian's.

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
  **macOS (2026-10-12, merging into main; macOS 27, Apple M3, threads pinned):** side 1.2 falls on
  the other side of that inversion — `scf_symmetry_status` returns `(True, -0.0764 Ha)` — so G6
  sees no `RuntimeWarning` and G7's `symmetric_scf=True` QKSD lands on B1g (`-4.3551` vs Ag
  `-4.2059` electronic), on the unmerged branch (`d9e9d2e`) as well as the merge. G6/G7 are
  therefore Linux-reference-only (`_LINUX_REFERENCE_ONLY`, as `SPEC_reachability_tolerance` §10);
  the other 21 gates pass on macOS. No macOS-native witness was searched for.
- **R3 — multithreaded BLAS nondeterminism.** Without thread-pinning, repeated default-`conv_tol`
  builds at the same geometry on this machine landed on visibly different SCF trajectories run to
  run (`chem-8tl`'s already-documented phenomenon, previously measured via Hamiltonian hash digest;
  here independently reproduced via HF population magnitude swinging ~1e-51 to ~1e-8 across runs at
  side 1.2). Not a correctness bug in this change; G6 sidesteps it by pinning
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`, the established
  convention from `test_hamiltonian_determinism.py`.
- **R4 — G8's thread leg can only fail where PySCF has OpenMP.** The macOS PySCF wheel is built
  without it: `lib.num_threads(n)` warns "OpenMP is not available" and the count stays 1. On macOS
  G8 therefore tests the floor only. Measured: a dipole that applies the floor but skips the pin
  passes G8 there. The pin half of chem-9ud reuses chem-7cx's pin, which was verified on Linux
  (`SPEC_scf_branch_determinism`). G8 has not yet been run on Linux.

## 9. Deliverables

- `hybrid_quantum_solver/molecular_hamiltonian.py` — `_point_group` helper, `symmetric_scf` kwarg.
- `tests/test_symmetric_scf_spec.py` — G1–G7 gates.
- Results: `uv run pytest tests/test_symmetric_scf_spec.py -v` — see `sandbox-handoffs/chem-n3w.md`
  for the verbatim pass/fail line and the full regression sweep over every spec file this change
  could plausibly move.
- **G8 (chem-9ud), macOS 27 arm64, 2026-11-05:** with the fix the file gives 27 passed and 2 skipped
  (G6/G7 are Linux-only). With the old dipole body put back (raw conv_tol, no pin), G8 fails at all
  6 sides. MO bits differ at 5 of 6: whole MO columns sign-flipped at 1.05/1.2/1.35 (max |ΔC|
  1.65/1.55/1.52), 6.9e-5 at 1.19 and 6.1e-9 at 1.4. At 1.1 the bits coincide and the conv_tol
  check (1e-9 vs 1e-13) catches it. The flips matter: the dipole that pre-fix callers got, compared
  with the one on the Hamiltonian's SCF, is off in spectral norm by ‖Δμ_x‖ = ‖Δμ_y‖ = 4.44, 4.99
  and 5.36 a.u. at 1.05, 1.2 and 1.35. At side 1.0 (not a G8 side), ‖Δμ_x‖ is 8.59 a.u.

## 10. chem-8ls (2026-11-05): the residue this flag removed is gone from the default build too

**What changed.** Every builder now follows internal SCF instabilities to a stable minimum (`SPEC_scf_stability_rule.md`). R1 above was half right.
- **What R1 got right.** Tightening `conv_tol` does not choose the RHF solution.
- **What R1 got wrong.** It called the symmetric solution and the broken one two branches the SCF may "land on". In fact the D2h-symmetric solution is a saddle point at every side, with orbital-Hessian eigenvalue −0.53 to −0.67 Ha.
- **Why R2's split looked inverted.** The free SCF stopped on that saddle at sides that depended on rounding: 1.15/1.19/1.20/1.25/1.30 on Linux against 1.05/1.10/1.19/1.35 on macOS, at `conv_tol` 1e-9.
- **Now.** The rule puts every side on the broken-symmetry minimum. Its lowest level is the B1g ground state, which HF genuinely overlaps (population 0.461 at side 1.2). So there is no lowest-level residue left for `symmetric_scf=True` to collapse, with or without the flag.

**Gates.**
- **G6 is re-pinned, and is no longer Linux-only.** At side 1.2, five builds each with the flag off and on: p₀ equal to 1e-4 and > 0.4, and no residue warning either way.
  - It was: without the flag, p₀ = 5.08e-9 and warned; with it, < 1e-20 (Linux).
  - Measured now: p₀ = 0.4612 both ways, on both platforms.
- **G7 is re-pinned, and is no longer Linux-only.** QKSD from either build sits on the B1g ground state (M = 28…32, within 1e-8; measured ≤ 1.2e-12 from M = 12 on macOS), and the reference starts there.
  - It was: the default build on B1g (called the artifact), and the flag on Ag.
- **G1–G5 and G8 are unchanged.**
  - G5 still pins what the flag is now for: a floor at `TIGHT_SCF_CONV_TOL` that removes the *excited*-level residue a default build can leave. `SPEC_scf_conv_tol` §11 shows one extra reachable level and Ritz leakage up to 3.7e-6 at the default, against ≤ 3.3e-14 tight.
  - G8 holds because the dipole and the Hamiltonian share `_driver_problem`, rule included.
  - G8's spy moved from `PySCFDriver.run` to `to_problem`. `_driver_problem` now calls `run_pyscf`, then the rule, then `to_problem`, so a `run` override recorded nothing and G8 failed 6/6 with the rule alone. `to_problem` sees the final orbitals, after the rule.

**What this does to the claim.**
- "QKSD lands on the physical level, not the artifact" identified the physical level with Ag, the level a *saddle-point* reference reaches. On the minimum the physical ground state is B1g, and every build reaches it.
- The flag's purpose narrows to the excited-level residue.
- §7's out-of-scope item "a true symmetry-constrained SCF" stays out of scope. It would select the saddle on purpose.

**Gate run.** See `SPEC_scf_stability_rule` §9.
