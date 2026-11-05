# SPEC: `_driver_orbitals`'s rebuild is the builder's own computation (chem-g14)

**Status:** REVISED 2026-11-05 (chem-g14, macOS arm64). Gate green: `tests/test_rebuild_gauge_search_spec.py`, 20/20.
- **What changed.** The rebuild now replays `build_args` through the builder's own SCF path. It reproduces `mh`'s operator bit-for-bit, so the reproducibility check is strict again (`_SAME_OPERATOR_ATOL = 1e-12`).
- **What it supersedes.** The first version (2026-10-12) worked around the symptom. It loosened the tolerance to 1e-4 and added an energy-filtered SCF retry and a per-orbital sign-gauge search. §10 records that design and why it went.
- **Why the file name stays.** It keeps chem-g14's history in one place. The original text is `git show 5460d0d:specs/SPEC_rebuild_gauge_search.md`.

---

## 1. Goal

**Claim:** on one host, `reachability._driver_orbitals(mh)` reproduces every non-identity coefficient of `mh.qubit_hamiltonian` **exactly**, not merely within an SCF tolerance, for any `mh` from `build_molecular_hamiltonian`. It does so with one SCF. The reason is that it calls the builder's own `molecular_hamiltonian._driver_problem` with the builder's own arguments, and maps with the builder's own `map_canonical`.

The claim is false if any of these happens:
- a gated configuration rebuilds to a nonzero residual;
- the rebuild takes more than one SCF;
- the rebuild stops reproducing `mh` after the builder's SCF path changes.

## 2. Background and honest framing

- **The bug (2026-10-09, Linux container, pyscf 2.13.1).** `_driver_orbitals` re-ran `PySCFDriver` and required the rebuilt operator to match to 1e-12. A second run landed about 4e-8 away, so every spatial-symmetry bit disappeared. That caused 17 failures in `tests/test_eigenstate_reachability_spec.py`.
- **Its source.** chem-7cx found it. PySCF's OpenMP integral evaluation is not bit-associative, so two runs on identical input followed different SCF trajectories. On square H4 they sometimes reached different branches. chem-7cx, chem-s6c, chem-0yh and chem-50w pinned `pyscf.lib.num_threads(1)` around every such SCF, the rebuild's included.
- **What was still wrong (this revision).** The rebuild was a hand-copy of the builder's recipe, not the recipe itself. It differed in three ways:
  - **(a)** It mapped with plain `mapper.map`, which sums terms in hash order. The builder uses `map_canonical`.
  - **(b)** Its `to_problem()` integral transform ran outside the thread pin.
  - **(c)** `build_args` did not record `symmetric_scf`.

  On this host, (a) alone leaves residuals up to 5.3e-15. (b) can reintroduce thread noise on an OpenMP host. And no change to the builder's SCF path would reach the copy. chem-8ls, the square-H4 branch choice, is the next bead and may well change that path.
- **The first version treated the symptom.** Its tolerance comment claimed 1e-4 had ">10x margin" over the 8.1e-5 it measured. The actual margin is 1.2x.
- **What we can claim if the gates pass.** On one host the rebuild reproduces the builder bit-for-bit for every gated configuration, with one SCF. It also follows changes to the builder's SCF path.
- **What we cannot claim.**
  - Cross-host or cross-version reproduction: an `mh` built on another machine is refused, so its diagnostic is not computed.
  - Linux: this revision has not been run there (§8 R1).
  - Which RHF branch square H4 *should* land on: that is chem-8ls.

## 3. Approach

1. `_driver_problem` returns its driver as well as its problem.
2. `_driver_orbitals` calls it with `build_args`: `atom, basis, charge, spin, conv_tol, symmetric_scf`.
3. It applies the builder's active-space transform and maps with `map_canonical`.
4. It compares non-identity coefficients against `_SAME_OPERATOR_ATOL = 1e-12`.

There is no energy check, no retry and no sign search. A deterministic recomputation cannot land on another branch or gauge, so those steps cannot change the answer.

**References.**
- The builder's own output, bit-for-bit.
- The separation measured in §9. Roundoff stays at or below 5.3e-15. A genuinely separate SCF trajectory to the same solution differs by at least 3.6e-10, whenever it differs above roundoff at all.

## 4. Public interface

No public change.
- `molecular_hamiltonian._driver_problem` (private) now returns `(problem, effective conv_tol, driver)`.
- `MolecularHamiltonian.build_args` gains `symmetric_scf`.
- Three private constants are removed from `reachability`: `_ENERGY_MATCH_ATOL`, `_MAX_SIGN_SEARCH_ORBITALS` and `_MAX_REBUILD_RETRIES`.

## 5. Acceptance criteria (validation gates)

The gates are in `tests/test_rebuild_gauge_search_spec.py`. It uses pyscf and qiskit-nature only, and runs in its own process.

- **G1: DEFINITION OF DONE, bit-for-bit.** With `_SAME_OPERATOR_ATOL` monkeypatched to 0, `_driver_orbitals(mh)` must not be None for each of 16 builder configurations:
  - square H4 at a ∈ {1.05, 1.10, 1.19, 1.35} × conv_tol ∈ {1e-6, 1e-9} (the G2 grid of the reachability gate file);
  - square H4 at a = 1.20 with the default conv_tol (that file's G4 broken-symmetry witness);
  - square H4 at a = 1.19 with conv_tol = 1e-13, and with `symmetric_scf=True`;
  - H2;
  - linear H4, singlet and triplet;
  - LiH CAS(2,5);
  - OH CAS((3,2),4).
- **G2: strict.**
  - `_SAME_OPERATOR_ATOL ≤ 1e-12`.
  - A 1e-10 edit to one coefficient is refused. 1e-10 is below the smallest separate-trajectory residual measured (3.6e-10).
  - An identity shift still verifies.
- **G3: one SCF.** A match and a refusal (a stale operator) each make exactly one `run_pyscf` call and one `to_problem` call. Before this revision a refusal made 5 and 40.
- **G4: the rebuild follows the builder.** Monkeypatch the builder's `PySCFDriver` to use `init_guess=ATOM`.
  - The changed build differs from the default build by more than the tolerance, so the gate is non-vacuous.
  - `_driver_orbitals` still verifies the changed build at the production tolerance.
- **G5: provenance is complete.** `build_args["symmetric_scf"]` records the argument.

> G1 is the definition of done. G4 is the claim chem-8ls relies on: its fix to the builder's SCF path needs no matching change here.

## 6. Implementation plan (as done)

1. Rewrote the gate file first. On the start point (`5460d0d`) it gave 20/20 RED:
   - every G1 case failed on the unsorted mapping's ≤ 5.3e-15 residual;
   - G2 failed on the 1e-4 tolerance;
   - G3 counted 5 SCFs and 40 transforms on a refusal;
   - G4 returned None, because the copy ignored the builder's initial guess;
   - G5 raised `KeyError`.
2. Changed the code:
   - `_driver_problem` returns the driver;
   - `build_args` records `symmetric_scf`;
   - `_driver_orbitals` rebuilds through both;
   - the tolerance is back to 1e-12;
   - the retry and the search are deleted.
3. Result: 20/20 GREEN in 31 s. The RED run took 7 min 18 s, most of it the retry and search loops grinding through refusals.

## 7. Out of scope

- **Which RHF branch square H4 lands on per platform.** That is chem-8ls. This revision makes `_driver_orbitals` follow whatever branch the builder picks. §9 has an initial-guess scan that bears on chem-8ls.
- **Cross-host reproduction of a serialized `mh`.** No such path exists in this repo.

## 8. Caveats and risks

- **R1: Linux not run.** This revision was verified on macOS arm64 only, and PySCF there has no OpenMP.
  - **The argument for Linux.** The rebuild is now literally a repeated builder call. chem-7cx measured repeated builds as bit-identical on Linux under the same pin (`SPEC_scf_branch_determinism.md` G1–G3). The `to_problem()` transform, previously unpinned in the rebuild, is now inside the pin.
  - **If G1 fails on Linux.** Then the bit-for-bit claim is false there, and that is the finding. Record it; do not loosen the tolerance to hide it.
- **R2: a refusal is now cheap and final.** One SCF settles every case that the builder's path does not reproduce: an edited operator, a different mapper, or an `mh` from another host. Each returns None, and the diagnostic disappears rather than being computed in a nearby basis. That is the intended contract: `orbital_parity_bounds` needs `mh`'s own MO basis.

## 9. Results (2026-11-05)

**Platform.** macOS 27.2, arm64, pyscf 2.13.1 (no OpenMP), qiskit-nature 0.8.0, numpy 2.4.6, scipy 1.15.3.

**Measurement.** This used a throwaway script, not committed, over 45 builder configurations:
- square H4 at a ∈ {1.00 … 1.40} (ten sides) × conv_tol ∈ {1e-6, 1e-9, 1e-13};
- three `symmetric_scf=True` builds;
- H2, stretched H2, linear H4, HeH+, LiH, LiH CAS(2,5), N2 CAS(6,6), the H4 triplet, OH, OH CAS((3,2),4), and the two half-detected geometries.

| Rebuild | Max non-identity residual vs `mh` |
|---|---|
| Start point (the copy, unsorted mapping) | 0 to 5.3e-15 (worst: OH). One SCF in all 45 cases: the retry and the search never fired. |
| This revision (`_driver_problem` + `map_canonical`) | **0.0 in all 45** |
| Separate SCF trajectory to the same solution (`init_guess=ATOM`, energy-matched to 1e-6, signs aligned in the AO metric) | conv_tol 1e-6: 1.6e-6 to 8.1e-5. conv_tol 1e-9: 3.6e-10 to 3.2e-6 (0.0 for H2 and stretched H2). conv_tol 1e-13: either ≤ 4e-16 or 4.2e-10 to 1.5e-9. |

Some configurations have no same-solution comparison:
- In 13 of the 33 square-H4 configurations, ATOM landed on the other RHF branch (ΔE ≈ 0.07 Ha).
- One more (a = 1.05 at 1e-6) missed the energy filter with ΔE = 1.5e-6.
- LiH, N2 and OH rotated within degenerate π shells.

The trajectory row reproduces the first version's Linux numbers: 8.1e-5 at 1e-6, and 2.4e-6 at 1e-9. Those were separate-trajectory differences, which the thread pin has since removed from the rebuild.

**For chem-8ls: the branch is an initial-guess property.** At conv_tol = 1e-9, the table gives the free RHF energy minus the D2h-symmetric RHF energy in Ha. A value near 0 means the symmetric branch.

| a (Å) | 1.00 | 1.05 | 1.10 | 1.15 | 1.19 | 1.20 | 1.25 | 1.30 | 1.35 | 1.40 |
|---|---|---|---|---|---|---|---|---|---|---|
| `minao` (the builder's guess) | −0.066 | 0 | 0 | −0.074 | 0 | −0.076 | −0.079 | −0.081 | 0 | −0.084 |
| `atom` | −0.066 | 0 | −0.072 | 0 | −0.076 | −0.076 | 0 | −0.081 | 0 | −0.084 |
| `1e` | −0.066 | 0 | −0.072 | 0 | −0.076 | −0.076 | −0.079 | 0 | 0 | 0 |

What the scan shows:
- At least one guess finds a broken-symmetry solution 0.066 to 0.084 Ha lower at every side except 1.05 and 1.35, where none of the three does.
- Whether the free SCF reaches it depends on the guess, and not monotonically in a.
- The reachability witnesses are exactly where `minao` stays symmetric (a = 1.10 and 1.19, with residue). With `atom` or `1e`, both break symmetry.

So a platform's rounding can plausibly tip the same guess into the other basin. A branch-pinning fix needs a rule, such as stability analysis or a symmetric SCF, not a different tolerance.

**Gate runs.** Each file ran in its own process (`uv run pytest <file> -q`), before (start point `5460d0d`) and after (this revision):

| Gate file | Before | After |
|---|---|---|
| `test_rebuild_gauge_search_spec` (this spec) | 5 passed (first version's G3–G5) | **20 passed** (G1–G5 above) |
| `test_eigenstate_reachability_spec` | 50 passed | 50 passed |
| `test_driver_orbitals_scf_determinism_spec` | 5 passed | 5 passed |
| `test_reachability_scf_determinism_spec` | 5 passed | 5 passed |
| `test_reachability_tolerance_spec` | 9 passed, 3 skipped | 9 passed, 3 skipped |
| `test_symmetry_reachability_spec` | 14 passed | 14 passed |
| `test_symmetric_scf_spec` | 27 passed, 2 skipped | 27 passed, 2 skipped |
| `test_scf_branch_determinism_spec` | 4 passed | 4 passed |
| `test_scf_conv_tol_spec` | 14 passed | 14 passed |
| `test_qksd_properties_spec`, `test_certified_dipole_spec`, `test_certified_dipole_noise_spec` | 4, 5, 4 passed | 4, 5, 4 passed |
| `test_certified_gaps_spec`, `test_excited_bounds_spec`, `test_hf_overlap_certificate_spec` | 4, 10, 23 passed | 4, 10, 23 passed |
| `test_odmd_spec`, `test_odmd_optical_spec` | 5, 4 passed | 5, 4 passed |
| Core: `test_reference_energies`, `test_krylov_convergence`, `test_noise_resilience`, `test_hamiltonian_determinism` | 21 passed | 21 passed |

- The skips are the same platform-conditional gates before and after.
- `uv run ruff check .` is clean.
- On this host the description's 17 Linux failures do not reproduce at the start point either. The revision does not change any outcome here. What it fixes is that the check is now correct by construction, and strict.

## 10. Superseded design (first version, 2026-10-12), kept for the record

**Its claim.** A rebuild's mismatch had two separable causes:
- **(A)** free RHF converging to a worse local minimum;
- **(B)** free RHF landing on an energy-degenerate alternate MO gauge that differs by per-orbital signs.

**Its mechanism.** Three pieces worked around those causes:
- an energy pre-filter with up to `_MAX_REBUILD_RETRIES = 5` fresh SCFs, for (A);
- a search over the 2^(k−1) sign patterns of the active MOs, up to `_MAX_SIGN_SEARCH_ORBITALS = 10`, for (B);
- the tolerance loosened to 1e-4 to absorb the residual left after gauge matching.

**Why it went.** chem-7cx and chem-50w later showed that (A) and (B), and the residual, all came from one source: unpinned OpenMP threads giving every rebuild its own SCF trajectory.
- With the threads pinned, the rebuild is a deterministic recomputation. On this host the copy never retried or searched in any of the 45 configurations.
- With the shared path, nothing is left for the mechanism to absorb.
- The mechanism also had costs:
  - 5 SCFs and 40 transforms per refusal at k = 4, and up to 5 × 1023 transforms at k = 11;
  - a tolerance 8 orders of magnitude looser than the residual it was guarding, with a mis-stated margin;
  - an energy filter that accepted any operator whose `e_tot` matched.

Its other findings stand as recorded at the time:
- the 14–15 then-remaining failures traced to chem-7cx and chem-3fd;
- the signed-identity AO-metric overlap of gauge twins.
