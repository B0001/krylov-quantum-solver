# SPEC: the builder's RHF reference is a stable minimum, on every platform (chem-8ls)

**Status:** IMPLEMENTED 2026-11-05. Gates: `tests/test_scf_stability_rule_spec.py`. Results: §9.

---

## 1. Goal

**Claim.** Square H₄ (STO-3G) has two RHF stationary points at every side from 1.00 to 1.40 Å:
- the D2h-symmetric solution, which is a **saddle point**;
- a broken-symmetry **minimum**, 0.066–0.084 Ha lower.

A free SCF started from PySCF's symmetric `minao` guess keeps the symmetry up to rounding. So whether it stops on the saddle or escapes to the minimum depends on the platform's rounding. The fix is a branch rule, not a tolerance: `molecular_hamiltonian._driver_problem` follows any internal SCF instability to a stable minimum before the integral transform. Every builder routes through it: `build_molecular_hamiltonian`, `build_dipole_operators`, the rebuild in `reachability._driver_orbitals`, and since this change `reachability.scf_symmetry_status`.

The claim is false if any of these happens:
- a side where the symmetric RHF is itself a minimum (the rule would then be choosing between two minima, which is a physics call);
- a platform whose rule-built reference differs from the recorded table;
- a reference that was already stable changes under the rule.

## 2. Background and honest framing

- **The bug (chem-8ls).** Four gate files pinned the RHF branch at square-H₄ witness sides. Those branches had been recorded on one platform and did not reproduce on another:
  - `test_reachability_tolerance_spec`, `test_eigenstate_reachability_spec`, `test_scf_conv_tol_spec` and `test_symmetry_reachability_spec`;
  - example: at a = 1.1 the free SCF gave `scf_symmetry_status = (False, -5e-15)` on macOS and `(True, -0.0716)` on Linux.
- **Earlier evidence.**
  - `SPEC_rebuild_gauge_search` §9 showed the branch is an initial-guess property, non-monotonic in a.
  - `SPEC_symmetric_scf` R2 recorded the Linux pattern as "essentially inverted" from macOS.
  - This spec explains both. It is one saddle and one minimum at every side; which one the free SCF stops on is a rounding accident.
- **Why stability analysis, not a symmetric SCF.**
  - The symmetric solution is not a valid RHF minimum. Its lowest orbital-Hessian eigenvalue is −0.53 to −0.67 Ha (§9), far beyond any rounding.
  - The minimum is variationally lower.
  - |HF⟩ on the minimum overlaps the true (B1g) ground state with population 0.44–0.48. On the saddle that overlap is zero by symmetry, apart from SCF residue (`SPEC_reachability_tolerance` §2b).
  - For a solver that starts Krylov from |HF⟩, the minimum is the right reference. Stability analysis with descent along the unstable mode is the standard way to reach it.
  - A symmetry-adapted SCF by default would deliberately select the saddle. It would also run in PySCF's re-oriented symmetry frame (the reason `orbital_parity_bounds` maps back to the driver's frame), which for any basis beyond s functions changes the AO basis that dipole components are expressed in.
  - `symmetric_scf=True` stays what it was, a `conv_tol` floor (`SPEC_symmetric_scf`). It now applies on top of the rule.
- **What we can claim if the gates pass.**
  - On square H₄, 1.00–1.40 Å, at `conv_tol` 1e-6, 1e-9 and 1e-13, the builder's reference is the same stable minimum on macOS arm64 and on Linux x86_64. E_HF agrees to 1e-10 Ha and p₀ to 8 digits at `conv_tol` ≤ 1e-9.
  - On the ordinary gated systems the rule changes nothing, bit for bit.
- **What we cannot claim.**
  - **Uniqueness.** Square H₄ has two mirror-image broken-symmetry minima. Measured at a = 1.10, stepping from the saddle in the two directions of its unstable mode gives bonds along the x edges (density-matrix elements P₀₁ = P₂₃ = 0.89) or along the y edges. Both are stable, both have E = −1.7825511826 Ha, and the densities differ by 1.3. Which one the rule reaches can depend on the sign of a Davidson eigenvector, so the operator's coefficients are not claimed bit-identical across platforms. They never were.
  - **Global minimality.** The rule is a local descent, so it reaches *a* stable minimum, not necessarily the global one.
  - **External instabilities.** RHF→UHF and real→complex instabilities are not followed: the qubit Hamiltonian is built on restricted, real orbitals.
  - **Excited-level residues.** A loose SCF can still leave a small residue on symmetry-forbidden *excited* levels of the minimum's own symmetry, and which geometries carry one is still trajectory-dependent (§9). The rule fixes the lowest level only.

## 3. Approach

`_follow_instabilities(mf, max_steps=5)`:
1. Runs PySCF's internal stability check (`stability.rhf_internal`, or `rohf_internal` for ROHF) with `with_symmetry=False`, so the Davidson start vector can break point-group symmetry.
2. If the check finds the solution unstable, re-runs the SCF from the rotated orbitals it returns.
3. Repeats until the solution is stable. It raises `RuntimeError` if that takes more than `max_steps`, rather than return a saddle silently.

A stable solution is only read, never moved, so on stable systems the build is bit-identical to the pre-rule builder. The call sits inside the `pyscf.lib.num_threads(1)` pin, between `run_pyscf()` and `to_problem()`.

**References.**
- The **dense** real orbital Hessian, from `newton_ah.gen_g_hop_rhf` applied to unit vectors. It is independent of the Davidson solver the rule uses.
- The symmetry-enforced RHF (`symmetry=True`). Its irrep-blocked SCF cannot leave the symmetric stationary point, so it is the saddle on every platform.

## 4. Public interface

No public signature changes.
- `molecular_hamiltonian._follow_instabilities(mf, max_steps=5)` is new and private.
- `reachability.scf_symmetry_status(atom, basis, conv_tol)` now reads the builder's SCF (`_driver_problem`) instead of running its own free RHF. Its answer describes the reference `build_molecular_hamiltonian` builds.

## 5. Acceptance criteria (validation gates)

`tests/test_scf_stability_rule_spec.py`, pyscf and qiskit-nature only, about 1 min.

- **G1 — one branch on every platform (DEFINITION OF DONE).** The sweep is a ∈ {1.00, 1.05, 1.10, 1.15, 1.19, 1.20, 1.22, 1.25, 1.30, 1.35, 1.40} × `conv_tol` ∈ {1e-6, 1e-9, 1e-13}. For each build:
  - E_HF equals the recorded minimum (§9) within max(1e-8, 10·conv_tol) Ha;
  - p₀ (the lowest level's HF population) equals the recorded value within 1e-4.
- **G2 — the reason for the rule.** At every side:
  - the symmetry-enforced RHF has lowest dense orbital-Hessian eigenvalue < −0.5 Ha (a saddle);
  - the builder's reference has lowest eigenvalue > +0.5 Ha (a minimum);
  - the saddle lies 0.06–0.09 Ha above the minimum.
- **G3 — what the rule does, and where it does nothing.**
  - **(a) Platform-independent non-vacuity.** At every side, a free SCF started *on* the saddle stays there to 1e-10 Ha, and `_follow_instabilities` takes it to the recorded minimum within 1e-8 Ha.
  - **(b) Stable references are untouched.** On H₂ (0.74, 2.0), linear H₄ singlet and triplet, HeH⁺, LiH, LiH CAS(2,5), N₂ CAS(6,6) and the OH doublet, a build with the rule gives the same operator as a build with `_follow_instabilities` replaced by a no-op (the pre-rule builder). The Pauli labels and the coefficient arrays must be identical, and so must E_HF.
- **G4 — failure is loud.** If an instability survives `max_steps` re-runs, `_follow_instabilities` raises `RuntimeError` after exactly that many re-runs.

G1 is the definition of done.

## 6. Implementation plan (as done)

1. Measured the stability picture on both platforms (§9) before choosing the rule.
2. `_follow_instabilities` in `_driver_problem`; `scf_symmetry_status` through `_driver_problem`.
3. Wrote this gate file. With the rule switched off (a pytest plugin replacing `_follow_instabilities` with a no-op) it gives **32 failed / 33 passed** on macOS:
   - G1 fails at exactly the five macOS saddle sides 1.05, 1.10, 1.19, 1.22 and 1.35, at all three `conv_tol`s;
   - G2 fails at the same five sides;
   - G3(a) fails at all 11 sides;
   - G4 fails.
   
   With the rule it gives **65 passed**.
4. Re-pinned the four gate files to the branch the rule selects, and recorded the old values next to the new ones in each spec:
   - `SPEC_reachability_tolerance` §11;
   - `SPEC_eigenstate_reachability` §9;
   - `SPEC_scf_conv_tol` §11;
   - `SPEC_symmetry_reachability` §10.
5. Re-pinned the other gates that pinned the old branch (§9).

## 7. Out of scope

- External (RHF→UHF, real→complex) stability.
- Choosing between symmetry-equivalent minima.
- Systems beyond the gated set. The rule runs on every build, so any other system whose free SCF stops on a saddle now changes. That is the intended effect, but no other system was measured.
- Excited-level SCF residues at loose `conv_tol` (§2): `TIGHT_SCF_CONV_TOL` still removes them.

## 8. Caveats and risks

- **R1 — the threshold.** PySCF calls a solution unstable below a Hessian eigenvalue of −1e-5. A system with a near-zero mode, such as a continuous degeneracy, sits near that line. None of the gated systems does: their eigenvalues are either ≥ +0.5 or ≤ −0.5 on square H₄, and the ordinary set is untouched (G3b).
- **R2 — cost.** One Davidson on the orbital Hessian per build, plus one SCF where the rule fires. This is negligible at the gated sizes. It has not been measured at scale.
- **R3 — chemcheck hashes.** The frozen tiers are H₂, linear H₄, LiH and N₂. All of them are in G3b's untouched set or are the same kind of stable closed-shell system, so no hash should move. The hashes are Linux-only, and Linux was not re-run on chemcheck here.

## 9. Results (2026-11-05)

**Platforms.**
- macOS 27.2 arm64: pyscf 2.13.1 without OpenMP, numpy 2.4.6, scipy 1.15.3 (Accelerate), qiskit-nature 0.8.0.
- Linux x86_64: the `claude-chem` image under emulation, with OpenMP and 8 CPUs. pyscf 2.13.1, numpy 2.4.6, scipy 1.15.3.

**The stability picture.** These are scratch scans, not committed; G2 and G3 re-derive what they assert. Free RHF from `minao` at `conv_tol` 1e-9 on one thread. λ_min is the lowest eigenvalue of the dense orbital Hessian.

| a (Å) | macOS free E (Ha) | macOS λ_min | Linux free E (Ha) | Linux λ_min | D2h saddle E | saddle λ_min | minimum E | minimum λ_min |
|---|---|---|---|---|---|---|---|---|
| 1.00 | −1.76107505 | +0.529 | −1.76107505 | +0.529 | −1.69488959 | −0.529 | −1.76107505 | +0.529 |
| 1.05 | **−1.70614482** | **−0.552** | −1.77510448 | +0.551 | −1.70614482 | −0.552 | −1.77510448 | +0.551 |
| 1.10 | **−1.71095268** | **−0.573** | −1.78255118 | +0.571 | −1.71095268 | −0.573 | −1.78255118 | +0.571 |
| 1.15 | −1.78463922 | +0.590 | **−1.71054711** | **−0.593** | −1.71054711 | −0.593 | −1.78463922 | +0.590 |
| 1.19 | **−1.70714947** | **−0.608** | **−1.70714947** | **−0.608** | −1.70714947 | −0.608 | −1.78312491 | +0.605 |
| 1.20 | −1.78236132 | +0.608 | **−1.70593123** | **−0.611** | −1.70593123 | −0.611 | −1.78236132 | +0.608 |
| 1.25 | −1.77653013 | +0.624 | **−1.69792788** | **−0.629** | −1.69792788 | −0.629 | −1.77653013 | +0.624 |
| 1.30 | −1.76781622 | +0.639 | **−1.68721696** | **−0.645** | −1.68721696 | −0.645 | −1.76781622 | +0.639 |
| 1.35 | **−1.67436361** | **−0.659** | −1.75677660 | +0.651 | −1.67436361 | −0.659 | −1.75677660 | +0.651 |
| 1.40 | −1.74387673 | +0.661 | −1.74387673 | +0.661 | −1.65983984 | −0.672 | −1.74387673 | +0.661 |

Bold marks the free SCF stopped on the saddle: macOS at 1.05/1.10/1.19/1.35, and Linux at 1.15/1.19/1.20/1.25/1.30. This scan omits 1.22, where it stops on the saddle on both platforms (the RED runs below).
- The saddle and minimum energies agree between platforms to 8 digits.
- Only the stopping point differs, and one step along the unstable mode always reaches the minimum.
- At 1.05 and 1.35 the minimum exists although `SPEC_rebuild_gauge_search` §9's three-guess scan never found it. That scan sampled stopping points, not stationary points.

**The rule's reference (G1's table).** E_HF and p₀ at `conv_tol` = 1e-13, macOS and Linux identical to 1e-10 Ha and 8 digits:

| a (Å) | 1.00 | 1.05 | 1.10 | 1.15 | 1.19 | 1.20 | 1.22 | 1.25 | 1.30 | 1.35 | 1.40 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| E_HF (Ha) | −1.7610750541 | −1.7751044780 | −1.7825511826 | −1.7846392221 | −1.7831249107 | −1.7823613168 | −1.7804148460 | −1.7765301318 | −1.7678162175 | −1.7567766045 | −1.7438767269 |
| p₀ | 0.47587995 | 0.47270840 | 0.46922809 | 0.46541280 | 0.46210250 | 0.46123748 | 0.45946120 | 0.45667855 | 0.45171429 | 0.44632542 | 0.44049595 |

- p₀ falls smoothly and monotonically with a.
- At `conv_tol` = 1e-6, E_HF stays within 3.1e-7 Ha of the table on both platforms (macOS at 1.00, Linux at 1.05), and p₀ within 1.2e-6.

**Where the rule fires.** This is a scratch spy over the builder.
- **Square H₄ at `conv_tol` ∈ {1e-6, 1e-9, 1e-13}.** On macOS, exactly one extra SCF at a = 1.05/1.10/1.19/1.22/1.35 and none elsewhere.
- **No-op (zero extra SCFs).** H₂ (0.74, 2.0), linear H₄ (0.9, 1.0, 2.0), the H₄ triplet, HeH⁺, LiH, LiH CAS(2,5), N₂ CAS(6,6), OH, OH CAS((3,2),4) and linear H₆.

**What is left at loose `conv_tol`: excited-level residues, trajectory-dependent.** Default-`conv_tol` Ritz leakage at depth 6 (the `SPEC_scf_conv_tol` table), macOS → Linux:

| a | 1.00 | 1.05 | 1.19 | 1.20 | 1.35 |
|---|---|---|---|---|---|
| default | 1.4e-14 → **1.6e-7** | **2.7e-7** → 8.1e-15 | **3.0e-7** → **3.0e-7** | **5.1e-7** → **2.9e-7** | **3.7e-6** → **2.0e-7** |
| tight (1e-13) | ≤ 3.3e-14 on both | | | | |

- On the macOS default build at a = 1.1, one excited level sits between the 1e-10 and 1e-8 cuts (6.3e-10); Linux has none.
- None of this touches the lowest level, and `TIGHT_SCF_CONV_TOL` removes all of it on both platforms. So the gates pin lowest-level facts and tight-`conv_tol` leakage, never a default-`conv_tol` residue.

**Gate runs (2026-11-05/06, tree `f6a7c3f`).**

*macOS 27.2 arm64, the full gate suite* (`GATE_NO_CACHE=1 GATE_JOBS=5 bash scripts/run_gates.sh`, i.e. `make gates-fresh`, each gate file in its own process):

| tree | passed | failed | skipped (optional dependency) |
|---|---|---|---|
| start point `bb2bba7` | 118 | 0 | 4 |
| the rule alone, gates not re-pinned | 111 | 7: eigenstate_reachability, chained_overlap, odmd, reachability_tolerance, scf_conv_tol, symmetric_scf, symmetry_reachability | 4 |
| this tree (122 files + this gate) | **119** | **0** | 4 |

The 7 failures are exactly the files re-pinned here (§6). Also on macOS:
- The non-spec half of `make test` gives 60 passed and 4 skipped.
- `tests/test_dmrg_reference.py`, in its own process, gives 2 passed and 1 skipped (the missing-dependency path, since block2 is installed).
- `uv run ruff check .` is clean.

*Repeats: no run-to-run variation.* These are the four bead files plus this gate.

| platform | runs | reachability_tolerance | eigenstate_reachability | scf_conv_tol | symmetry_reachability | scf_stability_rule |
|---|---|---|---|---|---|---|
| macOS (full sweep + 5 rounds, each file its own process) | 6 | 11 passed every run | 57 | 16 | 18 | 65 |
| Linux x86_64 (`claude-chem` image under emulation, OpenMP on, 8 threads; mayor's `chem-linux-tests.sh`, a fresh container and one pytest process per round) | 5 | 167 passed every run, the same five files together | | | | |

*Linux, once, the other re-pinned files.* The run covered symmetric_scf, odmd, scf_branch_determinism, reachability_scf_determinism, rebuild_gauge_search, driver_orbitals_scf_determinism, hf_overlap_subspace, chained_overlap and chemcheck, plus the four core tests (reference_energies, krylov_convergence, noise_resilience, hamiltonian_determinism). Result: 198 passed and 3 failed.
- **The 3 failures are pre-existing and unrelated to this change.** They are `test_chemcheck_spec`'s frozen-tier checks: the hash for T1 (H₄ chain) and T3 (N₂), and the CX count for T3.
- **Why they are not from the rule.** With the rule switched off, the same 3 fail. A spy shows the rule takes zero extra SCF runs on every chemcheck tier on Linux, so those builds are bit-identical with and without it.
- **Not pursued here.** The hashes are platform-specific by design (`chemcheck/tiers.py`). This container is not the platform that froze them.

*RED.* With the rule switched off this gate fails:
- **macOS, 32 of 65:** G1/G2 at 1.05, 1.10, 1.19, 1.22 and 1.35.
- **Linux, 36 of 65:** G1/G2 at 1.15, 1.19, 1.20, 1.22, 1.25 and 1.30.
- **On both:** all of G3(a), and G4.

Those are each platform's own saddle stops. They share only 1.19 and 1.22, which is the cross-platform divergence chem-8ls reported.
