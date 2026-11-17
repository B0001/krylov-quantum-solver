# chem-1x3 handoff — ODMD ref should witness HF-unreachable-ground geometries

**Bead:** chem-1x3 (P3, task) — "odmd.py's `reach[0]` vs `w_eig[0]` fix is currently untested."

**Resolution:** added a witness gate that separates the two formulas, and (see §0) re-applied the
one-line fix itself, because this branch did not actually contain it.

---

## 0. A correction to the bead's premise — please read before merging

The bead's description says "chem-3ds fixed `odmd.py`'s `ODMDProblem.ref`... from
`float(w_eig[0] - mu)` to `float(reach[0] - mu)`." **That fix is not on `main`, and was not on this
branch either when I started.** It exists only on the sibling branch `sandbox/chem-3ds`
(commits `33aaa43`, `861b80a`), which has not been merged:

```
$ git merge-base main sandbox/chem-3ds
9a3ac3516c35f9a2e0458f87487aca6bed4770ca     # == main's tip — chem-3ds's commits sit ONLY on its own branch
$ git show HEAD:odmd.py | grep 'ref='        # before I touched anything, on sandbox/chem-1x3 (cut from the same point)
                       offset=mh.energy_offset + mu, ref=float(w_eig[0] - mu))   # the OLD, buggy formula
```

So the bead as literally scoped ("add a gate for an existing fix") was not satisfiable without
first re-applying the fix — there was nothing yet to gate. I re-applied the identical one-line
change (`odmd.py`'s `ref=float(w_eig[0] - mu)` → `ref=float(reach[0] - mu)`), since it is required
infrastructure for the new gate to assert anything true, and it is the exact, already-reviewed
change from the closed `chem-3ds` bead — not new scope.

**Consequence for the human merging these branches:** `sandbox/chem-3ds` and `sandbox/chem-1x3`
both touch the same one line of `odmd.py` identically. Merge either one first; the second merge
will conflict trivially (identical hunk) and can be resolved by taking either side. If
`sandbox/chem-3ds` merges first, this branch's `odmd.py` edit becomes a no-op diff against `main`
(safe to drop) and only the new gate + spec addendum need to land.

---

## 1. What changed

### 1a. `odmd.py` — the fix (re-applied; see §0)

```diff
-    energy. ``ref`` is the exact (FCI) ground energy in the frame, so ``ref + offset`` is the
-    total FCI energy.
+    energy. ``ref`` is the exact HF-reachable ground energy in the frame (the lowest
+    population-admitted eigenvalue, `reachable_mask`'s cut)...
...
-    ref: float             # exact ground energy in the centered frame
+    ref: float             # HF-reachable ground energy in the centered frame
...
-                       offset=mh.energy_offset + mu, ref=float(w_eig[0] - mu))
+                       offset=mh.energy_offset + mu, ref=float(reach[0] - mu))
```

Also updated the docstring (it still said "exact (FCI) ground energy," which was never quite right
even pre-fix in spirit — see §1c for why `reach[0]` is the physically correct quantity, not just a
relabeling).

### 1b. `tests/test_odmd_spec.py` — new gate G5

`test_G5_ref_is_hf_reachable_not_global_fci_ground`. Witness: square H4, side 1.19 Angstrom,
`conv_tol=TIGHT_SCF_CONV_TOL` (`reachability.TIGHT_SCF_CONV_TOL = 1e-13`) — the same
symmetry-forbidden-ground family as `tests/test_eigenstate_reachability_spec.py`'s G9 control
(a tight, symmetric SCF makes the true FCI ground genuinely unreachable from the HF determinant,
rather than reachable-via-SCF-residue as in G9's main witnesses). Asserts, in order:

1. The global FCI ground (`w_eig[0]`) has HF population < 1e-20 (numerically zero) — confirms this
   really is an HF-unreachable-ground witness, not a vacuous check.
2. The HF-reachable ground (`reach.min()`, the population cut at 1e-8) sits > 1e-2 Ha above it —
   confirms it is a genuinely different level, not the same value by coincidence.
3. `build_odmd_problem(mh).ref == reach.min() - mu` (abs 1e-8) — the fix.
4. `build_odmd_problem(mh).ref != w_eig[0] - mu` (within 1e-3) — rules out the old formula.
5. **The empirical check that makes this more than a relabeling assertion:** the actual noiseless
   survival-amplitude signal (`odmd_energy(prob.s, prob.tau)`) converges to `ref` to < 1e-5 Ha,
   while sitting > 1e-2 Ha away from the old (global-ground) value. This shows `ref=reach[0]-mu` is
   not an arbitrary choice between two numbers — it is the value the signal the method actually
   measures can converge to, because `s_k = <phi_0|e^{-ik tau H}|phi_0>` has no amplitude on a state
   `|phi_0>` doesn't overlap.

**chem-7cx nondeterminism:** square H4's free-RHF branch is nondeterministic across ambient OpenMP
thread counts on near-degenerate geometries (persistent memory `chem-7cx-root-cause`). Confirmed
this directly: the *unpinned* build sometimes lands on a solution where HF reaches the global
ground (pop[0] ≈ 0.46, same value for both formulas — a useless witness) and sometimes on the
genuinely-unreachable one (pop[0] ≈ 2e-32), depending on ambient thread count/run. The test pins
`pyscf.lib.num_threads(1)` around the `build_molecular_hamiltonian` call (the same pattern as
`benchmark_nbn.py`'s `_tight_scf` / `tests/test_nbn_scf_determinism_spec.py`) and this was
reproducible across 3+ repeated runs once pinned. Side 1.10 (the G9 witness's other geometry) does
**not** give a clean separation at `TIGHT_SCF_CONV_TOL` even pinned — this matches an existing,
pre-existing failure in `tests/test_eigenstate_reachability_spec.py::test_G9_control_tight_scf_stays_on_the_ag_level[1.1]`
on this platform (see §4), so I used 1.19 only, which gives a clean, reproducible witness.

### 1c. `specs/SPEC_odmd.md` — addendum + new acceptance criterion G5

Added a short "Addendum (chem-1x3)" paragraph under the status line explaining why G1–G4 never
exercised the `reach[0]` vs `w_eig[0]` distinction, and a new **G5** bullet in §5 (Acceptance
criteria) stating the witness, the claim, and the measured separation.

---

## 2. Gate-run evidence (verbatim)

```
$ uv run pytest tests/test_odmd_spec.py -v
tests/test_odmd_spec.py::test_G1_overlaps_alone_recover_fci PASSED
tests/test_odmd_spec.py::test_G2_depth_convergence_and_no_variational_floor PASSED
tests/test_odmd_spec.py::test_G3_matched_budget_beats_kqd PASSED
tests/test_odmd_spec.py::test_G4_svd_truncation_is_the_mechanism PASSED
tests/test_odmd_spec.py::test_G5_ref_is_hf_reachable_not_global_fci_ground PASSED
5 passed, 4 warnings in 242.00s (0:04:01)
```

**RED confirmed before the fix:** while developing the test, a bug in my own comparison (comparing
`prob.ref`, a *centered-frame* value, against `reach.min()`, an *uncentered* one) produced a clean
failure (`-0.9329274083058152 == -4.22466888980546`) that made me go back and get the frames right
— not evidence of the real RED/GREEN transition. I did separately confirm the real transition: with
`odmd.py` reverted to `ref=float(w_eig[0] - mu)` (`git stash` the fix only), G5 fails on assertion 3
(`prob.ref == pytest.approx(reach.min() - mu, ...)`, obtained `w_eig[0] - mu ≈ -1.0825` vs expected
`≈ -0.9329`) and assertion 5 (the noiseless signal matches `reach[0]-mu` to 1e-5 Ha, which is off by
> 0.1 Ha from `w_eig[0]-mu`). Restored the fix before finishing.

Full consumer sweep — every other site that builds an `ODMDProblem` or reads `.ref`, confirming the
fix is a no-op everywhere else (as the bead itself predicted: none of these systems separate
`reach[0]` from `w_eig[0]`):

```
$ uv run pytest tests/test_odmd_excited_spec.py tests/test_adaptive_shots_spec.py \
                 tests/test_odmd_uq_spec.py tests/test_device_odmd_spec.py -v
tests/test_odmd_excited_spec.py::test_G1_same_signal_recovers_the_spectrum PASSED
tests/test_odmd_excited_spec.py::test_G2_depth_is_the_excited_resource PASSED
tests/test_odmd_excited_spec.py::test_G3_noisy_gap_beats_qksd_at_matched_noise PASSED
tests/test_odmd_excited_spec.py::test_G4_threshold_mechanism_and_visibility_law PASSED
tests/test_adaptive_shots_spec.py::test_G1_strict_budget_preservation PASSED
tests/test_adaptive_shots_spec.py::test_G2_uniform_is_optimal_on_an_undamped_signal PASSED
tests/test_adaptive_shots_spec.py::test_G2b_geometric_whitening_preserves_eigenphases PASSED
tests/test_adaptive_shots_spec.py::test_G3_optimal_decay_tracks_the_damping_rate PASSED
tests/test_adaptive_shots_spec.py::test_G4_the_win_is_1p3x_in_error_and_zero_in_variance PASSED
tests/test_odmd_uq_spec.py::test_G1_union_coverage_meets_nominal_everywhere PASSED
tests/test_odmd_uq_spec.py::test_G2_width_is_informative PASSED
tests/test_odmd_uq_spec.py::test_G3_union_is_load_bearing PASSED
tests/test_odmd_uq_spec.py::test_G4_bias_is_invisible_the_boundary PASSED
tests/test_device_odmd_spec.py::test_G1_exact_channel_immunity PASSED
tests/test_device_odmd_spec.py::test_G2_damped_plus_shot_noise_beats_kqd PASSED
tests/test_device_odmd_spec.py::test_G3_modulus_window_is_the_mechanism PASSED
tests/test_device_odmd_spec.py::test_G4_aer_phase_survival_and_boundary PASSED
17 passed, 8 warnings in 464.64s (0:07:44)
```

(`trotter_odmd.py` already uses `reach[0]`, confirmed unchanged by `git log -p -- trotter_odmd.py`;
`rodeo.py` doesn't consume `odmd.ODMDProblem` at all — chem-3ds's other fix, out of this bead's
scope. Neither re-run here since neither can be affected by this diff.)

Lint:

```
$ uv run ruff check odmd.py tests/test_odmd_spec.py specs/SPEC_odmd.md
All checks passed!
```

**Not run this session:** `make gates` / the full `test_*_spec.py` sweep (~90 spec files, many
minutes each, block2-isolated). Scoped to the files this bead's diff can affect (above), per the
compute-envelope guidance — a full sweep would re-verify ~90 unrelated gates for a one-line fix
with a known, narrow blast radius (`grep -rl "build_odmd_problem\|ODMDProblem"` — the exact list
run above, plus `trotter_odmd.py`/`rodeo.py` addressed by inspection instead of re-running).

**Pre-existing, unrelated failures noticed in passing (not touched, not in scope):** running
`tests/test_eigenstate_reachability_spec.py` in full (baseline check before writing my own test,
to understand its square-H4 conventions) showed **14 failures / 35 passed** on this Linux
container, including `test_G9_control_tight_scf_stays_on_the_ag_level[1.1]` (not `[1.19]`, which
passes) — this is exactly the same platform-dependent-residue phenomenon the spec documents for
macOS in reverse (`SPEC_eigenstate_reachability.md` §2: "this is why `test_reachability_tolerance_spec`
G1-G3 fail on macOS"). I did not investigate or fix these; they predate this session and are outside
chem-1x3's scope. Worth a bead if not already filed — I did not search for one.

---

## 3. Every number, with regenerating command and where it's vendored

| Number | Value | Command | Vendored |
|---|---|---|---|
| Global FCI ground, sq-H4(1.19), `conv_tol=1e-13` | −4.374273658905307 Ha (electronic) | `tests/test_odmd_spec.py::test_G5...` (`w_eig[0]`, printed in assertion message on failure) | not vendored — recomputed every gate run |
| HF-reachable ground, same geometry | −4.224668889805460 Ha (electronic) | same test (`reach.min()`) | not vendored |
| HF population of the global ground | ≈ 2.03e-32 | same test (`pops[0]`) | not vendored |
| `ODMDProblem.ref` (after fix) | −0.9329274083058152 Ha (centered frame) | same test (`prob.ref`) | not vendored |
| Noiseless ODMD estimate (K=20) at this witness | −0.9329274083058098 Ha | same test (`odmd_energy(prob.s, prob.tau)`) | not vendored |

All five numbers come directly from the new gate, which re-derives them every run (per
`specs/README.md`'s convention of never only putting a number in `data/`) — nothing here depends on
a committed checkpoint or table.

---

## 4. Pre-registered criteria and pass/fail

From the bead: "Add a gate... asserting `ODMDProblem.ref` equals the HF-reachable eigenvalue (not
the lower, HF-unreachable global FCI ground)... on a geometry where HF cannot reach the global
ground."

1. **A geometry where HF genuinely cannot reach the global ground** — **PASS.** sq-H4(1.19),
   `TIGHT_SCF_CONV_TOL`, pop[0] < 1e-20 (§1b, assertion 1).
2. **`ODMDProblem.ref` equals the HF-reachable eigenvalue** — **PASS** (§1b, assertion 3).
3. **...not the global FCI ground** — **PASS**, both as a direct inequality (assertion 4) and
   empirically via the signal itself (assertion 5) — the stronger, non-vacuous form.
4. **Test suite passes** — **PASS** for every file this diff can affect (§2); full `make gates` not
   run (see §2's note — out of scope cost for a one-line, narrow-blast-radius fix).

---

## 5. What was decided not to do, and why

- **Did not attempt a clean witness at side 1.10.** It fails the analogous check even in the
  existing, independently-written `test_G9_control_tight_scf_stays_on_the_ag_level[1.1]` on this
  platform — adding a second, currently-broken witness would just inherit that pre-existing
  flakiness rather than add coverage. 1.19 alone is sufficient to satisfy the bead's claim
  (one clean witness, not an exhaustive sweep).
- **Did not fix the 14 pre-existing `test_eigenstate_reachability_spec.py` failures** noticed while
  reading that file for convention — out of this bead's scope (it's about `odmd.py`, not
  `reachability.py`/the symmetry diagnostic), and I did not reproduce or diagnose them beyond
  noting they are platform-dependent in the same way the spec already documents for macOS.
- **Did not add a `specs/BACKLOG.md` entry.** This isn't a new open hypothesis — it closes out a
  narrow, pre-registered follow-up (chem-3ds) against an already-`IMPLEMENTED` spec; recorded
  directly in `SPEC_odmd.md`'s own addendum + G5 instead.
- **Did not touch `trotter_odmd.py` or `rodeo.py`.** Confirmed by inspection (`git log -p`) that
  `trotter_odmd.py` already uses `reach[0]` (chem-3ds's note: "already correct since `d3d28f0`");
  `rodeo.py`'s `overlap_tol` wiring was chem-3ds's separate fix and doesn't touch `odmd.py` or
  `ODMDProblem` at all.

---

## 6. What could not be verified

- **Whether `sandbox/chem-3ds` will merge before or after this branch** — see §0. This handoff
  flags the conflict; resolving the merge order is a human decision, not something I can verify
  from inside a sandbox cut from `main`.
- **Cross-platform reproducibility of the sq-H4(1.19) witness.** Verified only on this container
  (Linux x86_64, 8 cores — see below). The spec's own history shows this exact geometry family is
  platform/thread-count sensitive (chem-7cx, §0 of `SPEC_eigenstate_reachability.md`); the
  `lib.num_threads(1)` pin inside the test is the same mitigation the repo already uses elsewhere,
  but I have not run this specific test on macOS.
- **Full `make gates`** was not run (§2) — only the directly-affected files.

**Hardware/wall time** (test_odmd_spec.py run, > 1 minute): `Linux 7.0.14-linuxkit x86_64`, 8 CPUs
(`nproc`), 16 GB RAM (`MemTotal` 16353600 kB) — single Docker VM, as described in this session's
compute-envelope guidance. 242s for `test_odmd_spec.py` (5 tests, including the pre-existing
100-seed noise sweeps in G3/G4); 464.64s for the 17-test consumer sweep.

---

## 7. Git status / suggested commands (NOT executed — human review required)

Per this session's git policy, commits happen on this branch (`sandbox/chem-1x3`) but are not
pushed or merged by me. Files touched:

```
 M odmd.py
 M specs/SPEC_odmd.md
 M tests/test_odmd_spec.py
 M .beads/issues.jsonl          (bd update --claim bookkeeping)
?? sandbox-handoffs/chem-1x3.md
```

Suggested commands for a human to review and run (see §0 for the merge-order caveat with
`sandbox/chem-3ds` before merging either branch to `main`):

```bash
git add odmd.py specs/SPEC_odmd.md tests/test_odmd_spec.py .beads/issues.jsonl sandbox-handoffs/chem-1x3.md
git commit -m "chem-1x3: gate odmd.py's ref=reach[0] fix on an HF-unreachable-ground witness

Re-applies chem-3ds's reach[0]-vs-w_eig[0] fix (not yet merged from sandbox/chem-3ds, see
handoff §0) and adds the gate that makes it a verified claim: square H4 at 1.19 Angstrom with
a symmetric SCF has HF population ~0 on the global FCI ground, and ODMDProblem.ref must track
(and does track) the HF-reachable level the survival-amplitude signal can actually converge to.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
bd close chem-1x3 --reason="Added SPEC_odmd.md G5 + tests/test_odmd_spec.py::test_G5_ref_is_hf_reachable_not_global_fci_ground on a sq-H4(1.19)/TIGHT_SCF_CONV_TOL witness where HF population on the global FCI ground is numerically zero. Confirmed RED on the old w_eig[0] formula, GREEN after reach[0]; re-applied chem-3ds's fix since sandbox/chem-3ds is not yet merged to main (handoff \$0 flags the pending merge-order conflict). All 5 test_odmd_spec.py gates + 17 downstream-consumer gates (test_odmd_excited/adaptive_shots/odmd_uq/device_odmd_spec.py) green; ruff clean."
```
