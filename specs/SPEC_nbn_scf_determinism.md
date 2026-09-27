# SPEC: NbN SCF regeneration is deterministic regardless of ambient BLAS/OMP thread count

**Status:** IMPLEMENTED — gate green (`tests/test_nbn_scf_determinism_spec.py`).

---

## 1. Goal

**Claim:** `benchmark_nbn.ground_state_mf` (via `_tight_scf`'s Newton/SOSCF stability-following
loop) returns a **bit-identical** SCF energy for a fixed geometry, regardless of the ambient
`pyscf.lib.num_threads()` setting the caller leaves in place. Before the fix below, this was
false: chem-owr (`sandbox-handoffs/chem-owr.md`) measured 3 of 4 default-thread trials on the real
NbN geometry (`specs/nbn_mp-2634.cif`) landing on one locally-stable UHF minimum and 1 of 4 landing
on a different one, 0.548 mHa apart in SCF energy (4.84 mHa apart after CASCI(14,14)/FCI(10,4)) —
a coin flip driven by floating-point non-associativity in multi-threaded BLAS reductions inside
pyscf's Davidson-based `mf.stability()` call, because the real geometry's (10,4)/S=3 solution sits
essentially exactly on a stability-Hessian zero crossing. The claim is false if two fresh runs at
different ambient thread counts disagree.

## 2. Background and honest framing

- Filed as part of **chem-y6w** (P2): the committed `specs/nbn_scf_reference.chk` /
  `data/nbn_scf.chk` is anchored to one of two genuinely, separately-stable UHF minima in the
  (10,4)/S=3 sector, and which one a *fresh* SCF regeneration finds used to be non-deterministic.
  chem-y6w's decision (`SPEC_nbn_low_spin.md` §0) is to **keep the vendored chk as the reference**
  (re-deriving every downstream number in `SPEC_nbn_dmrg_reference.md` G1-G3b and
  `SPEC_nbn_low_spin.md` Tables 1-2/G4/G5 on the other minimum is out of scope for a P2 bug fix —
  filed separately, `chem-*` follow-up) — but a reference whose *regeneration path* is a coin flip
  is worse than one that is at least honestly, deterministically something else.
- **What you can claim if the gate passes:** any future call to `ground_state_mf`/`_tight_scf` —
  whether from `benchmark_nbn.py`'s own driver, `nbn_low_spin.py`'s `cif`/`scan` subcommands, or a
  fresh clone with no cached `data/nbn_scf.chk` at all — gives the same answer every time, on this
  host, independent of how many threads pyscf happens to be configured for when it's called.
- **What you cannot claim:** that the deterministic answer matches the vendored
  `specs/nbn_scf_reference.chk`. It does not (chem-owr measured it 0.548 mHa higher in SCF energy)
  — pinning determinism does not resolve *which* minimum is "the" reference, it only makes the
  regeneration path stop depending on scheduler/BLAS noise. Nor does this claim the fix generalizes
  to a different host's BLAS backend/build (only tested on this container's OpenBLAS) — see R1.

## 3. Approach

`pyscf.lib.num_threads(n)` is pyscf's own recommended runtime lever for its OpenMP-threaded C
extensions (its docstring explicitly prefers it over `os.environ['OMP_NUM_THREADS']`, which pyscf
reads once at import and cannot be changed after). `mf.stability()`'s Davidson iterations run
through these same pyscf-threaded contraction routines, so pinning `lib.num_threads(1)` for the
duration of `_tight_scf` (restoring the ambient count on exit, so callers elsewhere in the same
process are unaffected) removes the non-associativity. **Reference for this claim:** repeated,
varied-ambient-thread reruns of the exact function under test, each compared bit-for-bit against
the others — no physics reference is needed here (this is a determinism claim about the code, not
an energy claim about the molecule), so the FCI/DMRG ground-truth requirement in `CLAUDE.md`/
`AGENTS.md` does not apply; the "reference" is reproducibility itself, checked directly.

## 4. Public interface

No new public interface. `benchmark_nbn._tight_scf` now pins/restores `pyscf.lib.num_threads`
internally; callers (`ground_state_mf`, `nbn_low_spin.py`) are unchanged.

## 5. Acceptance criteria (validation gates)

`tests/test_nbn_scf_determinism_spec.py` — pyscf only (imports `benchmark_nbn`, which imports
`qiskit`/`qiskit_nature`; **must not** share a process with block2/DMRG tests — own file, per the
`test_*_spec.py` isolation convention).

- **G1 — DEFINITION OF DONE: cross-ambient-thread determinism.** Two fresh `ground_state_mf` runs
  on `specs/nbn_mp-2634.cif` (distinct chkfile paths, so neither restores the other's cache), with
  `pyscf.lib.num_threads()` left at 8 before one call and 2 before the other, give **bit-identical**
  `e_tot` (`==`, not `pytest.approx`). This is exactly the two ambient settings chem-owr's own
  4-default-thread-trial experiment varied implicitly (OS scheduling); pinning explicit values here
  makes the check itself deterministic rather than relying on getting unlucky in CI.
- **G2 — repeat-call determinism at one ambient setting (regression floor).** Three fresh runs, all
  at the default ambient thread count, are pairwise bit-identical. Directly reproduces (in
  miniature, and now made to pass) chem-owr's 4-trial experiment that originally found 3/4 vs 1/4.

## 6. Implementation plan (test-first)

1. `tests/test_nbn_scf_determinism_spec.py` encoding G1/G2 — RED before the fix (verified manually
   by reverting the `lib.num_threads` pin and rerunning; not committed in a permanently-flaky form
   since G1/G2 must be stable CI gates, not coin flips — see §8 R2).
2. `benchmark_nbn._tight_scf`: pin `lib.num_threads(1)` around the SCF+Newton+stability loop,
   restore the ambient count in a `finally`.
3. Green (below).

## 7. Out of scope

- Re-vendoring `specs/nbn_scf_reference.chk` from the deterministic (lower-energy) minimum, and
  re-deriving `SPEC_nbn_dmrg_reference.md` G1-G3b / `SPEC_nbn_low_spin.md` Tables 1-2/G4/G5 against
  it. chem-y6w's decision was to keep the current chk with an explicit caveat instead (see
  `SPEC_nbn_low_spin.md` §0); **chem-czw (2026-09-27) did the re-derivation and confirmed the
  decision**: the other minimum's orbitals give a strictly worse (higher) CASCI(14,14) energy in
  every sector (`SPEC_nbn_low_spin.md` §0b, Table 4) — the vendored chk stays. This spec's own
  scope remains just the *regeneration path's* determinism, unchanged.
- Explaining *why* the real geometry's (10,4)/S=3 UHF solution sits so close to a stability
  boundary in the first place (a near-dissociated 3.7255 Å Nb-N separation vs MP's ~2.25 Å
  reported equilibrium is the leading candidate, per `SPEC_nbn_low_spin.md` §0 — not investigated
  further here).
- Cross-BLAS-backend / cross-host verification (see R1).

## 8. Caveats and risks

- **R1:** verified only on this container (x86_64, 8 cores, OpenBLAS via `scipy-openblas64`,
  pyscf's own `_np_helper` OMP threading). `pyscf.lib.num_threads` controls pyscf's own compiled
  extensions specifically, which is what `mf.stability()`'s contractions run through on this
  build; a pyscf build linked differently, or a stability solver path that instead bottlenecks on
  raw numpy/BLAS matrix multiplies outside pyscf's OMP control, could still be non-deterministic
  there. Not tested on a second machine.
- **R2:** G1/G2 must never be allowed to regress into flaky gates. If a future change reintroduces
  ambient-thread sensitivity, these should FAIL deterministically (not flip pass/fail run to run) —
  if a rerun ever shows either gate flip, that is itself evidence the fix has regressed, not a
  reason to retry and move on.
- ~~A fresh regeneration will now reliably land on the *other* minimum every time.~~ **Platform-
  dependent (2026-09-27, macOS arm64 re-verification, chem-czw).** Deterministic *within* one BLAS
  build, but which minimum the build finds is not portable: on macOS the pinned path lands on
  e_tot = −110.02758538272664 Ha, the **vendored** minimum (−110.0275853827266), not the Linux
  container's −110.02813374576796. G1/G2 still pass there (2/2) because they assert only
  within-build repeatability.
- **R3 — loop-body bug (fixed 2026-09-27).** chem-y6w's rewrite of `_tight_scf` dedented
  `mf.kernel(dm0=dm)` out of the stability loop: no re-convergence after following an instability,
  and an `UnboundLocalError` whenever the first stability check passes (always, on macOS). The
  Linux lower-minimum results (this spec's G1/G2 pins, `scripts/gen_nbn_lower_minimum_chk.py`,
  chem-czw's FCI/DMRG table) were produced by the buggy loop and had not been regenerated with the
  fixed one at the time R3 was first written. **Confirmed 2026-09-27 (chem-ivu, same Linux
  container, host `796434596efe`):** re-ran `scripts/gen_nbn_lower_minimum_chk.py` with the fixed
  loop — lands on `e_tot = -110.0281337457679` (diff 5.68e-14 Ha from the pinned
  `-110.02813374576796`, i.e. the same minimum, not a new one) — then re-ran the cheap (10,4)/S=3
  exact FCI on the freshly regenerated chk and got `-110.04119080464493`, bit-identical to the
  literal chem-czw pinned from the buggy-loop run
  (`tests/test_nbn_czw_lower_minimum_spec.py::LOWER_MINIMUM_FCI[(10, 4)]`). The buggy loop
  happened not to have corrupted these Linux numbers: `mf.stability()` found this minimum
  internally stable on its first check, so the dedented `mf.kernel(dm0=dm)` (only reachable after
  an *unstable* verdict) was simply never executed on this path — the bug was latent here, not
  triggered. No re-derivation needed; the S=1/S=2/(7,7)-DMRG rows in the same table were not
  independently re-run (out of scope for this confirmation — see
  `sandbox-handoffs/chem-ivu.md`).

## 9. Deliverables

- `benchmark_nbn.py` — `_tight_scf` pins/restores `pyscf.lib.num_threads`.
- `tests/test_nbn_scf_determinism_spec.py` — G1, G2.
- `specs/SPEC_nbn_low_spin.md` §0 — chem-y6w's explicit decision (keep the vendored chk, caveat the
  dependents, fix the regeneration path's determinism) recorded there, not just here.
