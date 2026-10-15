# SPEC: build_molecular_hamiltonian's free RHF is deterministic regardless of ambient BLAS/OMP thread count

**Status:** IMPLEMENTED — gate green (`tests/test_scf_branch_determinism_spec.py`).

---

## 1. Goal

**Claim:** `hybrid_quantum_solver.molecular_hamiltonian.build_molecular_hamiltonian` returns a
**bit-identical** `MolecularHamiltonian` (same `hf_energy`, same qubit Hamiltonian) for a fixed
`(atom, basis, charge, spin, conv_tol, active_electrons, active_orbitals)`, regardless of the
ambient `pyscf.lib.num_threads()` setting the caller leaves in place or of repeated calls at a
fixed setting. Before the fix below, this was false on square H4 (sto-3g, a near-degenerate,
symmetry-rich geometry): chem-7cx measured 9 of 10 default-thread trials at a = 1.1 landing on
`hf_energy = -1.7825510081...` and 1 of 10 on a different locally-stable RHF solution,
`hf_energy = -1.7825511822...`, differing in the 7th digit — a coin flip driven by floating-point
non-associativity in PySCF's multi-threaded integral evaluation, exactly the same mechanism
`specs/SPEC_nbn_scf_determinism.md` (chem-y6w) found and fixed for `benchmark_nbn._tight_scf`. The
claim is false if two fresh builds at different ambient thread counts, or two fresh builds at the
same ambient count, disagree.

## 2. Background and honest framing

- Filed as **chem-7cx** (P2), found while fixing chem-g14 (`reachability._driver_orbitals`
  rebuild-reproducibility — a *different* mechanism, the rebuild's own operator-equality check,
  not the original SCF). chem-7cx's own repro: 10 fresh `build_molecular_hamiltonian` calls on
  square H4 (`"H 0 0 0; H 1.1 0 0; H 1.1 1.1 0; H 0 1.1 0"`, `conv_tol=1e-6`) under default
  (multi-threaded) ambient settings. Because `tests/test_eigenstate_reachability_spec.py`'s own
  `_system()` caches one build per `(atom, conv_tol)` per pytest process, the flip was invisible
  within a single test run but visible across repeated full-file invocations: `test_G1` and
  `test_G9_control` pass or fail depending on which branch that process's one build happened to
  land on — never both, never predictably.
- **Root cause, isolated directly (not inferred):** pinning only `OMP_NUM_THREADS=1` (PySCF's own
  OpenMP-threaded integral evaluation, read once at interpreter start) made 10/10 repeated builds
  bit-identical; pinning only `OPENBLAS_NUM_THREADS=1` (numpy/scipy's BLAS, used by the SCF's own
  linear algebra) did **not** — the builds still varied. So the noise enters through PySCF's
  compiled integral/Fock-build extensions, the same subsystem `SPEC_nbn_scf_determinism.md` found
  for NbN's `mf.stability()` Davidson contractions, not through numpy's linear algebra.
  `pyscf.lib.num_threads(n)` is the runtime lever that actually controls that subsystem (its own
  docstring prefers it over the `OMP_NUM_THREADS` env var, which PySCF reads only once at import
  and cannot change afterward) — confirmed: wrapping just the `PySCFDriver(...).run()` call in
  `lib.with_omp_threads(1)` reproduces the `OMP_NUM_THREADS=1` result exactly, with the ambient
  setting left untouched for every other caller in the process.
- **What you can claim if the gate passes:** any future call to `build_molecular_hamiltonian` for
  a fixed geometry/argument set returns the same answer every time, on this host, independent of
  how many threads PySCF happens to be configured for when it is called. The flakiness in
  `test_G1`/`test_G9_control` across repeated full-suite runs stops being a coin flip.
- **What you cannot claim:** that the now-deterministic branch is the "physically correct" one, or
  that it matches `specs/SPEC_eigenstate_reachability.md`'s macOS-measured witness values. It does
  not always: `test_G9_control_tight_scf_stays_on_the_ag_level[1.1]` and three `test_G1` cases
  still fail, deterministically, on this host/pyscf version after this fix — the SAME four cases,
  every run (5/5 verified, §9). That is chem-g14's territory (a pre-existing, already-filed,
  separate cross-platform/version divergence in which RHF solution this host's PySCF reaches at
  all, not a reproducibility defect) — this spec's claim is narrowly about reproducibility, not
  about matching another platform's physics. Nor does this fix touch `reachability.py`'s own raw
  `scf.RHF(...)` calls (`scf_symmetry_status`, `hf_orbital_irreps`, used by `test_G4` and others) —
  those run outside `build_molecular_hamiltonian` entirely and are a separate, filed follow-up
  (chem-7cx's notes document the same mechanism there; out of scope here, see §7).

## 3. Approach

Pin `pyscf.lib.num_threads(1)` around the `PySCFDriver(...).run()` call inside
`build_molecular_hamiltonian`, restoring the ambient count in a `finally` so callers elsewhere in
the same process are unaffected — the identical pattern `benchmark_nbn._tight_scf` already uses
(chem-y6w). **Reference for this claim:** repeated, varied-ambient-thread reruns of the exact
function under test, compared bit-for-bit against each other — this is a determinism claim about
the code, not an energy claim about the molecule, so no FCI/DMRG ground truth is needed (same
reasoning `SPEC_nbn_scf_determinism.md` §3 used).

## 4. Public interface

No new public interface and no signature change. `build_molecular_hamiltonian` pins/restores
`pyscf.lib.num_threads` internally around its one PySCF driver call.

## 5. Acceptance criteria (validation gates)

`tests/test_scf_branch_determinism_spec.py` — pyscf + qiskit only (imports
`hybrid_quantum_solver.molecular_hamiltonian`, which imports `qiskit`/`qiskit_nature`); **must
not** share a process with block2/DMRG tests — own file, per the `test_*_spec.py` isolation
convention (CLAUDE.md).

- **G1 — DEFINITION OF DONE: cross-ambient-thread determinism.** Two fresh
  `build_molecular_hamiltonian` calls on square H4 (`a = 1.1`, `conv_tol=1e-6`) in separate
  subprocesses, one with ambient `OMP_NUM_THREADS=2` and one with `OMP_NUM_THREADS=8` set before
  the interpreter starts, give **bit-identical** `hf_energy` (`==`, not `pytest.approx`) and
  identical qubit-Hamiltonian digests (SHA1 of sorted Pauli labels + coefficient bytes, the same
  check `tests/test_hamiltonian_determinism.py` uses for the hash-seed claim).
- **G2 — repeat-call determinism at one ambient setting (regression floor).** Ten fresh builds, all
  in separate subprocesses at the same (default) ambient thread count, are pairwise bit-identical.
  Directly reproduces, in full, chem-7cx's own 10-trial experiment that originally found 9/1 — now
  made to pass every time, not most of the time.
- **G3 — same-process repeat-call determinism (the `_system`-cache-masked case).** Ten fresh,
  uncached builds of the same geometry *within a single process* (no subprocess spawn — this is
  what a single pytest process doing ten builds looks like, as opposed to G2's ten separate
  processes) are pairwise bit-identical, and leave the process's ambient
  `pyscf.lib.num_threads()` unchanged after each call (the `finally`-restore is load-bearing: a
  caller elsewhere in the same process must not find its own thread count silently pinned to 1).
- **G4 — the two known-flipping witnesses stabilize, pass or fail.** `test_eigenstate_reachability_
  spec.py::test_G1_residue_level_identified_warned_and_kept[1.1-1e-06-1e-08]` and
  `::test_G9_control_tight_scf_stays_on_the_ag_level[1.1]` are run five times as separate `pytest`
  subprocess invocations; this gate asserts the pass/fail OUTCOME is identical across all five (it
  does not assert which outcome — chem-g14 governs whether it should be GREEN on this host, this
  spec only requires it stop being a coin flip).

## 6. Implementation plan (test-first)

1. `tests/test_scf_branch_determinism_spec.py` encoding G1-G4 — RED before the fix (verified
   manually by reverting the `lib.num_threads` pin and rerunning G2/G3, which showed exactly the
   chem-7cx-reported split; not committed in a permanently-flaky form).
2. `hybrid_quantum_solver/molecular_hamiltonian.py`: pin `lib.num_threads(1)` around the
   `PySCFDriver(...).run()` call in `build_molecular_hamiltonian`, restore the ambient count in a
   `finally`.
3. Green (§9).

## 7. Out of scope

- **chem-g14** — `reachability._driver_orbitals`'s rebuild-reproducibility check and the resulting
  17 pre-existing `test_eigenstate_reachability_spec.py` failures on this host/pyscf version. This
  spec's G4 deliberately does not require those two witnesses to turn green, only to stop flipping.
- **`reachability.py`'s own raw `scf.RHF` calls** (`scf_symmetry_status`, `hf_orbital_irreps`) —
  chem-7cx's notes record the identical multi-threaded-BLAS mechanism there too (e.g. square H4
  a = 1.20's "broken RHF" classification flipping, and an exactly-degenerate orbital pair's parity
  bit flipping at a = 1.00 with no energy change), but those functions build their own `gto.M`/
  `scf.RHF` objects directly, outside `build_molecular_hamiltonian` — fixing them is the same
  three-line pattern but a separate function, a separate blast radius, and arguably a different
  bead (test_G4 and `SPEC_symmetry_reachability`'s per-system availability check live there, not
  here). Filed separately: **chem-s6c**.
- Re-deriving which branch is "correct" relative to `specs/SPEC_eigenstate_reachability.md`'s
  macOS-measured witness numbers — not this spec's claim (see §2).
- Cross-host / cross-BLAS-backend verification (see §8 R1).

## 8. Caveats and risks

- **R1:** verified only on this container (x86_64, 8 cores, OpenBLAS via `scipy-openblas64`,
  PySCF's own `_np_helper`/integral-evaluation OMP threading). A PySCF build linked differently
  could still be non-deterministic through a code path `lib.num_threads` does not control — not
  tested on a second machine (`SPEC_nbn_scf_determinism.md` R1 records the same caveat, and its
  macOS re-verification found the pin IS still deterministic there, just lands on a different
  branch than Linux — platform-dependent WHICH branch, not whether pinning gives determinism).
- **R2:** G1-G4 must never be allowed to regress into flaky gates. If a future change
  reintroduces ambient-thread sensitivity, these should FAIL deterministically (not flip pass/fail
  run to run) — a rerun that shows any of them flip is itself evidence of a regression, not a
  reason to retry and move on (same discipline as `SPEC_nbn_scf_determinism.md` R2).
- **R3:** this does not change `effective_conv_tol` or the `symmetric_scf` opt-in's behavior in any
  way — it only makes whichever branch the existing arguments were already going to select
  reproducible. The blast-radius analysis in `specs/SPEC_symmetric_scf.md` (why `symmetric_scf`
  defaults to `False`) is unaffected: this fix changes zero recorded numbers on any system where
  the free SCF was already stable (i.e. everywhere except near-degenerate, exact-point-group
  geometries like square H4 at these bond lengths) — G2/G3 would have passed on an ordinary system
  even before this fix; they are only interesting, and were only RED, on the witness geometry.

## 9. Results

**Linux x86_64 container, 8 cores, pyscf (via qiskit-nature's `PySCFDriver`), OpenBLAS
(`scipy-openblas64`), 2026-10-14.** `uname -a`: see handoff. Wall time for the full gate file:
recorded in the handoff (`sandbox-handoffs/chem-7cx.md`), well under a minute (ten repeated cheap
sto-3g H4 builds plus five `pytest` subprocess spawns for G4).

- `tests/test_scf_branch_determinism_spec.py -q`: **G1-G4 pass** (see handoff for the verbatim
  line).
- Before the fix (`lib.num_threads` pin reverted, manual check): G2's 10-trial loop reproduced
  chem-7cx's split exactly — 9 of 10 fresh builds at `hf_energy = -1.7825510081366662`, 1 of 10 at
  `-1.7825511822577162` (one single outlier run, consistent with a rare race rather than a 50/50
  coin flip); G4's five-rerun check showed `test_G9_control[1.1]` and `test_G1[1.1-1e-06-1e-08]`
  each flip pass/fail at least once across five reruns (not committed in this flaky form — see
  implementation plan).
- After the fix: G2/G3 are 10/10 bit-identical in every observed run; G4's five reruns give the
  identical FAIL/FAIL outcome for `test_G9_control[1.1]`/`test_G1[1.1-1e-06-1e-08]` every time
  (chem-g14's territory, not this spec's to turn green — see §2, §7).

## 10. Deliverables

- `hybrid_quantum_solver/molecular_hamiltonian.py` — `build_molecular_hamiltonian` pins/restores
  `pyscf.lib.num_threads` internally around its PySCFDriver call.
- `tests/test_scf_branch_determinism_spec.py` — G1-G4.
- `chem-s6c` (new, filed, open) — the same fix for `reachability.py`'s own `scf.RHF` calls.
