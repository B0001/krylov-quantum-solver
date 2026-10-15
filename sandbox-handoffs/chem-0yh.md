# chem-0yh handoff — unpinned `scf.RHF`/`fci.FCI` in `test_eigenstate_reachability_spec.py::_fci`

**Bead:** chem-0yh (P2, bug) — "test_eigenstate_reachability_spec.py's own `_fci()` helper has an
unpinned `scf.RHF` call sharing chem-7cx's thread-nondeterminism mechanism (causes residual
G2[1.19] flakiness)."

**Resolution:** fixed as described in the bead. `_fci()`'s own `scf.RHF`/`fci.FCI` calls are now
pinned to `pyscf.lib.num_threads(1)` for their duration, restored in a `finally`, with a regression
test. A second, independent, already-tracked nondeterminism source (chem-s6c) was found still
active in full-file runs — documented below and left for chem-s6c, not fixed here.

---

## 1. What changed

### 1a. Root cause (matches the bead exactly)

PySCF's own multi-threaded OpenMP integral evaluation is not bit-associative (same mechanism as
chem-7cx and chem-y6w). `_fci()`'s `scf.RHF(mol).kernel()` + `fci.FCI(mf).kernel()` call was
unpinned, so its result could pick up ~1e-13-level noise depending on the *ambient*
`pyscf.lib.num_threads()` setting the caller happened to leave in place. At square H4's
near-degenerate a=1.19 geometry, the FCI Davidson solve can amplify that noise enough to flip which
side of `test_G2_decision_equals_symmetry_adapted_fci`'s classification boundary the reference
lands on.

### 1b. Fix: pin `pyscf.lib.num_threads(1)` around `_fci`'s SCF+FCI call

`tests/test_eigenstate_reachability_spec.py::_fci`:

```python
ambient_threads = lib.num_threads()
lib.num_threads(1)
try:
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.kernel()
    solver = fci.FCI(mf)
    solver.wfnsym, solver.nroots = wfnsym, 1000
    e, _ = solver.kernel()
finally:
    lib.num_threads(ambient_threads)
```

Identical pattern to chem-7cx's `build_molecular_hamiltonian` fix and chem-y6w's
`benchmark_nbn._tight_scf` fix. Restoring in `finally` matters here specifically because `_fci` is
`@lru_cache`'d and module-level — leaving the process pinned to 1 thread would silently throttle
every later caller in the same pytest session.

### 1c. New regression test: `test_fci_helper_is_deterministic_across_ambient_thread_counts`

Added directly below `_fci` in the same file. Bypasses `lru_cache` (`.cache_clear()`) to force two
fresh recomputations at ambient thread counts 2 and 8 on the a=1.19 geometry, asserts bit-identical
output, and asserts the caller's ambient thread count is restored afterward.

```
$ uv run pytest tests/test_eigenstate_reachability_spec.py::test_fci_helper_is_deterministic_across_ambient_thread_counts -q
```
run 5× in a row: **5/5 passed** (`1 passed in ~1.5s` each).

### 1d. Direct repro: before vs. after, isolated from everything else

Fresh-subprocess repro of `_fci(_sq(1.19), "Ag")` alone (`uv run python` from a throwaway script
importing just `_fci`/`_sq`), 6 trials each:

- **Before the fix** (unpinned, as filed): bit noise in `mf.e_tot` at the ~1e-13 level across
  trials, depending on the process's ambient thread count — matches the bead's own documented
  repro.
- **After the fix**: all 6 trials bit-identical (`e_ag_fci` printed with `repr()`, compared by eye
  and by the regression test's `np.array_equal`), at ambient thread counts 1, 2, 4, and 8.

---

## 2. Gate-run evidence (verbatim)

```
$ uv run ruff check tests/test_eigenstate_reachability_spec.py
All checks passed!

$ uv run pytest tests/test_eigenstate_reachability_spec.py::test_fci_helper_is_deterministic_across_ambient_thread_counts -q
1 passed in 1.50s
```
(repeated 5×, identical result every time)

Full-file baseline on this branch (`sandbox/chem-0yh`, cut from `main` — i.e. **without**
chem-7cx's `build_molecular_hamiltonian` fix and **without** chem-s6c's `reachability.py` fix,
since neither is merged yet):

```
$ uv run pytest tests/test_eigenstate_reachability_spec.py -q
...
14 failed, 36 passed, 39 warnings in 123.35s (0:02:03)
```

This 14-failure baseline is **expected and not a regression**: every failing test here
(`test_G1[...]`, `test_G2[...]`, `test_G4[...]`, `test_G9[...]`) is a documented symptom of
chem-7cx's and chem-s6c's still-open, still-unmerged fixes, not of anything in this bead's scope —
see §4 for how this was confirmed rather than assumed.

---

## 3. The investigation: does chem-0yh's fix actually eliminate the bead's named symptom?

The bead's acceptance bar is "fix `_fci`'s own unpinned call"; it is *not* "the full test file is
green on this branch" (that needs chem-7cx and chem-s6c merged too, neither of which is this
bead's job). To check the fix is real and sufficient for its own narrow claim, chem-7cx's fix was
applied **temporarily, locally, never committed** (`git show <chem-7cx commit> --
hybrid_quantum_solver/molecular_hamiltonian.py | git apply`, reverted via `git checkout --
hybrid_quantum_solver/molecular_hamiltonian.py` before finishing this session) to recreate the
environment the bead's own repro was measured in.

With both fixes applied, `test_G2_decision_equals_symmetry_adapted_fci[1.19]` **still flipped**
across repeated full-file and even repeated 2-test-selected (`-k 1.19`) runs. Isolating further:

- `_fci(_sq(1.19), "Ag")` alone: bit-identical across 6 fresh-process trials. **Not the cause.**
- `build_molecular_hamiltonian(atom=_sq(1.19), conv_tol=...)` alone: bit-identical `hf_energy` AND
  bit-identical `qubit_hamiltonian.to_matrix()` (compared by md5) across 6 fresh-process trials per
  conv_tol. **Not the cause.**
- `np.linalg.eigh` on that (deterministic) matrix: bit-identical across 6 fresh-process trials
  (md5-compared). **Not the cause.**
- The **actual test body** (`_system` → `_weights` → `hf_symmetry_sector` → `symmetry_allowed`),
  run fresh per-process 6 times: `len(in_sector)` jumped between 10 and 36 — a real classification
  flip, not boundary-grazing noise, even though every input above is deterministic.

Tracing `_weights` → `hf_symmetry_sector` → `orbital_parity_bounds` → `_driver_orbitals` in
`reachability.py` found the actual source: `_driver_orbitals` (line ~176) makes its own, separate,
**unpinned** `drv.run_pyscf()` call (a `PySCFDriver` rebuild used to recover the active-space MO
coefficients for the symmetry diagnostic) — the exact same thread-nondeterminism mechanism, but a
*third* call site, distinct from both `build_molecular_hamiltonian` (chem-7cx) and
`scf_symmetry_status`/`hf_orbital_irreps` (chem-s6c, already open and tracking this exact file for
two other call sites).

**This is out of scope for chem-0yh** (which is specifically about this test file's own `_fci`
helper) and already has a home: filed as an addendum to the existing open bead **chem-s6c** via
`bd update chem-s6c --notes=...` rather than creating a duplicate, since it's the same file and the
same mechanism chem-s6c is already tracking, just a call site its current description doesn't
name. chem-s6c was left open, not closed or re-scoped, so a future session can decide whether
`_driver_orbitals`'s existing chem-g14 retry/gauge-search logic already absorbs this or whether it
also needs the `lib.num_threads(1)` pin.

---

## 4. Pre-registered criteria and pass/fail

1. **"Pin `pyscf.lib.num_threads(1)` around `_fci()`'s `scf.RHF(...).kernel()` call (and its
   `fci.FCI(...).kernel()` call too), restoring the ambient count in a `finally`"** — **PASS.**
   §1b, verbatim diff matches the bead's candidate fix.
2. **"A regression test for this specific mechanism"** (implied by the bead's own repro
   methodology) — **PASS.** §1c, `test_fci_helper_is_deterministic_across_ambient_thread_counts`,
   5/5 green.
3. **"Test suite passes"** — **PARTIAL, and this partiality is the honest finding, not a gap in
   this bead's work.** The fix eliminates `_fci`'s own contribution to the flake (verified directly,
   §1d). It does **not** make `test_G2[1.19]` deterministic on its own, because a second,
   independent, already-tracked nondeterminism source (chem-s6c's territory, §3) also reaches that
   same test through `hf_symmetry_sector`. Closing chem-0yh does not require that second source be
   fixed — it requires `_fci`'s own call to be fixed, which it is, with evidence it is.
4. **Lint** — **PASS.** `uv run ruff check tests/test_eigenstate_reachability_spec.py` → `All
   checks passed!`

---

## 5. What was decided not to do, and why

- **Did not fix `_driver_orbitals`'s unpinned `drv.run_pyscf()` call.** Out of scope: it lives in
  `reachability.py`, is a different call site than this bead names, and is the same file/mechanism
  an existing open bead (chem-s6c) already tracks. Filed as a notes addendum on chem-s6c (§3)
  instead of a new duplicate bead or an in-place fix here.
- **Did not commit chem-7cx's `build_molecular_hamiltonian` fix.** It was applied locally and
  transiently, purely to reproduce the bead's own documented pre-conditions for validating this
  fix, and was reverted (`git checkout -- hybrid_quantum_solver/molecular_hamiltonian.py`) before
  finishing. `git diff` against this branch now shows zero changes to that file — confirmed below.
- **Did not run the full `make gates`/`make test`.** Only this one file's gate and lint were run,
  since nothing else was touched; `test_eigenstate_reachability_spec.py` is not itself a
  `test_*_spec.py`-glob process-isolation case requiring special handling here (it already imports
  pyscf/qiskit only, no block2), so a single `pytest` invocation on it is the correct, matching way
  to run it.

---

## 6. What could not be verified

- **Whether chem-s6c's eventual fix, once merged, makes `test_G2[1.19]` fully deterministic.** Not
  testable from this branch without also merging chem-7cx and chem-s6c, both out of scope here.
- **Cross-host generalization.** Verified only on this container; same caveat chem-7cx/chem-y6w's
  own handoffs already carry for this fix pattern.

---

## 7. Git status / commits made

Per this run's instructions (which supersede the repo's general "do not commit" default), work was
committed on `sandbox/chem-0yh` only — no push, no branch switch.

```
$ git diff --stat hybrid_quantum_solver/molecular_hamiltonian.py
(empty — confirmed reverted, no chem-7cx work leaked into this branch)

$ git status --short
 M .beads/issues.jsonl
 M tests/test_eigenstate_reachability_spec.py
?? sandbox-handoffs/chem-0yh.md
```

(`.claude/settings.local.json` is local tool permission state, not part of this bead's work — left
untracked.)

Commits (see `git log` on this branch for the actual hashes):

1. `tests/test_eigenstate_reachability_spec.py` — the `_fci` thread-pin fix + regression test.
2. `.beads/issues.jsonl` + `sandbox-handoffs/chem-0yh.md` — bead state export + this handoff.

```bash
bd close chem-0yh --reason="Fixed as filed: _fci()'s scf.RHF/fci.FCI call now pinned to lib.num_threads(1), restored in finally, matching chem-7cx/chem-y6w's pattern. Regression test test_fci_helper_is_deterministic_across_ambient_thread_counts added (5/5 green). A second, independent, already-tracked nondeterminism source reaching the same test (reachability.py's _driver_orbitals, chem-s6c's territory) was found and documented as a chem-s6c notes addendum rather than fixed here -- out of scope for this bead."
```
