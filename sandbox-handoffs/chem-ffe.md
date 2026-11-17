# chem-ffe — Linux verification pass for the 2026-10-02 batch PRs (#58, #59)

## Scope (from the bead description, verbatim)

1. chem-8tl (PR #58): the revised `centered_pds_roots` G1/G1b gates (worst pre-registered
   summation/term-order path) were only run on macOS/Accelerate; §6's Linux table came from the
   older, pre-revision authoring run (commit 1ff123a).
2. `tests/test_reachability_tolerance_spec.py` still built its Hamiltonian with a bare
   `mapper.map(...)` instead of `map_canonical` (hash-seed dependent), which PR #58 made the
   norm elsewhere.

Both are done. Nothing else was in scope — see "What I deliberately did not do" below.

## Hardware / environment

`uname -a`: `Linux fd56fe9691fe 7.0.14-linuxkit #1 SMP PREEMPT Fri Sep 18 10:19:32 UTC 2026 x86_64
GNU/Linux`. 8 vCPU (`nproc`), one worker. Python 3.13.15, numpy 2.4.6, scipy 1.15.3
(scipy-openblas 0.3.31.188.0, `DYNAMIC_ARCH` Haswell, `USE64BITINT`), qiskit 2.5.0, qiskit-nature
0.8.0, pyscf 2.13.1 — same package versions as the macOS run `SPEC_centered_pds_roots.md` already
had, different BLAS build/container.

## 1. `centered_pds_roots` — Linux re-run of the revised gates

```
uv run pytest tests/test_centered_pds_roots_spec.py -v
```
→ **4 passed in 37.96s** (G1, G1b, G2, G3). Re-ran twice more, same result.

Captured the pinned-configuration JSON directly:
```
PYTHONHASHSEED=2 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 uv run python tests/test_centered_pds_roots_spec.py
```
run twice → byte-identical output (G1's own fresh-process check). Ran again at
`PYTHONHASHSEED=0` → byte-identical to `PYTHONHASHSEED=2` (G1b). This is the data vendored into
`specs/SPEC_centered_pds_roots.md` §6's new "Revised (chem-8tl), Linux x86_64 / OpenBLAS ... (chem-ffe,
2026-10-15)" table and surrounding prose; the old 17-path, commit-1ff123a table is kept underneath,
explicitly marked superseded, for the record.

**The finding that matters beyond "gates pass":** the worst of the 33 bit-paths is platform- *and*
container-specific, more sharply than the spec previously said. Three runs, three different
answers, same package versions where shared:

| run | worst path | raw K=7 margin |
|---|---|---|
| chem-8tl authoring run, Linux/OpenBLAS (commit 1ff123a, 17 paths) | shuffle 15 | −5.9075 mHa |
| chem-8tl, macOS 27/Accelerate (33 paths) | shuffle 7 | −11.9992 mHa |
| **chem-ffe, this Linux/OpenBLAS container (33 paths)** | **shuffle 14** | **−9.8664 mHa** |

Shuffle 15 (the old Linux pin) and shuffle 7 (the macOS pin) both *respect* the bound here
(+0.1838 mHa and +0.0607 mHa). G1 already took "worst of the pre-registered set" rather than a
named path for exactly this reason (R1); this run is further evidence for it, now recorded in R1's
text. The claim G1–G3 actually test (centering removes the sensitivity) held on all three runs: the
centered margins agree to 6.1e-14 Ha here, well inside the 1e-8 Ha gate, and every centered root
stayed ≥ the reachable target.

Updated: `specs/SPEC_centered_pds_roots.md` — top status line, the G1 bullet in §5, a new §6 Linux
table (chem-ffe) with the superseded chem-8tl-authoring table kept below it, and R1 in §8.

## 2. `test_reachability_tolerance_spec.py` — `map_canonical` fix

`_p0_at_conv_tol` (the helper behind gate G6) built its qubit Hamiltonian with
`JordanWignerMapper().map(prob.hamiltonian.second_q_op())`. Per chem-8tl's own finding (now in
`molecular_hamiltonian.map_canonical`'s docstring), that is `PYTHONHASHSEED`-dependent in the low
coefficient bits and the term order. Changed it to `map_canonical(mapper, ...)`, matching
`build_molecular_hamiltonian` and every other site chem-8tl touched. Diff:

```python
-    H = mapper.map(prob.hamiltonian.second_q_op()).to_matrix()
+    H = map_canonical(mapper, prob.hamiltonian.second_q_op()).to_matrix()
```

Grepped the rest of PR #59's touched test files (`test_scf_conv_tol_spec.py`,
`test_chemcheck_spec.py`, `test_chained_overlap_spec.py`) for the same pattern — none found, so
this was the only site named by the bead and the only one needing it in that PR's diff.

**Verified the fix changes nothing it shouldn't:** ran the file's full gate set before and after
the edit (`git stash` / `git stash pop` to bracket it). Identical pass/fail set both times (6
failed, 6 passed) — `test_G1`–`test_G4[1.2/1.3]`, `test_G6` fail; `test_G4[1.0/1.3→passed cases]`,
`test_G5*`, `test_G7`, `test_G8` pass. Those 6 failures are the RHF branch-selection divergence
chem-8ls is already tracking for this exact file (its description lists
`test_reachability_tolerance_spec.py` G1/G2/G3/G4/G5-or-G6 by name) — not caused or fixed by this
change, and not touched here. Added a short maintenance note to
`specs/SPEC_reachability_tolerance.md` §2b pointing at `map_canonical` and at chem-8ls for the
unrelated failures, so a future reader doesn't reopen either investigation.

## Gate-run lines (verbatim pass/fail counts)

```
$ uv run ruff check .
All checks passed!

$ uv run pytest tests/test_centered_pds_roots_spec.py -v
4 passed in 37.96s

$ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  uv run pytest tests/test_reachability_tolerance_spec.py -v
6 failed, 6 passed in 5.03s     # identical before/after the map_canonical edit (chem-8ls territory)

$ uv run pytest tests/test_reference_energies.py tests/test_krylov_convergence.py \
  tests/test_noise_resilience.py tests/test_hamiltonian_determinism.py -q
21 passed in 41.33s

$ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  uv run pytest tests/test_centered_pds_roots_spec.py tests/test_hamiltonian_determinism.py \
  tests/test_reference_energies.py tests/test_krylov_convergence.py tests/test_noise_resilience.py -v
25 passed in 88.96s
```

Did not run the full 114-file `make gates` sweep. This bead's two action items are both covered by
the files above; a full sweep's remaining Linux-divergence findings (square-H4 witness families in
`test_scf_conv_tol_spec.py`, `test_eigenstate_reachability_spec.py`,
`test_symmetry_reachability_spec.py`, and this file) are chem-8ls's explicit, already-in-progress
scope (owned by "sandbox", status `in_progress` as of this session) — reproduced one of them
(`test_scf_conv_tol_spec.py`: 4 failed / 10 passed, `H4 square 1.35`/`1.19` witnesses land on the
broken-symmetry branch here) only to confirm it's the same family chem-8ls already describes, then
left it alone rather than duplicate that investigation.

## Pre-registered criteria vs. outcome

| Criterion (from the bead) | Outcome |
|---|---|
| Run the revised G1/G1b gates on Linux, record §6 | Done — 4/4 pass, bit-reproducible, data vendored in the spec |
| Fix `test_reachability_tolerance_spec.py`'s un-canonicalized map | Done — now uses `map_canonical`; verified no change to pass/fail set |

## What I deliberately did not do, and why

- **Did not touch the square-H4 witness non-reproducibility** (`test_scf_conv_tol_spec.py`,
  `test_eigenstate_reachability_spec.py`, `test_symmetry_reachability_spec.py`, the other failures
  in `test_reachability_tolerance_spec.py` itself). chem-8ls already owns this exact family,
  in_progress, with its own repro and hardware notes — duplicating it here would fork the
  investigation.
- **Did not run the full `make gates` (114 spec files).** The bead names two specific, narrow
  fixes; both are verified in isolation plus the core validation suite. A full sweep would mostly
  re-surface chem-8ls's already-filed findings at a much higher time/compute cost on an 8 vCPU / 16
  GB shared box.
- **Did not change `reachability.py`'s own bare `mapper.map(...)` call** (`_driver_orbitals`,
  line ~211) even though it's the same anti-pattern. Checked: that site's divergence (~4e-8 per
  coefficient, per chem-g14's investigation) is loose-SCF-convergence noise, not hash-seed/term-order
  bit noise (~1e-16) — `map_canonical` would not fix chem-g14's issue, and chem-g14 is already open
  and tracking it. Left as-is; not filing a new bead since it's not a live problem map_canonical
  would solve.

## What I could not verify

- Whether `map_canonical` changes any bit of the numbers already recorded in
  `SPEC_reachability_tolerance.md` §2b's conv_tol table — could not regenerate that table
  end-to-end on this container because `test_G6` already fails here for the unrelated chem-8ls
  branch-selection reason (same geometry). Argued from the physics (mapping order can't change
  which RHF branch PySCF converges to, and this Hamiltonian is nowhere near PDS's cond~1e17
  sensitivity regime) rather than measuring it directly — flagged as such in the spec's new note.

## Files changed

- `specs/SPEC_centered_pds_roots.md` — Linux re-run data (§6), status line, G1 bullet, R1.
- `specs/SPEC_reachability_tolerance.md` — maintenance note on the `map_canonical` fix.
- `tests/test_reachability_tolerance_spec.py` — `_p0_at_conv_tol` now uses `map_canonical`.
- `.beads/issues.jsonl` — bead state (claim).

## Proposed commands (not run — git policy for this session is commit-only, no push)

```
git add specs/SPEC_centered_pds_roots.md specs/SPEC_reachability_tolerance.md \
        tests/test_reachability_tolerance_spec.py .beads/issues.jsonl sandbox-handoffs/chem-ffe.md
git commit -m "chem-ffe: Linux-verify centered_pds_roots revised gates; map_canonical in reachability_tolerance G6"
```

`bd close chem-ffe` after the above, once reviewed.
