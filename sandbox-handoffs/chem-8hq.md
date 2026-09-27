# chem-8hq: excited_bounds.py __main__ Be2 demo: unchecked fci.direct_spin0.kernel(nroots=3)

## What changed

`excited_bounds.py`'s `__main__` block (the Be₂ CAS(4,8)/cc-pVDZ showcase, line ~317) called
`fci.direct_spin0.kernel(h1, eri, 8, ne, ecore=float(e_core), nroots=3)` with no convergence check
on any of the 3 Davidson roots — same unchecked-Davidson bug class chem-3f8 fixed for
`direct_spin1.kernel` call sites, but out of that bead's grep scope (`direct_spin0` + `nroots`
shape).

Replaced the bare `.kernel()` call with an explicit `fci.direct_spin0.FCI()` instance
(`max_cycle=1000`, `conv_tol=1e-10`, `nroots=3`), then `raise RuntimeError` if
`not np.all(fci_solver.converged)` — the same pattern as `nbn_low_spin.py:98-115`'s `spin_fci` and
`hybrid_quantum_solver/dmrg_reference.py`'s `fci_energy`. `fci_energy` itself was not reused because
it does not accept `nroots` (noted in the bead), so this is a local, minimal fix rather than a new
shared wrapper — only one call site in the repo needed it.

Diff (only change made):

```diff
-    roots = fci.direct_spin0.kernel(h1, eri, 8, ne, ecore=float(e_core), nroots=3)[0]
+    fci_solver = fci.direct_spin0.FCI()
+    fci_solver.max_cycle, fci_solver.conv_tol, fci_solver.nroots = 1000, 1e-10, 3
+    roots, _ = fci_solver.kernel(h1, eri, 8, ne, ecore=float(e_core))
+    roots = np.atleast_1d(roots)
+    if not np.all(fci_solver.converged):
+        raise RuntimeError(f"Be2 CASCI FCI did not converge: converged={fci_solver.converged}, "
+                            f"roots={np.round(roots, 6)}")
```

## Verification that the recorded Be2 numbers did not silently change

Before touching the file, I ran the old unchecked call and the new checked call side-by-side
(`/tmp/test_be2_fci.py`, same CAS(4,8)/cc-pVDZ Be₂ integrals at R=2.45 Å):

```
old (unchecked) roots: [-29.18412726 -29.08503265 -29.08503265]
converged flags:       [ True  True  True]
new roots:             [-29.18412726 -29.08503265 -29.08503265]
```

All 3 roots were already converged at the default settings — the fix changes nothing numerically,
it only makes the previously-implicit convergence explicit and would now raise instead of silently
printing a bad number if it ever weren't converged. This matches the recorded reference values
already in the tree:
- `specs/SPEC_excited_state_certification.md:194`: "E₁ = E₂ = −29.085033 Ha"
- `specs/BACKLOG.md:570`: "E₁ = E₂ = −29.085033 Ha"

Both match `-29.08503265` rounded to 6 places — no document update needed, no number silently
changed.

## Gate run (verbatim, pass/fail counts)

`excited_bounds.py` is a demo script (`if __name__ == "__main__":` guarded), not imported by
`tests/test_excited_bounds_spec.py` (confirmed by reading the test file's imports — it only imports
the library functions, not the `__main__` block), so the fix is not gated by that spec today. Ran it
anyway to confirm nothing else regressed:

```
$ uv sync --extra test
$ uv run pytest -q tests/test_excited_bounds_spec.py
..........                                                               [100%]
10 passed, 4 warnings in 101.16s (0:01:41)
```

Full end-to-end run of the modified script (`uv run python excited_bounds.py`, exit code 0):

```
Be2 CAS(4,8)/cc-pVDZ at R=2.45 A: singlet CASCI roots [-29.184127 -29.085033 -29.085033] --
E_1 = E_2 (degenerate pi pair), so NO separator theta_1 < beta <= E_2 exists
   M | theta_1 - E_1^CASCI (mHa) | L_1        | rank
   6 |                   207.798 |   -29.1870 | 3
   8 |                   207.767 |   -29.1871 | 3
  10 |                   146.802 |   -29.1849 | 3
  12 |                   147.011 |   -29.1849 | 3
  14 |                   147.257 |   -29.1849 | 3
Be2 refuses at every depth: ...
wrote data/excited_certification_bench.csv (69 rows)
```

Identical to a baseline run against the pre-fix code (same roots, same CSV row count) — confirmed
by diffing the printed root line only, since `data/` is gitignored and not a tracked artifact.

Lint:
```
$ uv sync --extra dmrg --extra test --extra dev
$ uv run ruff check excited_bounds.py
All checks passed!
```

Repo-wide grep confirming the acceptance criterion's underlying class of bug is now fully closed
(chem-3f8 fixed all `direct_spin1.kernel` sites; this bead was the last `direct_spin0.kernel` site):

```
$ grep -rn "direct_spin0.kernel\|direct_spin1.kernel" --include=*.py .
(no output)
```

## Pre-registered criteria (from the bead)

> `fci.direct_spin0.kernel` call in `excited_bounds.py`'s `__main__` block is replaced with a solver
> whose `.converged` is checked for all requested roots (raises or reports on failure), and any
> recorded Be2 root values that change are noted, not silently updated.

- Replaced with a checked solver: **PASS** (see diff above).
- Recorded Be2 root values that change are noted, not silently updated: **PASS** — no values
  changed (bit-identical to 8 decimal places between old and new code paths), so there was nothing
  to note beyond confirming the match, which is done above.

## What I decided not to do, and why

- Did not build a shared `nroots`-aware wrapper in `hybrid_quantum_solver/dmrg_reference.py`
  alongside `fci_energy`. This is the only call site in the repo needing `nroots`-aware convergence
  checking (confirmed by the repo-wide grep above), so a new shared abstraction would be unused
  elsewhere — inlining the ~6-line check at the one call site follows the existing
  `nbn_low_spin.py:spin_fci` precedent without adding speculative surface area.
- Did not add a new `tests/test_*_spec.py` gate. The bead notes this code path is demo-script
  printing under `__main__`, not imported by any spec test, and the bead's acceptance criteria do
  not ask for a new gate — only that the call site be fixed. Filing a "add regression coverage for
  excited_bounds.py's Be2 demo" bead was considered but the bead's own text frames this as low
  urgency, so I did not create additional scope beyond what was asked.

## What I could not verify

- Nothing outstanding. Both the isolated FCI comparison and the full script run reproduce the exact
  pre-fix numbers, `test_excited_bounds_spec.py` passes 10/10, ruff is clean, and the repo-wide grep
  confirms no unchecked `direct_spin{0,1}.kernel` call sites remain anywhere.

## Hardware / wall time

Run on `Linux cd1439ac6e30 7.0.12-linuxkit x86_64`, 8 cores. Full `excited_bounds.py` run
(LiH/H4/N2/H2 sweeps + Be2 CASCI): ~2 min wall time. `test_excited_bounds_spec.py`: 101 s
(10 passed). Both under the "report wall time if over a minute" threshold — reported above.

## Git status — left dirty per policy (do not commit/push)

Only file touched by this bead: `excited_bounds.py` (diff above). All other modified/untracked
files in the working tree predate this session (other beads' in-progress work) and were left
untouched.

Suggested commands for a human to run:
```bash
git add excited_bounds.py
git commit -m "excited_bounds.py: check FCI Davidson convergence for Be2 nroots=3 CASCI (chem-8hq)"
```
