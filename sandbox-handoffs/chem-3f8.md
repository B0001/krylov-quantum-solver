# chem-3f8 handoff: unchecked FCI Davidson convergence outside `dmrg_reference.py`

## Bead (verbatim, `bd show chem-3f8`)

> chem-rhw found that PySCF's direct_spin1 default (100 Davidson iterations, no convergence
> check) silently returned H_n FCI energies 0.28 mHa (n=10) and 91 mHa (n=12) too high in a site
> basis. It fixed hybrid_quantum_solver/dmrg_reference.fci_energy (max_cycle=1000, raise on not
> converged). These still call fci.direct_spin1.kernel unchecked: lambda_ladder.py,
> df_factorization.py, thc_factorization.py, shift_both_sides.py, tests/test_scdf_lambda_spec.py,
> tests/test_thc_lambda_spec.py. They are canonical-basis today, so they are probably converged,
> but "exact reference" is only as good as the convergence check. Route them through fci_energy,
> or add the same check.

**Acceptance criteria (verbatim):** "No unchecked direct_spin1 kernel call remains (grep);
affected spec gates still pass; any recorded reference energy that changes is reported, not
silently updated."

## What I changed

Every listed call site, plus one the bead's file list missed but its own grep criterion catches
(`cross_check.py`), routed through the already-vetted `hybrid_quantum_solver.dmrg_reference.fci_energy`
(raises `RuntimeError` if Davidson doesn't converge in 1000 cycles, instead of PySCF's silent
100-cycle default):

- **`lambda_ladder.py`** — `fci_energy_error()` (used by `lambda_ladder()`'s print loop and by
  `tests/test_lambda_ladder_honest_caveat_spec.py`).
- **`df_factorization.py`** — `rank_for_accuracy()` and the `__main__` demo's per-rank FCI sweep.
- **`thc_factorization.py`** — the `__main__` demo's THC-reconstruction FCI check.
- **`shift_both_sides.py`** — `spectrum_preserved()` (the G4 spectrum-invariance gate for the
  symmetry shift).
- **`cross_check.py`** — the CASCI reference in `cross_check()`. **Not in the bead's file list**
  (chem-rhw's audit apparently predated or missed this file), but it matched
  `fci.direct_spin1.kernel(h1, eri, norb, nelec)` unchecked and is exactly the same bug: this is
  the CASCI arm of the four-method cross-check harness, i.e. the thing the *other* three solvers
  are being judged against. Found by running the acceptance criterion's own grep
  (`grep -rn "direct_spin1" --include=*.py .`) rather than trusting the file list literally.
- **`tests/test_scdf_lambda_spec.py`** — `test_G2_shift_preserves_spectrum`.
- **`tests/test_thc_lambda_spec.py`** — `test_G1_exact_reconstruction_and_fci`.

Each edit is mechanical: `e, _ = fci.direct_spin1.kernel(h1, eri, norb, nelec[, ecore=e_core])`
→ `e = fci_energy(h1, eri, nelec[, e_core=e_core])`, dropping the manual `+ e_core` where it was
added after the call (since `fci_energy` takes `ecore` itself and returns the total energy
already). `fci_energy` derives `norb` from `h1.shape[0]` and accepts `nelec` as either an int or
an `(na, nb)` pair, so no call site needed its `nelec`/`norb` plumbing changed beyond that.

`git diff -- lambda_ladder.py df_factorization.py thc_factorization.py shift_both_sides.py
cross_check.py tests/test_scdf_lambda_spec.py tests/test_thc_lambda_spec.py` shows the exact
diff (7 files, ~40 lines).

I did **not** touch `hybrid_quantum_solver/dmrg_reference.py`, `nbn_dmrg_reference.py`,
`be2_cbs.py`, `specs/*`, or `tests/test_hubbard_lieb_wu_spec.py` /
`tests/test_nbn_dmrg_reference_spec.py` / `tests/test_nbn_low_spin_spec.py` — those were already
modified in the working tree when I claimed this bead (a prior session's uncommitted work,
unrelated to this bug). I left them alone.

## Acceptance criterion 1: "No unchecked direct_spin1 kernel call remains (grep)"

```
$ grep -rn "direct_spin1\.kernel\|direct_spin0\.kernel" --include=*.py .
excited_bounds.py:317:    roots = fci.direct_spin0.kernel(h1, eri, 8, ne, ecore=float(e_core), nroots=3)[0]
```

That is the only remaining hit, and it does **not** match the bead's literal grep target
(`direct_spin1`) — it's `direct_spin0` with `nroots=3`, inside `excited_bounds.py`'s `__main__`
demo block, not imported or exercised by any `tests/test_*_spec.py`. `fci_energy` doesn't support
`nroots`, so fixing it needs a small nroots-aware wrapper, not a one-line swap. **Filed as a new
bead, chem-8hq, rather than fixed here** — it's the same bug class but a different call shape and
outside this bead's explicit scope (the bead's file list and its grep target both say
`direct_spin1`).

`grep -rn "direct_spin1\.kernel" --include=*.py .` (the bead's exact pattern) now returns nothing
outside `hybrid_quantum_solver/dmrg_reference.py`'s own `fci_energy`, which doesn't match that
text pattern anyway (it builds `fci.direct_spin1.FCI()` then calls `.kernel(...)` on the instance,
already checked).

## Acceptance criterion 2: "affected spec gates still pass"

Every gate file that imports or exercises a function I changed, run one-per-process per this
repo's isolation rule (none of these touch block2, but I kept the discipline):

```
$ uv run pytest -q tests/test_scdf_lambda_spec.py
....                                                                     [100%]
4 passed in 4.48s

$ uv run pytest -q tests/test_thc_lambda_spec.py
....                                                                     [100%]
4 passed in 5.68s

$ uv run pytest -q tests/test_lambda_ladder_honest_caveat_spec.py
1 failed, 3 passed in 14.11s   # test_G3 -- PRE-EXISTING flake, see below, not caused by this change

$ uv run pytest -q tests/test_shift_both_sides_spec.py
1 failed, 13 passed, 261 warnings in 62.12s   # test_G5[H2] -- PRE-EXISTING flake, see below

$ uv run pytest -q tests/test_cross_check_trust_semantics_spec.py
.....                                                                    [100%]
5 passed in 89.04s (0:01:29)

$ uv run pytest -q tests/test_validate_and_cost_composition_spec.py   # exercises cross_check.py transitively via validate_and_cost.py
....                                                                     [100%]
4 passed in 289.80s (0:04:49)
```

`uv run ruff check lambda_ladder.py df_factorization.py thc_factorization.py shift_both_sides.py
cross_check.py tests/test_scdf_lambda_spec.py tests/test_thc_lambda_spec.py` → `All checks
passed!`.

### The two failures are pre-existing, not caused by this change

`test_lambda_ladder_honest_caveat_spec.py::test_G3_df_rank_truncation_accuracy_is_not_monotonic`
and `test_shift_both_sides_spec.py::test_G5_identity_term_inflates_the_gain[H2]` are flaky across
repeated fresh-process runs — same code, same inputs, different numeric results each process
(non-monotonic vs monotonic error sequences for G3; `inflated`/`honest` ordering flips for G5). I
confirmed both are **unrelated to my edit**: `git stash push -- lambda_ladder.py` (resp.
`shift_both_sides.py`) and re-running against the *original, unmodified* file reproduces the same
failure — e.g. `test_G3` failed 2/5 repeated runs and `test_G5[H2]` failed on the very first
stashed run too, both before any of my changes existed in those files. Root cause is likely RHF/
CASCI orbital-gauge nondeterminism for near-degenerate systems (N2's pi pair, possibly H2 in a
minimal basis), not the Davidson-convergence bug this bead is about. **Filed as a new bead,
chem-1yr** (P2, since it makes two spec gates unreliable) — investigating and fixing that
nondeterminism is out of this bead's scope.

To be doubly sure my `spectrum_preserved` edit itself is solid (that's the function in
`shift_both_sides.py` I actually touched, not the flaky `test_G5`), I ran just its gate (`test_G4`,
which calls `spectrum_preserved`) 3x back to back: `3 passed` every time.

## Acceptance criterion 3: "any recorded reference energy that changes is reported, not silently
updated"

No committed/vendored reference energy changed. None of the touched functions have a value
checked into a spec or table — `fci_energy_error`, `rank_for_accuracy`, and `spectrum_preserved`
are all called live from `__main__` demos or from spec-gate fixtures that recompute everything
from scratch each run (no cached numbers in `specs/*.md` depend on them). I verified the FCI
values themselves are numerically unchanged where PySCF's default already converged (checked
directly for the N2 CAS(4,3) sweep used by the honest-caveat spec: old
`fci.direct_spin1.FCI().kernel(...)` and new `fci_energy(..., e_core=e_core)` gave bit-identical
energies at every rank, `err_old == err_new` to full precision) — the fix only changes behavior
when Davidson would otherwise have silently failed to converge, which none of these canonical-
basis call sites hit in the systems the gates use.

## What I decided not to do, and why

- **`excited_bounds.py`'s `direct_spin0.kernel(nroots=3)`** — same bug class, different call
  shape (multi-root, no `fci_energy` support), not in the bead's file list or its literal grep
  target, not gated by any spec test. Filed as **chem-8hq** instead of fixing inline.
- **The two flaky pre-existing gate failures** (`test_G3`, `test_G5[H2]`) — confirmed unrelated to
  this bug via `git stash`, root cause not yet isolated (orbital-gauge nondeterminism, not FCI
  convergence). Filed as **chem-1yr** instead of trying to fix blind.
- **`nbn_low_spin.py`'s `spin_fci()`** — uses `fci.direct_spin0`/`direct_spin1` but already has
  its own `if not np.all(solver.converged): raise RuntimeError(...)` check (line ~114). Already
  compliant; left untouched.
- **`krylov_subspace_solver.py`'s `fci.direct_spin1.absorb_h1e`/`contract_2e`** — these build a
  dense CI Hamiltonian matrix explicitly (for `np.linalg.eigh`), not a Davidson `.kernel()` call;
  there's no iterative convergence to check. Out of scope by construction, not an oversight.

## What I could not verify

- I did not re-derive or independently sanity-check `hybrid_quantum_solver.dmrg_reference.fci_energy`
  itself (chem-rhw's fix) — I only relied on it being the vetted primitive per `CLAUDE.md` and the
  bead text, as intended.
- I did not run the full `make gates` / `make test` suite (would also re-run the ~90 other specs,
  including slow DMRG ones, well beyond this bead's scope) — only the gates that import or
  exercise the functions this bead's diff touches.

## Verify these commands yourself

```bash
grep -rn "direct_spin1\.kernel" --include=*.py .   # expect: no hits outside dmrg_reference.py's own vetted use
uv run pytest -q tests/test_scdf_lambda_spec.py tests/test_thc_lambda_spec.py tests/test_cross_check_trust_semantics_spec.py
uv run ruff check lambda_ladder.py df_factorization.py thc_factorization.py shift_both_sides.py cross_check.py tests/test_scdf_lambda_spec.py tests/test_thc_lambda_spec.py
```

## Git state

Tree left dirty per the sandbox git policy. My changes are isolated to:
`lambda_ladder.py`, `df_factorization.py`, `thc_factorization.py`, `shift_both_sides.py`,
`cross_check.py`, `tests/test_scdf_lambda_spec.py`, `tests/test_thc_lambda_spec.py`. Suggested
commit (does not include the pre-existing unrelated dirty files from before I started):

```bash
git add lambda_ladder.py df_factorization.py thc_factorization.py shift_both_sides.py \
        cross_check.py tests/test_scdf_lambda_spec.py tests/test_thc_lambda_spec.py
git commit -m "Route unchecked FCI Davidson calls through vetted fci_energy (chem-3f8)"
```

Follow-up beads filed: **chem-8hq** (excited_bounds.py's unchecked `direct_spin0.kernel`,
out of scope), **chem-1yr** (two pre-existing flaky spec gates from orbital-gauge
nondeterminism, unrelated to this bug, discovered while verifying this fix).
