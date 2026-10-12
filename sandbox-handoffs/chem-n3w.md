# chem-n3w — fix the SCF residue at its source: symmetric/tight SCF in `build_molecular_hamiltonian`

Status: **closed, fixed**. Acceptance criteria (pre-registered in `specs/SPEC_symmetric_scf.md`
before any gate ran) met exactly as stated.

## What changed

- `hybrid_quantum_solver/molecular_hamiltonian.py` — added `_point_group(atom, basis, charge,
  spin)` (PySCF `gto.M(..., symmetry=True).groupname`, geometry-only, no SCF) and an opt-in
  `symmetric_scf: bool = False` keyword on `build_molecular_hamiltonian`. When `True` and the
  geometry's detected point group is not `"C1"`, `conv_tol` is floored to
  `reachability.TIGHT_SCF_CONV_TOL` (1e-13, duplicated as a literal — the package does not import
  the top-level `reachability.py`) before `PySCFDriver` runs. `build_args["conv_tol"]` now records
  the *effective* tolerance actually used, so `reachability._driver_orbitals`'s rebuild-from-args
  stays faithful. The flag defaults to `False`, so every existing call site (~49, all keyword-arg)
  is byte-for-byte unaffected — verified directly, see below.
- `specs/SPEC_symmetric_scf.md` — the pre-registered spec (goal, approach, G1–G7 gates, out of
  scope, caveats).
- `tests/test_symmetric_scf_spec.py` — the gate file, G1–G7 (below).
- `specs/BACKLOG.md` — linked the new spec from the chem-obf "CLOSED" entry's "Follow-up: fix the
  residue where it belongs" line.

This is the exact fix the bead and `SPEC_eigenstate_reachability.md` §0 asked for: "Where the
residue really belongs — a tight or symmetric SCF in the builder — is a follow-up."

## Pre-registered acceptance criteria and pass/fail

All seven gates in `specs/SPEC_symmetric_scf.md` §5, encoded in `tests/test_symmetric_scf_spec.py`
before being run to green:

- **G1 — defaults pinned.** `conv_tol` default `1e-9`, `symmetric_scf` default `False`. ✅
- **G2 — point-group detection matches PySCF.** `Dooh`/`D2h`/`C1` battery, no SCF run. ✅
- **G3 — no-op on C1.** Asymmetric 3-D H₄ geometry: `symmetric_scf=True` leaves
  `build_args["conv_tol"]` exactly equal to the input and `E0`/`E_HF` exactly equal to
  `symmetric_scf=False`. ✅
- **G4 — harmless on the clean set.** `SPEC_scf_conv_tol.CLEAN`'s 9 geometries: `|ΔE0| < 1e-9` Ha,
  `|ΔE_HF| < 1e-8` Ha between flag on/off. ✅
- **G5 — equivalent to the manual `TIGHT_SCF_CONV_TOL` workaround.** Square H₄, sides 1.05–1.40 Å
  (includes both `SPEC_scf_conv_tol` witnesses): `symmetric_scf=True` from a loose `conv_tol=1e-6`
  input is bit-identical (`ΔE0 = ΔE_HF = Δoverlap = 0`) to an explicit `conv_tol=TIGHT_SCF_CONV_TOL`
  build, because both resolve to the same literal `1e-13` passed to `PySCFDriver`. ✅
- **G6 — residue collapse at a Linux-native witness (definition of done, 1/2).** Thread-pinned,
  square H₄ side 1.2 Å, default `conv_tol=1e-9`: 5 repeated builds without the flag give a
  **repeatable** HF population `p0 = 5.0769158288758025e-09` (`> REACHABLE_TOL_CERTIFIED` =
  `1e-10`, admitted by `reachable_mask` and flagged with a `RuntimeWarning` every time); 5 repeated
  builds with `symmetric_scf=True` give `p0 = 8.479726128617941e-34` (`< 1e-20`, no warning, every
  time). ✅
- **G7 — QKSD lands on the physical level, not the artifact (definition of done, 2/2 — the bead's
  actual claim).** At the same witness, independent PySCF symmetry-adapted FCI gives `E_Ag =
  -4.20594019906833` Ha, `E_B1g = -4.355115240509856` Ha (electronic, B1g is the *lower*,
  HF-symmetry-forbidden level here). Without the flag, `QuantumKrylovSolver` from `|HF⟩` converges
  (M ∈ [28,32], the pre-registered G9 window from `SPEC_eigenstate_reachability`) to within `1e-9`
  Ha of `E_B1g` and more than `0.1` Ha below `E_Ag` — reproducing chem-obf's G9 finding at this
  geometry. With `symmetric_scf=True`, every sampled QKSD energy in that window matches `E_Ag` to
  `< 1e-8` Ha, and `reachable_eigenpairs(mh)[0][0]` equals `E_Ag`, not `E_B1g`. ✅

Verbatim result:
```
$ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    uv run pytest tests/test_symmetric_scf_spec.py -v
...
======================= 23 passed, 45 warnings in 48.68s =======================
```
(23 = 1 G1 + 4 G2 + 1 G3 + 9 G4 + 6 G5 + 1 G6 + 1 G7 parametrizations. The lone warning that
matters — `reachability.py:310` `RuntimeWarning` — fires exactly once, inside G7's `off` case, as
pre-registered; the rest are the pre-existing `qiskit_nature` deprecation notice.)

Thread-pinning is required for G6/G7 for the same already-documented reason
`tests/test_hamiltonian_determinism.py` (chem-8tl) pins it: multithreaded BLAS reduction order is a
separate source of run-to-run bit noise in this exact SCF (see "What could not be verified" below
for the un-pinned behavior I measured directly).

## Every number produced, and where it comes from

- `p0` (HF population), `E0`, `E_HF`, `build_args["conv_tol"]` for every geometry above: produced
  by `tests/test_symmetric_scf_spec.py` itself (`_build`, `_dense_hf_projection`,
  `exact_reachable_overlap`) — rerun the file to regenerate.
- `E_Ag`, `E_B1g` (independent FCI reference): `test_eigenstate_reachability_spec._fci`, imported
  directly by the new gate file (not duplicated) — this is the same reference chem-obf's G9 used.
- The "bit-identical to pre-change" numbers in the next section are not vendored in any file; they
  were a scratch check during this session, reproducible with the commands shown there.

## Regression check: default-argument call sites are provably unaffected

Since `symmetric_scf` defaults to `False`, the new code's only effect on the default path is that
`_point_group` is *not called at all* (short-circuit in `if symmetric_scf and ...`) — so
`effective_conv_tol == conv_tol` unconditionally and the rest of the function body is character-for
-character what it was before. Verified directly rather than argued: built 4 representative
geometries (H2 0.74, the two `SPEC_scf_conv_tol` witnesses — square H4 1.35/1.19 — and H4 chain 1.0)
with the pre-change code (`git stash` the one changed file) and again with the post-change code,
both thread-pinned, both via default arguments only:

```
H2 0.74          identical: True
H4 square 1.35   identical: True
H4 square 1.19   identical: True
H4 chain 1.0     identical: True
```
(`(E0, E_HF, build_args["conv_tol"])` tuples compared with Python `==`, not a tolerance — thread
-pinned SCF is deterministic, so this is exact-equality, not "close enough".)

In addition, ran the three core validation gates plus the determinism gate (unrelated to this
change, run for confidence that nothing else moved):
```
$ uv run pytest tests/test_reference_energies.py tests/test_krylov_convergence.py \
    tests/test_noise_resilience.py tests/test_hamiltonian_determinism.py -q
...
21 passed, 14 warnings in 460.21s (0:07:40)
```

Lint: `uv run ruff check .` → `All checks passed!` (whole repo, including the two new files).

## What I decided not to do, and why

- **Did not flip the global `conv_tol` default, or default `symmetric_scf=True`.** Explicitly out
  of scope per the bead and `specs/SPEC_scf_conv_tol.md`'s prior verdict (chem-ayr): the default
  residue's *presence* is pinned by that spec's own G1–G3 and by `SPEC_eigenstate_reachability`'s
  G9, both deliberately. The opt-in design makes this a strict superset of capability with zero
  blast radius, rather than reopening that verdict.
- **Did not implement a true symmetry-constrained SCF** (PySCF irrep-occupation control). A
  tolerance floor is sufficient for the stated goal (residue → machine zero) and reuses the
  already-validated `TIGHT_SCF_CONV_TOL` mechanism verbatim (G5 proves the two are literally the
  same number reaching the same driver call) rather than introducing a new, unvalidated SCF path.
- **Did not run the full unpinned `make gates` / all five potentially-affected spec files
  (`test_scf_conv_tol_spec.py`, `test_eigenstate_reachability_spec.py`,
  `test_symmetry_reachability_spec.py`, `test_reachability_tolerance_spec.py`,
  `test_chained_overlap_spec.py`) to a clean pass/fail diff.** Attempted this first; it is already
  known (chem-dy4, chem-ffe) that these files flake on Linux without thread-pinning, and I
  independently reproduced that directly: the *same* unmodified code gave different fail-sets on
  two consecutive unpinned runs of `test_scf_conv_tol_spec.py` (first run: `test_G2[1.35]`,
  `test_G2[1.19]`, `test_G3[1.35]`, `test_G3[1.19]`, + 2 more, 6 failed; second run on the
  *post-change* tree: `test_G1[1.2]`, `test_G2[1.35]`, `test_G2[1.19]`, `test_G3[1.35]`, 4 failed —
  a different set, including one (`G1[1.2]`) that only failed the second time on a file my change
  cannot touch via default arguments). A run-to-run diff of flaky baselines is not evidence about
  this change either way, so I replaced it with the exact-equality thread-pinned check above, which
  actually isolates the one variable that matters (my code change) from the one that doesn't
  (BLAS thread count). Did not spend further session time chasing a clean unpinned run of
  `test_chained_overlap_spec.py` either (it alone took >10 minutes and was still running when I
  moved on) — logged as still-open flakiness, not blocking this bead's gates.
- **Did not open a new bead for the Linux platform finding.** Added it as notes on `chem-ffe`
  instead (`bd update chem-ffe --notes=...`, done this session) — that bead's stated purpose is
  exactly a Linux verification pass over these witnesses, so a new bead would duplicate it. The
  finding: RHF spontaneous symmetry-breaking at square H4 is essentially *inverted* between this
  Linux container and the macOS numbers the two specs were measured on (here, sides
  1.05/1.10/1.19/1.35/1.40 break symmetry and 1.00/1.20 stay symmetric — roughly the opposite
  split), which is why this bead's own witness is side 1.2, not the macOS witnesses 1.19/1.35.

## What could not be verified

- **macOS behavior is not re-verified in this session** (container is Linux-only). The fix's
  mechanism (`effective_conv_tol = min(conv_tol, TIGHT_SCF_CONV_TOL)`) is platform-independent by
  construction, and G5 shows it is literally identical to the existing, macOS-validated
  `TIGHT_SCF_CONV_TOL` manual workaround, but no gate in this change ran on macOS hardware.
- **The exact unpinned (non-thread-pinned) pass/fail state of `test_scf_conv_tol_spec.py`,
  `test_eigenstate_reachability_spec.py`, `test_symmetry_reachability_spec.py`,
  `test_reachability_tolerance_spec.py`, `test_chained_overlap_spec.py` on this tree** — not
  pinned down to a single clean number, because (as documented above and on `chem-ffe`) that
  number is not stable run-to-run on this container even without this change. What *is* verified
  is that my change does not alter the default-argument code path at all (byte-for-byte, thread
  -pinned), which is the strongest available guarantee given that pre-existing instability.
- **`.beads/issues.jsonl` auto-export.** `bd update chem-ffe --notes=...` printed `Warning:
  auto-export skipped: failed to check existing JSONL: invalid character '<' looking for beginning
  of value`. `git status` shows no changes to `.beads/` despite the Dolt-side note having been
  added (confirmed via `bd show chem-ffe`), so the JSONL export is stale relative to Dolt. This
  predates this session's work (not caused by anything here) and is out of this bead's scope; left
  as-is per the git policy (no `.beads/` export-regeneration command was given).

## Commands to reproduce the key numbers

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    uv run pytest tests/test_symmetric_scf_spec.py -v
uv run pytest tests/test_reference_energies.py tests/test_krylov_convergence.py \
    tests/test_noise_resilience.py tests/test_hamiltonian_determinism.py -q
uv run ruff check .
```
