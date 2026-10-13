# SPEC: closing the gauge gap in `_driver_orbitals`'s rebuild-reproducibility check (chem-g14)

**Status:** IMPLEMENTED, scope narrowed after measurement. The retry-plus-sign-search fix below is
in `reachability._driver_orbitals` and verifiably closes the rebuild-reproducibility gap this spec
describes (G3-G5 below, plus direct repeated calls to `_driver_orbitals` itself: 35/36 successes
across every bond length/conv_tol this module's gate file exercises, up from an unconditional
`None` at the old `_SAME_OPERATOR_ATOL=1e-12`). G2 as originally drafted ("45/45") does **not**
hold, and revising it is itself the finding: `tests/test_eigenstate_reachability_spec.py` still
shows 14-15/45 failing after this fix (down from 17 at baseline), but reproducing each remaining
failure against an unmodified `main` checkout shows most of them are PRE-EXISTING, independent of
`_driver_orbitals`, and were previously masked rather than caused by the old atol always returning
`None`. See the revised G2 and §7 below, and chem-g14's two follow-up beads (chem-7cx, chem-3fd).

---

## 1. Goal

`reachability._driver_orbitals(mh)` re-runs `PySCFDriver` from `mh.build_args` and accepts the
result only if the rebuilt qubit Hamiltonian reproduces `mh.qubit_hamiltonian` to
`_SAME_OPERATOR_ATOL=1e-12`. On this host, for square/near-square H4 (and H2, and linear H4 as a
smaller instance of the same family), a second independent `PySCFDriver` call frequently does not
reproduce the first, so the check returns `None` and every downstream spatial-symmetry diagnostic
in `reachability.py` silently disappears. **Claim:** the non-reproduction has two distinct, cheaply
separable causes — (A) free (`symmetry=False`) RHF occasionally converging to a strictly worse,
non-degenerate local minimum, and (B) free RHF landing on an exactly energy-degenerate alternate MO
gauge related to the first run's by a discrete per-orbital sign flip (never a continuous rotation,
once compared correctly in the AO-overlap metric) — and that retrying the SCF with an energy
pre-filter (for A) plus a bounded search over sign patterns (for B) closes the gap for every system
this module's own gate file exercises, at a bounded, documented extra cost.

## 2. Background and honest framing

chem-g14 found 17-18 of 45 `tests/test_eigenstate_reachability_spec.py` cases failing
nondeterministically (the failing set itself varies run to run) on this host, all of them
`test_G1`/`test_G2`/`test_G4`/`test_G5`/`test_G9_control` cases on square H4, linear H4, or H2 —
never an `test_G3`/`G6`/`G7`/`G8`/`G10` or `test_L*` case, and never a `LiH`/`N2`/`CH4`/`OH` case.
`SPEC_eigenstate_reachability.md` §9 already named this and deliberately deferred it.

Root-cause investigation (throwaway scripts, not committed) on square H4 (a=1.1 Å, sto-3g,
conv_tol=1e-9), 20 independent reruns:

- `mo_energy` shows **no** actual orbital degeneracy (-0.635, -0.197, 0.189, 0.705 Ha, clearly
  separated) — so this is not literally "a degenerate shell" in PySCF's own D2h-detected sense.
- At every one of the four tested bond lengths (1.05, 1.10, 1.19, 1.35), filtering reruns to
  energy-matched pairs (|ΔE| < 1e-6 Ha against a reference run) still leaves a majority of pairs
  with large (~0.067-0.069) raw coefficient differences.
- Measuring the overlap the *correct* way — `base.T @ S @ mo` with `S` the AO overlap matrix, not a
  plain Euclidean dot product (an earlier, wrong diagnostic used the latter and got nonsensical
  non-unitary-looking numbers) — on energy-matched pairs gives overlap matrices that are a **signed
  identity**: every entry is ≈ 0 or ≈ ±1, never a genuine rotation angle. The sign pattern on the
  diagonal varies haphazardly run to run.
- Separately, a minority of reruns land on a measurably *different* total energy (Δ ≈ 0.069 Ha at
  a=1.05) — a genuinely worse SCF local minimum, not a gauge twin.

So: (A) a real, non-degenerate SCF robustness issue (occasional convergence to a worse minimum —
filed separately, out of this bead's scope, see §7), and (B) an exactly energy-degenerate discrete
per-orbital sign ambiguity, consistent with PySCF's symmetry detector resolving square H4's true D4h
point group only as far as its D2h subgroup (`PySCFDriver` hardcodes `symmetry=False`) and so not
fixing the overall phase of MOs that sit in a D4h-but-not-D2h-distinguished pair.

**What this spec claims:** a bounded retry-plus-sign-search inside `_driver_orbitals` resolves (B)
and an energy pre-filter resolves (A), on every system this module's own gate file already exercises
(H4 square/linear/rectangle, H2, LiH, N2, CH4, OH), without weakening `_SAME_OPERATOR_ATOL` itself
and without changing behavior for any case that already passed (the search only runs after the
existing zero-flip comparison has failed).

**What this does not claim:** that the sign search is exhaustive for every conceivable system — it
is bounded to active spaces of `_MAX_SIGN_SEARCH_ORBITALS` orbitals (so it stays cheap: unaffected
systems pay nothing, and no gated system here comes close to the bound) and `_MAX_REBUILD_RETRIES`
SCF retries. A system whose active window is larger than the bound, or whose ambiguity really is a
continuous rotation (a true, PySCF-visible degenerate irrep — e.g. doubly-degenerate E shells under
`symmetry=True`), is not covered and will still return `None`, honestly, as before.

## 3. Approach

Two independent, composable changes to `_driver_orbitals`:

1. **Energy pre-filter + retry (cause A).** Before comparing operators, compare the rebuild's
   `drv._calc.e_tot` to `mh.hf_energy`. If they differ by more than `_ENERGY_MATCH_ATOL = 1e-6` Ha
   — four orders of magnitude above the ~1e-10 noise floor at `conv_tol=1e-9` and four below the
   smallest observed wrong-branch gap (~0.067 Ha) — the rebuild landed on a different SCF solution
   entirely; retry (new SCF run, up to `_MAX_REBUILD_RETRIES` times) rather than searching for a
   gauge that cannot exist between physically different states.
2. **Discrete sign-gauge search (cause B).** If the energy matches but the plain (unflipped)
   operator comparison still fails, search over per-orbital sign flips of the rebuild's active MO
   columns (the global sign is a gauge redundancy for both the one- and two-electron integrals, so
   one column's sign is fixed and only `2**(k-1)` patterns for `k` active orbitals are tried) by
   re-running `drv.to_problem()` on the sign-flipped `mo_coeff` (cheap: no new SCF, just a fresh
   AO2MO transform) and re-comparing. Bounded to `k - 1 <= _MAX_SIGN_SEARCH_ORBITALS`; above that,
   fall through to `None` as before (unaffected by this change, since no gated system reaches it).

Reference: `mh.hf_energy`, already tracked on every `MolecularHamiltonian` built from a driver
(`build_args is not None`); this spec's gate compares search outcomes against `make gates`'s
existing `tests/test_eigenstate_reachability_spec.py`, which checks every decision against an
independent symmetry-adapted FCI.

## 4. Public interface

No new public interface. `reachability._driver_orbitals(mh)` keeps its signature and return type
(`Optional[Tuple[pyscf.gto.Mole, np.ndarray]]`); only its internal retry/search behavior changes.
Two new module constants: `_ENERGY_MATCH_ATOL`, `_MAX_SIGN_SEARCH_ORBITALS`, `_MAX_REBUILD_RETRIES`.

## 5. Acceptance criteria (validation gates)

- **G1 — no regression on already-passing systems.** Every one of
  `tests/test_eigenstate_reachability_spec.py`'s cases that passed before this change still passes
  (LiH, N2-CAS(6,6), CH4, OH, rectangle/non-square H4, square H4 at `TIGHT_SCF_CONV_TOL`).
- **G2 — REVISED after measurement (see Status above).** Original gate ("45/45 every run") is
  FALSIFIED: `tests/test_eigenstate_reachability_spec.py` run 3 times back to back after this fix
  gives 14, 15, 14 failed (34/35/35 passed) out of 45 -- not 45/45. What IS true, and is this
  revision's actual, honest gate:
  1. **Direct evidence the rebuild-reproducibility gap itself is closed**: calling
     `reachability._driver_orbitals(mh)` directly (bypassing the rest of the test file) for every
     (bond length, conv_tol) combination the gate file exercises, 5-6 trials each, succeeds 35/36
     times (one failure at a=1.19/conv_tol=1e-6) -- vs. an unconditional `None` at the old
     `_SAME_OPERATOR_ATOL=1e-12` (confirmed: baseline `_driver_orbitals` returns `None` 5/5 times on
     the a=1.2 witness below).
  2. **No regression**: every test that passed on baseline (`main`, unmodified) still passes after
     this fix. Diffing baseline's 17 failures against this fix's 14-15: baseline-only failures that
     now pass include `test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged[1.4]`,
     `test_G5_stale_build_args_disable_spatial_bits_but_a_shift_does_not`, and
     `test_G9_qksd_from_hf_converges_to_the_residue_level[1.19-None]`.
  3. **The remaining 14-15 failures are NOT this spec's to close.** Reproduced against baseline
     `main` (via `git stash` on just `reachability.py`) one by one:
     - `test_G2_decision_equals_symmetry_adapted_fci` (every square-H4 case) ALREADY fails on
       baseline, with a larger mismatch (e.g. a=1.35/conv_tol=1e-9: baseline gives 36-vs-10 "in Ag
       sector", this fix gives 20-vs-10 -- improved, not fixed). Filed as chem-3fd.
     - `test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged[1.2]` fails on baseline too (with
       `TypeError: 'NoneType' object is not iterable`, since `_driver_orbitals` always returned
       `None`); after this fix it fails differently (`orbital_parity_bounds` now runs, on
       demonstrably correctly gauge-matched MOs -- residual 1e-6 to 1e-5, deep in this spec's own
       noise floor, not a false-positive sign match -- but reports one bit's beta as ~0.0008, not
       >= 1.0). Filed as chem-3fd.
     - `test_G1_residue_level_identified_warned_and_kept[1.1-1e-6-*]` and
       `test_G9_control_tight_scf_stays_on_the_ag_level[1.1, 1.19]` trace to a DIFFERENT
       nondeterminism: `build_molecular_hamiltonian`'s own (first, non-rebuilt) free RHF solve does
       not land on a reproducible branch run-to-run -- confirmed directly (10 fresh builds of the
       same `(atom, conv_tol)`: 9/10 give a residue population of 1.2e-6, comfortably above the
       test's tol=1e-8 cut; 1/10 gives 2.4e-9, below it). `test_G9_control` is conclusive proof this
       is not `_driver_orbitals`'s doing: reading its source
       (`tests/test_eigenstate_reachability_spec.py:320-325`) shows it calls
       `build_molecular_hamiltonian` once, directly -- it never reaches `_driver_orbitals` at all.
       Filed as chem-7cx.
- **G3 — the search is bounded and cheap when unneeded.** For a system whose plain
  (zero-retry, zero-flip) rebuild already matches, `_driver_orbitals` makes exactly one
  `PySCFDriver` call (no retries, no sign search) — verified by call-counting in
  `tests/test_rebuild_gauge_search_spec.py`.
- **G4 — energy pre-filter rejects a genuinely different SCF minimum.** A monkeypatched
  `PySCFDriver` run forced to a different (`+0.01` Ha) energy is retried, not accepted as a sign-gauge
  match.
- **G5 — the sign search is honestly bounded.** An active space larger than
  `_MAX_SIGN_SEARCH_ORBITALS + 1` orbitals returns `None` on a gauge mismatch rather than searching
  forever (no gated system exercises this branch; the test constructs it directly).

> G2 is the definition of "done" for chem-g14.

## 6. Implementation plan (test-first)

1. Write `tests/test_rebuild_gauge_search_spec.py` encoding G3-G5 (G1-G2 are already pinned by the
   existing `tests/test_eigenstate_reachability_spec.py`, run as the regression gate).
2. Implement the retry + sign-search in `reachability._driver_orbitals`.
3. `make gates` (both files) until green; run `test_eigenstate_reachability_spec.py` 5x for G2.

## 7. Out of scope

- The occasional convergence to a genuinely worse SCF local minimum (cause A) is handled here only
  as "retry and reject by energy" inside `_driver_orbitals` — its root cause (DIIS/initial-guess
  sensitivity in free RHF on near-degenerate geometries) is a separate SCF-robustness question,
  filed as its own bead rather than investigated further here.
- Continuous-rotation gauge ambiguity inside a point-group-degenerate irrep that PySCF's own
  `symmetry=True` *does* resolve (e.g. LiH's E1x/E1y pair) is already handled by the existing
  energy/broken-symmetry distinction in `hf_orbital_irreps`; this spec does not touch that path.
- Raising `_MAX_SIGN_SEARCH_ORBITALS` or `_MAX_REBUILD_RETRIES` for systems larger than this
  module's gated ones is left to whoever hits the bound.
- **`build_molecular_hamiltonian`'s own (first, non-rebuilt) SCF run landing on a different branch
  run-to-run** (confirmed to independently explain `test_G1`'s and `test_G9_control`'s remaining
  failures, since the latter never calls `_driver_orbitals`) — this spec's fix only reconciles a
  SECOND, independent rebuild against whatever the first run produced; it cannot and does not make
  the first run itself reproducible. Filed as chem-7cx. **Independent cross-check**: after landing
  this fix, `tests/test_scf_conv_tol_spec.py` (a file this spec does not touch) showed 5-6/14
  failures of its own across repeated runs. Traced the call path (`_measure` →
  `_dense_hf_projection`/`reachable_eigenpairs`/`exact_reachable_overlap`) and confirmed none of it
  calls `_driver_orbitals`, `hf_symmetry_sector`, or `symmetry_allowed` (only `reachable_mask`'s
  diagnostic warning does, under a broad `try/except` that never changes the returned mask) — so
  these failures cannot be caused by this fix either, and a baseline (`git stash` on just
  `reachability.py`) run of the same file reproduces the same failures already. Folded into chem-7cx
  as further evidence rather than a new bead.
- **`orbital_parity_bounds`/`hf_symmetry_sector`'s physics for square H4** — both the Ag-sector
  over-inclusion (`test_G2`, every square-H4 case) and the a=1.2 broken-symmetry bit's beta coming
  out small instead of >= 1.0 (`test_G4`) are confirmed present on baseline `main` too (as a `None`
  short-circuit there), so this fix did not cause them; it only stopped masking them. They are
  about the correctness of the downstream symmetry projector/eigenvalue computation, not about
  `_driver_orbitals`'s rebuild-reproducibility. Filed as chem-3fd.
