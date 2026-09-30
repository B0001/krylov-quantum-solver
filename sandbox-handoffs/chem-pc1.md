# Handoff: chem-pc1 — bound the higher roots of PDS's P_K(E) as excited-state estimates

**Status: bead closed.** Definition-of-done gate (G1) passes; the pre-registered usefulness check
(G2) and the pre-registered degeneracy check (G4) both come back **killed**, exactly as the bead's
author predicted — plus one *unplanned* finding (R3 below) discovered while writing G1 that goes
beyond what the bead asked for.

## What changed

- `moment_expansion.py`: added `pds_roots(moments, order, offset=0.0) -> np.ndarray` (all real roots
  of `P_K(E)`, ascending, each lifted by `offset`). `pds_energy` is now `pds_roots(...)[0]` — a pure
  refactor, no behavior change (verified: `tests/test_moment_pds_spec.py` G1-G4 still pass, `python
  moment_expansion.py` reproduces the exact numbers documented in `SPEC_moment_pds.md`'s status
  line: H4 67.8/10.3/2.0/0.43 mHa over K=1..4).
- `specs/SPEC_pds_excited_roots.md` (new): the spec for this work, including the R1-R3 caveats.
- `tests/test_pds_excited_roots_spec.py` (new): gates G1, G1b, G2, G3, G4.
- `specs/SPEC_moment_pds.md`: §7's "excited states from the higher roots" out-of-scope line now
  points at the new spec instead of being a dangling unimplemented follow-up.
- `specs/BACKLOG.md`: the "higher roots of P_K(E)" Method-rungs entry flipped `[ ]` → `[x]` with the
  measured outcome (mirrors the spec status line).
- `scratch_pds_roots_scout.py` (untracked, throwaway): the exploratory script that produced the
  numbers below before the gates were written (pre-registration, per repo convention).

## Gate-run lines (verbatim)

```
$ uv run pytest -v tests/test_pds_excited_roots_spec.py
tests/test_pds_excited_roots_spec.py::test_G1_variational_bound_holds_for_every_root PASSED
tests/test_pds_excited_roots_spec.py::test_G1b_bound_can_break_at_K7_on_near_degenerate_stretched_h4 PASSED
tests/test_pds_excited_roots_spec.py::test_G2_K8_first_excited_root_misses_chemical_accuracy PASSED
tests/test_pds_excited_roots_spec.py::test_G3_ground_root_converges_orders_of_magnitude_faster_than_first_excited PASSED
tests/test_pds_excited_roots_spec.py::test_G4_index_matching_skips_a_level_on_near_degenerate_stretched_h4 PASSED
5 passed, 6 warnings in 2.8s   (verified stable over >=5 independent fresh-process runs)

$ uv run pytest -q tests/test_moment_pds_spec.py   # regression check on pds_energy refactor
4 passed, 10 warnings in 47.1s

$ uv run ruff check moment_expansion.py tests/test_pds_excited_roots_spec.py scratch_pds_roots_scout.py
All checks passed!
```

Both gate files run standalone in their own process (pyscf/qiskit only, no block2 — `make gates`
already isolates every `test_*_spec.py`). Hardware: `Linux 0c74eec296aa 7.0.14-linuxkit x86_64`,
8 vCPU (`nproc`), all runs completed in well under a minute — no long-running-job caveat applies.

## The pre-registered criteria, and what happened to each

From the bead:

1. **`pds_roots` (or equivalent) returns the full sorted real-root array.** Done —
   `moment_expansion.pds_roots`.
2. **Validity: `root_j >= E_j - 1e-9 Ha` for all j, checked on H4 across K=3..8.** Checked. **Result:
   more nuanced than "holds everywhere."** It holds on H4 equilibrium-ish (`1.0/2.0/3.0 Å`) across
   the full K=3..8, verified over 8+ independent fresh-process trials with BLAS pinned to 1 thread
   (`OMP_NUM_THREADS=1` etc., same pattern as `tests/test_shift_both_sides_spec.py`). On the
   near-degenerate stretched H4 (`2.0/4.0/6.0 Å`, chosen for the next criterion), it holds K=3..6
   but **is reproducibly violated at K=7 on about half of fresh-process runs, by -7.49 mHa** —
   `cond(M) ≈ 1.2e17` there. This is not roundoff dust (1e-9 Ha scale); it's the variational
   guarantee, a theorem about exact arithmetic, actually failing in float64 once conditioning and
   near-degenerate orbitals combine. G1 is scoped to the regime where the bound is empirically
   reliable (K=3..8 equilibrium-ish, K=3..6 stretched); G1b pins `cond(M) > 1e14` at the excluded
   cell as the falsifiable reason it's excluded, so the exclusion itself can be checked, not just
   asserted in prose. **This was not predicted by the bead** — the bead's own killable check assumed
   the bound survives and only usefulness would die. See `SPEC_pds_excited_roots.md` R3.
3. **Usefulness: checked against the 1.6 mHa-at-K=8 threshold.** Checked, **killed as predicted.**
   PDS(8)'s root_1 (first excited state) on H4 equilibrium-ish misses reachable E_1 by 12-14 mHa
   (varies a few mHa across fresh interpreter processes at this conditioning — see R1), 8-9x over the
   bar. `test_G2` asserts the miss with a 5 mHa margin, safely clear of that run-to-run float noise.
   `test_G3` quantifies the accompanying "much slower than the ground root" claim: at K=7, the ground
   root is ~11400x more converged than root_1 from the *same* linear solve (0.0017 mHa vs 19.8 mHa)
   — sharper than the backlog's own "~1000x" estimate.
4. **Degeneracy/index-matching on near-degenerate stretched H4: explicitly tested and reported,
   including any index-skip failure mode.** Done, **killed exactly as the bead's author predicted.**
   At K=6 (chosen: `cond(M) ≈ 1.5e14`, inside the reliable double-precision range, unlike K=7/8's
   `>1e17`), nearest-neighbor assignment of the first 6 roots to their closest reachable level is
   **not injective** — reachable level index 2 (`E=-1.3689`) is skipped; the corresponding root
   (`-1.3500`) lands closer to level 3 (`-1.3542`) instead. `test_G4` asserts the skip is reproducible
   (verified over 5 independent fresh-process, thread-pinned trials, identical every time) rather
   than asserting the matching is accurate — that inaccuracy *is* the finding.

## Every number produced, command to regenerate, where it's vendored

All numbers below come from `scratch_pds_roots_scout.py` (exploratory, pre-registration) and the
gates in `tests/test_pds_excited_roots_spec.py` (the numbers actually asserted on). Regenerate with:

```
uv run python scratch_pds_roots_scout.py          # scout numbers (not gated)
uv run pytest -v tests/test_pds_excited_roots_spec.py   # gated numbers
```

- H4 eq K=7: err_root0 = 0.0017 mHa, err_root1 = 19.8 mHa, ratio ≈ 11400 — vendored in
  `SPEC_pds_excited_roots.md` §5 G3 and `test_G3`'s docstring/assertion (`>= 100`).
- H4 eq K=8: err_root1 = 12-14 mHa (process-dependent, always >5 mHa) — vendored in
  `SPEC_pds_excited_roots.md` §5 G2/R1 and `test_G2`'s assertion (`> 5e-3`).
- H4 eq K=3..8 `cond(M)`: 3.4e6 → 2.5e20 — vendored in `SPEC_pds_excited_roots.md` §2.
- H4 stretched K=6 index-skip: `assigned=[0,1,3,7,7,7]`, `skipped=[2,4,5]` (first 6 roots vs first 8
  reachable levels) — vendored in `SPEC_pds_excited_roots.md` §5 G4 and `test_G4`'s assertion.
- H4 stretched K=7 bimodal outcome: `+196 μHa` margin (3/6 trials) vs `-7.49 mHa` violation (3/6
  trials), `cond(M) ≈ 1.2e17` — vendored in `SPEC_pds_excited_roots.md` R3 and `test_G1b`'s
  `cond > 1e14` assertion (the bimodality itself is documented in prose, not gated, since it is
  inherently process-state-dependent — see "what I could not verify" below).

## What I decided not to do, and why

- **Did not attempt the centered-moment-frame conditioning fix.** `specs/BACKLOG.md` has a separate,
  still-open entry for that (chem-70c: "PDS's accuracy floor is a removable artifact of the
  uncentered moment frame"). The bead's own caveat explicitly says this higher-roots work should be
  "ordered after (or alongside)" that entry since raw K=7 is already conditioning garbage — I did it
  alongside, in the sense of measuring and reporting the raw-frame behavior honestly (including how
  bad it gets), but did not implement the centering fix itself. That is out of scope for chem-pc1 and
  would double the surface area of this bead. Filed nowhere new since chem-70c already exists and
  covers it.
- **Did not build a general degeneracy-robust matcher** (e.g. Hungarian assignment with an
  abstain/no-match option). The bead asked to test and report the index-skip failure mode, not to
  fix it. `SPEC_pds_excited_roots.md` §7 names this as a possible follow-up, not filed as a new bead
  (small, and not clearly wanted — the finding itself, "roots track density not index," may be the
  more valuable output than a patched matcher).
- **Did not probe LiH or H2 for the degeneracy effect.** H2's reachable space at this active space is
  only 2-dimensional — too small to show an index-skip. LiH wasn't probed; H4 alone was sufficient to
  answer both pre-registered checks and stayed within "cheap" per the bead's own cost estimate.
- **Did not run the full `make gates`/`scripts/run_gates.sh` suite.** Out of scope for verifying this
  one bead's two relevant files (`test_moment_pds_spec.py`, `test_pds_excited_roots_spec.py`), and
  the DMRG/block2 gates are unrelated to this change (no block2 import anywhere in
  `moment_expansion.py` or the new test file).

## What I could not verify

- **The exact root cause of the K=7 stretched-geometry bimodality (R3) is not pinned down.** I
  confirmed it's reproducible (same split, same two outcomes, across >10 fresh-process trials, both
  with and without BLAS thread-pinning), and that it correlates with `cond(M) > 1e14` combined with
  near-degenerate orbitals, but I did not isolate whether the split comes from PySCF's CASCI orbital
  gauge on degenerate orbitals (the same effect `tests/test_shift_both_sides_spec.py`'s own docstring
  describes for a different quantity), `np.roots`' companion-matrix eigensolve picking up
  process-dependent floating-point paths, or both. Stated as an open question in
  `SPEC_pds_excited_roots.md` R3, not asserted as solved.
- **Did not verify behavior on macOS or with a GPU BLAS backend** — this container is Linux x86_64
  only; the conditioning numbers (and possibly the K=7 bimodality) could differ on Accelerate/MKL.
  Not claimed either way.

## Git

Tree is dirty, not committed, per this session's git policy. Files touched:

```
 M moment_expansion.py
 M specs/BACKLOG.md
 M specs/SPEC_moment_pds.md
?? specs/SPEC_pds_excited_roots.md
?? tests/test_pds_excited_roots_spec.py
?? scratch_pds_roots_scout.py
?? sandbox-handoffs/chem-pc1.md
```

Suggested commit (not run):

```
git add moment_expansion.py specs/BACKLOG.md specs/SPEC_moment_pds.md \
        specs/SPEC_pds_excited_roots.md tests/test_pds_excited_roots_spec.py \
        scratch_pds_roots_scout.py
git commit -m "pds: add pds_roots for excited-state bound; kill usefulness+index-matching claims (chem-pc1)"
```

Note: several other files were already modified/untracked at session start (from a prior session,
unrelated to this bead — e.g. `certified_dipole.py`, `lambda_ladder.py`, `sandbox-handoffs/chem-52i.md`,
various `scratch_*_52i.py`). I did not touch or include those in the suggested commit above.
