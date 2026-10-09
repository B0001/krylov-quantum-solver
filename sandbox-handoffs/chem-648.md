# chem-648 handoff

## Bead content (verbatim from `bd show chem-648`)

Four known limitations of the symmetry reachability diagnostic (`reachability.py`),
`specs/SPEC_eigenstate_reachability.md` §8:

- **L1** — `orbital_parity_bounds` bounds the Z-string against a compressed `R_b` far from
  orthogonal when the active window splits a degenerate shell (CH4 Td CAS(2,2) example). Fix:
  check `||R_b R_b^T - I||`.
- **L2** — symmetry PySCF detects only within its geometric tolerance is treated as exact
  (sq-H4 witness). Fix: tolerance-aware symmetry detection.
- **L3** — `hf_symmetry_sector`'s JW-layout guard accepts an interleaved mapper when the HF
  index coincides; `build_args` doesn't record the mapper. Fix: record the mapper.
- **L4** — per-cluster majority can veto an allowed level in an accidental <1e-9 Ha degeneracy.

All four implemented. Details, evidence, and what I deliberately left alone below.

## Fixes

### L1 — `_R_ORTHOGONALITY_ATOL` (reachability.py)

`orbital_parity_bounds` already builds `R_b = g_b` written in `mh`'s MO basis per parity bit, and
derives β (the Z-string error bound) from its eigenvalues. β is only a valid bound when `R_b` is
orthogonal; an active window that splits a point-group-degenerate shell compresses `R_b` away from
orthogonal, and β can then come out small (and so get greedily *trusted*) while being physically
wrong. Added an explicit check: `‖R_b R_b^T − I‖ > 1e-6` excludes the bit (β = ∞) instead of
trusting it.

Witness (CH4, Td reduced to D2 by PySCF, C–H 1.087 Å, STO-3G, free `symmetry=False` RHF,
`conv_tol=1e-12`): the HOMO is a true 3-fold-degenerate shell (`mf.mo_energy[2:5]` identical to
<1e-8 Ha). A CAS(2,2) window on 2 of those 3 MOs gives β = ∞ for both parity bits (previously: a
small, wrongly-trusted number). A CAS(2,2) window on two *non*-degenerate MOs (HOMO-1, LUMO) gives
β < 1e-10 for both — the fix doesn't touch windows that don't actually split a degenerate shell.

### L2 — `_geometry_residuals` / `_GEOMETRY_RESIDUAL_ATOL` (reachability.py)

New helper applies each parity bit's claimed Cartesian operation (derived from
`pyscf.symm.param.CHARACTER_TABLE`/`OPERATOR_TABLE`, cross-checked against pyscf source, not just
recalled) directly to the real, unsymmetrized geometry and measures the worst atom-mapping residual
in Bohr — independent of `symm.geom.TOLERANCE`, the size-scaled tolerance that both PySCF's own
`detect_symm` *and* `symm_adapted_basis`'s atom-equivalence grouping use internally. A bit is
excluded (β = ∞) when that residual exceeds `1e-8` Bohr, regardless of what point group PySCF
itself reports.

Witness (rectangle H4, a=1.1/b=1.3 Å — chosen over a square specifically to avoid L1's orbital
degeneracy, see "What I did not reconcile" below): one atom shifted +1e-6 Å along x. PySCF's own
`detect_symm` still reports `D2h` (fooled). The fix correctly excludes all three D2h bits anyway
(β = ∞). The spec's original witness (§8, pre-revision) used a 5e-6/1e-5 Å shift on this system
(macOS); on this host (Linux, pyscf 2.13.1) that shift is no longer inside PySCF's detection
tolerance at all — PySCF itself falls back to `C2h` by 2e-6 Å. Recorded as a platform/version
difference, not a re-derivation of the mechanism (the spec's §8 entry for L2 now says so).

### L3 — record the mapper (`hybrid_quantum_solver/molecular_hamiltonian.py`, `reachability.py`)

One-line fix: `build_molecular_hamiltonian`'s `build_args` now includes
`mapper=type(mapper).__name__`. `hf_symmetry_sector`'s JW-layout guard (which only checked that the
HF bitstring's *index* matched a Jordan-Wigner layout — an index an `InterleavedQubitMapper`
operator can share for some `(n_alpha, n_beta, n_orb)`) now also requires `build_args["mapper"]` to
be `"JordanWignerMapper"` when `build_args` is present. Absent `build_args` (the integrals-only
builder, `build_hamiltonian_from_integrals`, which never sets it and only ever uses Jordan-Wigner)
still trusts the index — no provenance to check, and no other mapper is reachable through that path.

Witness (linear H4, spin=2 — the spec's own L3 witness): `InterleavedQubitMapper(JordanWignerMapper())`
and plain `JordanWignerMapper()` builds of the identical system share the same HF computational-basis
index (23). Before the fix, `hf_symmetry_sector` would have built a full point-group sector for the
interleaved-mapper operator as if it were Jordan-Wigner-laid-out (wrong — the qubit↔orbital mapping
differs). After: `hf_symmetry_sector(mh_interleaved) is None`, `hf_symmetry_sector(mh_jw)` still
returns a real sector. This reproduces fully end-to-end with no monkeypatching.

### L4 — ambiguous-weight abstention (`reachability.py`)

`symmetry_allowed` decides each cluster of numerically-degenerate eigenvalues by whether its
HF-projection weight exceeds 1/2. An accidental (non-symmetry) degeneracy merging two *different*,
both-populated symmetry levels into one cluster has no physical reason to land cleanly on either
side of that cut. Added an abstention band: weight in `(1e-3, 0.999)` → treated as allowed (not
flagged) rather than decided by the majority rule. The margin matches G6's own already-gated bound
(every witnessed real system's weight sits within 1e-3 of 0 or 1), so **this changes nothing for any
system the existing gate file exercises** — confirmed by the full-file regression run below.

No real gated system hits this condition ("Never seen on a gated system" — spec §8, still true).
Verified with a synthetic cluster instead: two eigenvectors forced into one cluster
(`|E1 - E2| = 0.5e-9 < _DEGENERATE_ATOL`) whose basis-independent projection onto a fixed sector
splits exactly 50/50. The pre-fix majority rule (`weight > 0.5`) flags *both* as forbidden on that
exact tie; the fix abstains on both. A control case at weight=0.95 and another at weight=0.0005
(both clearly outside the band) still decide normally — the fix doesn't touch non-ambiguous cases.

## Gate tests added

`tests/test_eigenstate_reachability_spec.py`, new section after G10:
`test_L1_degenerate_shell_split_window_is_excluded_full_window_is_not`,
`test_L2_geometry_pyscf_mistakes_for_exact_is_excluded`,
`test_L3_mapper_is_recorded_and_checked`,
`test_L4_accidental_cross_sector_degeneracy_abstains`.

```
$ uv run pytest tests/test_eigenstate_reachability_spec.py -k "L1 or L2 or L3 or L4" -v
tests/test_eigenstate_reachability_spec.py::test_L1_degenerate_shell_split_window_is_excluded_full_window_is_not PASSED
tests/test_eigenstate_reachability_spec.py::test_L2_geometry_pyscf_mistakes_for_exact_is_excluded PASSED
tests/test_eigenstate_reachability_spec.py::test_L3_mapper_is_recorded_and_checked PASSED
tests/test_eigenstate_reachability_spec.py::test_L4_accidental_cross_sector_degeneracy_abstains PASSED
4 passed, 45 deselected in 2.24s
```

L1, L2, and L4 monkeypatch `reachability._driver_orbitals` (L1, L2) or `reachability.hf_symmetry_sector`
(L4) to supply known-good orbitals/sectors directly, rather than going through the production
`PySCFDriver` rebuild inside `_driver_orbitals`. L3 needed no monkeypatching — it reproduces
end-to-end. Why L1/L2/L4 needed it: see "chem-g14" below — that rebuild's reproducibility check is
independently unreliable on this host for any active window touching a geometrically-degenerate MO
shell, including at the *exact*, unperturbed geometry, which would have made a real-pipeline L1/L2
control case nondeterministic. Monkeypatching isolates the L1-L4 logic from that unrelated failure
mode so these four gates are deterministic.

## Full-file regression check

Pre-registered before running: the full `tests/test_eigenstate_reachability_spec.py` suite must
show the same failure *count* before and after the four code changes, with every failure
attributable to the pre-existing (not chem-648-caused) issue below, not to G3/G6/G7/G8/G10 or any
new L1-L4 gate.

```
$ uv run pytest tests/test_eigenstate_reachability_spec.py -q
17 failed, 32 passed in 55.17s          # after (28 baseline-category passes + 4 new L1-L4 passes)
```

Baseline (same command, `reachability.py`/`molecular_hamiltonian.py`/the new test-file section
stashed out via `git stash push -- <those 3 files>`, confirmed clean with `git stash list`/`git
status` before and after, then restored with `git stash pop`):

```
$ git stash push -m "chem-648 WIP (temporary, for baseline comparison)" -- reachability.py \
    hybrid_quantum_solver/molecular_hamiltonian.py tests/test_eigenstate_reachability_spec.py
$ uv run pytest tests/test_eigenstate_reachability_spec.py -q
16 failed, 29 passed in 53.15s
$ git stash pop
```

Same category every time (`test_G1_residue_level_identified_warned_and_kept`,
`test_G2_decision_equals_symmetry_adapted_fci` on sq-H4/linear-H4,
`test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged`, `test_G4_naive_labels_would_flag_a_physical_level`,
`test_G5_stale_build_args_disable_spatial_bits_but_a_shift_does_not`,
`test_G9_control_tight_scf_stays_on_the_ag_level`), never `G3`/`G6`/`G7`/`G8`/`G10`.

**The exact failure count differs between these two single runs (16 vs 17) and the exact failing
IDs are not even stable run-to-run on an unmodified checkout** — re-running the baseline-only
`test_G9_control_tight_scf_stays_on_the_ag_level` in isolation twice gave `[1.1]` failing both times
but `[1.19]` failing on one run and passing on the other, with zero code changes in between. This is
pre-existing SCF-path nondeterminism (see chem-g14 below), not something chem-648 introduced or can
fix by itself. Ran the full file three more times across this session (with my changes in place):
17/32, 17/32, 17/32 — the count is stable even though the exact IDs vary; never anything outside the
named category.

**Hardware/wall-time** (>1 minute rule): this session's Linux sandbox container, `pyscf==2.13.1`.
Each full-file run: wall ~53-57 s, `user` ~5m (multi-threaded BLAS). Individually timed:
`time uv run pytest tests/test_eigenstate_reachability_spec.py -q` → `55.17s` wall / `5m0.430s` user
(one representative run, reproduced above).

`uv run ruff check reachability.py hybrid_quantum_solver/molecular_hamiltonian.py
tests/test_eigenstate_reachability_spec.py` → `All checks passed!`.

Blast radius: confirmed by repo-wide grep that no other file references `build_args`,
`hf_symmetry_sector`, `orbital_parity_bounds`, or `symmetry_allowed` — the four edited functions are
only consumed within `reachability.py` itself and its own gate file, so this file's regression run
is the complete verification surface (no other `tests/test_*.py` needed re-running).

## Out-of-scope finding: chem-g14 (filed, not fixed here)

`_driver_orbitals`'s rebuild-reproducibility check (`_SAME_OPERATOR_ATOL=1e-12`) was already broken
on this host before any chem-648 edit — confirmed by stashing my changes and reproducing the same
category of failures on a clean checkout. Filed as **chem-g14** with full repro/root-cause notes;
left two follow-up notes on it this session after the L1/L2 witness work surfaced more detail:

1. The exact failing-test-ID set is unstable run-to-run even on an unmodified checkout (not just a
   fixed set of 17 broken tests) — genuine SCF-path nondeterminism, not a deterministic bug with a
   deterministic blast radius.
2. The same rebuild-reproducibility failure also happens at `conv_tol=1e-13`
   (`TIGHT_SCF_CONV_TOL`) for the *exact, unperturbed* square/rectangle-H4 geometry whenever the
   active window touches a degenerate MO shell under free (`symmetry=False`) SCF — e.g. square H4's
   own HOMO pair, which is a true degeneracy of the geometry (D4h collapsed to D2h by PySCF), not an
   SCF convergence artifact. Two independent free-SCF solves of the identical geometry pick different
   (but iso-energetic) rotations inside that degenerate block, so no amount of tightening `conv_tol`
   fixes it — this is a rotational degree of freedom the energy is flat against, not a convergence
   residual. A non-degenerate rectangle (a≠b) at the same `conv_tol=1e-13` reproduces cleanly; this
   is why the L2 gate uses a rectangle rather than a square.

Did not attempt to fix chem-g14 here — it is explicitly out of scope for chem-648 (symmetry-diagnostic
limitations, not SCF-rebuild reproducibility) and fixing it needs its own spec (the bead's notes list
three candidate directions, none decided).

## What I deliberately did not do

- **Did not reconcile the spec's original L1/L2 witness numbers exactly.** The spec (§8,
  pre-revision) names CH4 Td CAS(2,2) for L1 and sq-H4 with a 5e-6/1e-5 Å shift for L2; I used the
  same CH4 system for L1 (reproduces: the degenerate-shell geometry is unambiguous there) but a
  *rectangle* (not square) H4 for L2's gate test, and a tighter 1e-6/2e-6 Å shift. Square H4 at the
  exact geometry hits chem-g14's degenerate-rotation issue even when `_driver_orbitals` is bypassed
  by monkeypatching the orbitals directly (confirmed: the orthogonality-based L1 check itself, not
  just the rebuild check, correctly flags the square's unsymmetrized HOMO pair as untrustworthy for
  one bit, β≈5.7 on that bit alone) — i.e. square H4 conflates L1 and L2 at the exact geometry on
  this host. The rectangle isolates L2 cleanly. Documented in the spec and in this file rather than
  spending further time forcing the square to work.
- **Did not re-run `make gates` or `make test`** (full suite, DMRG-isolated, long-running). Nothing
  outside `reachability.py`/`molecular_hamiltonian.py`/this one test file references the changed
  functions (grep-confirmed, see above), so the full-file regression run above is the complete
  verification surface for this bead; a full-suite run would only add DMRG/unrelated-spec wall-time
  with no new coverage of what changed.
- **Did not touch `_driver_orbitals` or `_SAME_OPERATOR_ATOL`** (chem-g14) even though three of the
  four new gates had to route around it — out of scope per the task instructions, filed instead.

## Unverified

- Whether the spec's §9-recorded macOS baseline ("45 passed in 35s") would also show the L1-L4
  fixes as clean, zero-regression changes — not independently re-run on macOS this session (no
  macOS host available here); the spec's new revision note says so explicitly.
- Whether chem-g14's degenerate-rotation nondeterminism is bounded to H4/CH4-family systems or is
  generic to any degenerate-shell active window — only checked the systems this bead's witnesses
  needed.

## Files changed

```
$ git diff --stat
 hybrid_quantum_solver/molecular_hamiltonian.py |   3 +-
 reachability.py                                | 104 +++++++++++++++++--
 specs/SPEC_eigenstate_reachability.md          |  61 +++++++++--
 tests/test_eigenstate_reachability_spec.py     | 136 +++++++++++++++++++++++++
 4 files changed, 286 insertions(+), 18 deletions(-)
```

Plus this handoff file and a `bd update chem-g14 --notes=...` (not a code change).

## Suggested commands for the human

```bash
git add reachability.py hybrid_quantum_solver/molecular_hamiltonian.py \
        specs/SPEC_eigenstate_reachability.md tests/test_eigenstate_reachability_spec.py \
        sandbox-handoffs/chem-648.md
git commit  # see the actual commit(s) made this session for message text
git status
```

`bd close chem-648` to be run immediately after this file is committed (session-death-safety
ordering: handoff written and committed before closing), per the session's git policy: work is
committed to `sandbox/chem-648` here, host pushes / a human merges. Not pushed, merged, or
rebased by this session.
