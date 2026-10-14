# chem-3fd handoff

## Disposition: closed as not-a-bug (investigation complete, no code change)

Both symptoms chem-3fd describes are fully explained by issues already filed elsewhere, or by a
genuine, provable mathematical limit of the current symmetry bit set for this geometry/basis
class — not by a defect in `orbital_parity_bounds`/`hf_symmetry_sector`'s projector or budget
math. No fix is possible within this bead's stated scope ("the budget selection, the bit-to-
operation mapping, or the projector math") because none of the three is at fault.

## Bead content (summary, see `bd show chem-3fd` for the full verbatim text)

1. **test_G2** (square H4, every bond length × `conv_tol`): after chem-g14's fix,
   `orbital_parity_bounds`/`hf_symmetry_sector` succeed but over-include — e.g. a=1.35,
   `conv_tol`=1e-9 gives 20 levels in the "Ag sector" vs. the independent symmetry-adapted FCI's
   10.
2. **test_G4[1.2]**: square H4's free RHF at a=1.2 is claimed to be genuinely symmetry-broken
   (0.076 Ha lower); the test expects every bit to report β ≥ 1.0 ("excluded"), but one bit
   computes β ≈ 0.0008 ("trustworthy") instead.

## Pre-registered criteria and outcome

| Criterion | Outcome |
|---|---|
| Identify whether the over-inclusion (symptom 1) is a bug in the projector/budget math | **Not a bug** — proved it is a provable, exact group-theoretic degeneracy of `_BIT_OPERATION["D2h"]`'s three names for this geometry/basis class (below). |
| Identify whether the β≈0.0008 case (symptom 2) is a detection failure | **Not reproducible as described** — a=1.2's free RHF is not symmetry-broken on this host at all (20/20 trials, see below), so there is nothing for any bit to detect. |
| If a fixable bug survives, implement it + add/update a gate test | N/A — no fixable bug found. |
| `make gates` / relevant pytest passes after any change | No production code changed; pre-existing failures persist (see "Current gate state" below) and are chem-g14's/chem-7cx's responsibility, not introduced or fixable here. |

## Why chem-g14 has to be bypassed to see chem-3fd's actual symptom

chem-g14 (`reachability._driver_orbitals`'s rebuild-reproducibility check) is still unmerged on
this branch (cut from `main`). Without its fix, `orbital_parity_bounds` returns `None`
unconditionally here, so the *currently observable* failures are chem-g14's documented shape
(36-vs-N "Ag sector" counts, `TypeError` on G4), not the 20-vs-10 / β≈0.0008 shape chem-3fd
describes (which is specifically the *post*-chem-g14-fix symptom, per the bead's own text).
Verbatim, this branch today:

```
$ uv run pytest "tests/test_eigenstate_reachability_spec.py::test_G2_decision_equals_symmetry_adapted_fci" -v
...
FAILED [...][H 0 0 0; H 1.35 0 0; H 1.35 1.35 0; H 0 1.35 0-1e-09]
  AssertionError: (36, 10)
9 failed, 1 passed in 7.56s

$ uv run pytest tests/test_eigenstate_reachability_spec.py -k "G4" -v
test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged[1.2] FAILED  (TypeError: 'NoneType' object is not iterable)
test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged[1.4] FAILED  (same)
test_G4_naive_labels_would_flag_a_physical_level FAILED               (same)
test_G4_lih_c2_bit_is_exact PASSED
3 failed, 1 passed in 11.77s
```

To investigate chem-3fd's actual claim, I monkeypatched only `reachability._driver_orbitals` to
call `PySCFDriver` directly and skip the `_SAME_OPERATOR_ATOL=1e-12` reproducibility gate — the
same isolation technique chem-648's own L1/L2/L4 gates use (see `sandbox-handoffs/chem-648.md`).
Everything downstream (`orbital_parity_bounds`, `hf_symmetry_sector`, `scf_symmetry_status`) is
the real, unmodified production code. This matches the shape the bead itself predicts (20 vs 10
at a=1.35), which is independent confirmation the bypass is faithful.

## Symptom 1 (test_G2 over-inclusion): root cause found, it is exact and structural

sq-H₄/sto-3g is planar (every atom at z = 0) with an all-s-type basis. On this geometry+basis,
σ_xy (the xy-plane reflection) acts as the **identity** on every AO — an s orbital centered on a
z = 0 atom is exactly invariant under z → −z. D2h's group relation `i = σ_xy · C2z` then forces
`rep(i) = rep(σ_xy) · rep(C2z) = rep(C2z)` **exactly**. `_BIT_OPERATION["D2h"] = (C2y, C2z, i)`
therefore carries only 2 independent bits for this system, not 3 — confirmed directly: bit1's and
bit2's `R_b` matrices are byte-identical at every bond length tested.

```
$ uv run python /tmp/debug_anomaly2.py   # a=1.00, direct inspection of r (= R_b) per bit
bit 1 (C2z): r = [[1,0,-0,-0],[0,-1,-0,-0],[-0,-0,-1,-0],[-0,-0,-0,1]]
bit 2 (i):   r = [[1,0,-0,-0],[0,-1,-0,-0],[-0,-0,-1,-0],[-0,-0,-0,1]]   # identical to bit 1
```

Consequence for `hf_symmetry_sector`'s greedy budget: whenever the C2y-type bit (bit 0) is
excluded for *any* reason, only one effective constraint beyond (N_α, N_β) remains from the
other two (redundant) bits — mathematically insufficient to cut 20 candidate levels down to the
true 10-level Ag sector. This is **not** a defect in the budget or projector math; both are
correct given the bits available. D2h genuinely does not have a third independent generator for
a planar, s-orbital-only geometry.

a=1.35, `conv_tol`=1e-9, 5 fresh-process trials (bypass in place), confirming this is deterministic
once the RHF branch is fixed (not noise):

```
trial=0 len(in_sector)=20 len(ref)=10  betas=['5.6569','0.0000','0.0000']  broken=True de=-8.24e-02
trial=1 len(in_sector)=20 len(ref)=10  betas=['5.6567','0.0000','0.0000']  broken=True de=-8.24e-02
trial=2 len(in_sector)=20 len(ref)=10  betas=['5.6568','0.0000','0.0000']  broken=True de=-8.24e-02
trial=3 len(in_sector)=20 len(ref)=10  betas=['5.6566','0.0000','0.0000']  broken=True de=-8.24e-02
trial=4 len(in_sector)=20 len(ref)=10  betas=['5.6568','0.0000','0.0000']  broken=True de=-8.24e-02
```

bit0 is correctly excluded (β ≈ 5.657 ≥ 1.0, genuine symmetry breaking, confirmed stable); bit1/
bit2 are the redundant pair above. 20 is the correct, inevitable output of the current mechanism
given this bit set — there is no "fix" that keeps `_BIT_OPERATION["D2h"]`'s three names and also
resolves this.

A second, related instability makes this worse at some bond lengths: sq-H₄ is a **perfect
square** at every tested `a` (true point group D4h, non-abelian; PySCF only ever detects the D2h
subgroup it's abelian). Some orbital pairs are then exactly degenerate under the untracked
symmetry, and the free (`symmetry=False`) SCF's arbitrary tie-breaking within that 2D degenerate
space is itself unstable — even for the identical, fully-converged, lowest (non-broken) energy.
At a=1.00 (15 fresh-process trials, same code, same geometry, same `conv_tol`):

```
trial=0  broken=False de=-7.1e-14  betas=['5.6569','0.0000','0.0000']
trial=2  broken=False de=-3.8e-13  betas=['0.0000','0.0000','0.0000']
trial=3  broken=True  de=-6.6e-02  betas=['5.6568','0.0000','0.0000']
trial=4  broken=True  de=-6.6e-02  betas=['0.0000','0.0000','0.0000']
...
```

The C2y-type bit's β flips between ≈0 and ≈5.66 with **no change in energy or broken-status** —
confirming this is the arbitrary-degenerate-rotation effect already documented and accepted for
LiH's π pair in `test_G4_lih_c2_bit_is_exact`'s docstring (PR #62 dropped an analogous assertion
for exactly this reason), now also present for square H4. Both mechanisms are documented together
in `specs/SPEC_eigenstate_reachability.md` §8 R2 (extended this session, see "What changed"
below).

## Symptom 2 (test_G4[1.2] β≈0.0008): premise not reproducible on this host

`scf_symmetry_status`'s docstring says square H4's broken solution is lower by "~0.08 Ha" in
general, not specifically at a=1.2; the bead over-attributed that generic docstring number to
a=1.2 in particular. Directly measured, repeated, fresh-process trials at a=1.2 (the test's own
default `conv_tol`), 20 trials across two separate scripts:

```
$ uv run python /tmp/debug_a14.py   # a=1.20, default conv_tol, 10 trials
trial=0 broken=False de=-1.776e-15 betas=['0.0002','0.0000','0.0000']
trial=1 broken=False de=-1.688e-14 betas=['0.0006','0.0000','0.0000']
trial=2 broken=False de=-8.882e-16 betas=['0.0002','0.0000','0.0000']
trial=3 broken=False de=2.665e-15  betas=['0.0001','0.0000','0.0000']
trial=4 broken=False de=-6.217e-14 betas=['0.0008','0.0000','0.0000']   # matches the bead's own cited β
... (10/10 trials: always broken=False, de within 1e-10 of zero)
```

a=1.2's free RHF **never** lands on a symmetry-broken solution on this host (de always ~1e-10 to
1e-16, far inside `SYMMETRY_BREAK_TOL`=1e-6) — contradicting the test's premise and the
`scf_symmetry_status` docstring's generic witness as applied to this specific bond length. The
β≈0.0008 the bead measured (and which I reproduced, trial 4 above, on a genuinely non-broken
reference!) is exactly what `orbital_parity_bounds` *should* report for a non-broken, trustworthy
MO set — small, because there is nothing to exclude. There is no detection failure here: a
broken-symmetry reference to detect was never produced in the first place, on this platform.

For contrast, a=1.4 (the test's other parametrize value) is stable broken=True across 10/10
trials with the C2y-type bit correctly at β≈5.657 every time — the bead's own text already notes
"a=1.4 ... DOES pass", consistent with this.

```
$ uv run python /tmp/debug_a14.py   # a=1.40, default conv_tol, 10 trials
trial=0..9: broken=True de=-8.404e-02  betas=['5.6567'..'5.6569', '0.0000', '0.0000']  (10/10 stable)
```

This is the same root-cause family chem-7cx already documents (the free RHF solve is not pinned
to a consistent branch) — just a much larger-magnitude manifestation (a full symmetric ↔ broken
basin flip at a=1.05/1.10/1.19/1.20, not chem-7cx's previously-documented 7th-digit energy drift).
Added this evidence to chem-7cx's notes (`bd update chem-7cx --notes=...`, done this session).

## What changed

- `specs/SPEC_eigenstate_reachability.md` §8, the existing **R2** limitation entry: extended with
  the exact mechanism above (σ_xy = identity ⟹ rep(i) = rep(C2z), the degenerate-rotation
  instability, and why this is out of scope to fix without a basis/geometry-dependent bit set
  beyond `_BIT_OPERATION["D2h"]`'s three names). This is the "revise the spec and record why" path
  CLAUDE.md's SDD section calls for when a gate turns out to be unsatisfiable by design — it is
  not a new claim, it names the mechanism behind a limitation the spec already flagged as open.
- No changes to `reachability.py`, `hybrid_quantum_solver/molecular_hamiltonian.py`, or the test
  file — no fixable bug found in any of them.

## Out-of-scope findings filed separately (not fixed here)

- **chem-7cx** (already filed, updated with notes): a=1.05/1.10/1.19/1.20's RHF branch selection
  (broken vs. symmetric) is platform-dependent, with a much larger blast radius than chem-7cx's
  original 7th-digit-drift evidence — includes full basin flips that change `test_G2`/`test_G4`'s
  behavior, not just `test_G1`/`test_G9_control`'s.
- **chem-0m8** (new): `test_G2`/`test_G4` in `test_eigenstate_reachability_spec.py` hardcode bond
  lengths whose broken/symmetric branch chem-7cx shows is platform-dependent, and assert a
  specific numeric outcome without checking which branch was actually reached. This pollutes
  `make gates`'s signal on non-macOS hosts for reasons that have nothing to do with code
  correctness. Filed as its own task since fixing it means touching the gate file's semantics,
  which is a different and riskier change than the investigation chem-3fd asked for.

## Current gate state (unmodified code, this branch, this host)

```
$ uv run pytest tests/test_eigenstate_reachability_spec.py -q
18 failed, 31 passed in 54.54s
```

All 18 failures are `test_G1`/`test_G2`(sq-H4 or linear-H4)/`test_G4`/`test_G5`/`test_G9_control`
cases — the chem-g14 category (confirmed: `orbital_parity_bounds` returns `None` for every case,
`TypeError` on G4) plus chem-7cx's `test_G9_control[1.1]`/`[1.19]` nondeterminism. No
`test_G3`/`G6`/`G7`/`G8`/`G10`/`test_L*` case fails. This is unchanged by this session (no
production code touched); recorded here for the next worker's baseline.

```
$ uv run ruff check .
All checks passed!
```

## What could not be verified

- Whether chem-g14's actual (unmerged) fix, applied for real rather than simulated via the
  `_driver_orbitals` bypass, reproduces exactly the same 20-vs-10 and β-pattern numbers reported
  here. The bypass is the same isolation technique chem-648's own L1/L2/L4 gates use and the
  resulting shape (20 vs 10 at a=1.35) matches what the bead's own text already predicted for
  after chem-g14 lands, which is strong but not direct confirmation.
- Whether the original worker who filed chem-3fd (via chem-g14's handoff) was on a platform
  where a=1.2 *does* reliably land on the broken branch — I have no access to that host. Given
  chem-7cx's already-established cross-run nondeterminism and my own 20/20-trial result here, I
  believe the most likely explanation is platform/BLAS/LAPACK divergence (chem-7cx's scope), not a
  difference in `reachability.py` between hosts, but this is inference, not direct comparison.
- Whether `_BIT_OPERATION`'s D2h generator choice could be replaced by a different triple of D2h
  elements that avoids the σ_xy-identity collapse for *this* system while not introducing an
  analogous collapse for some other planar/axial geometry already gated elsewhere (e.g. linear
  H4, which is also D2h). Did not attempt this — it is a redesign of `_BIT_OPERATION`, well beyond
  "investigate whether ... is at fault," and would need its own spec revision and gate coverage if
  ever pursued.

## Files changed

```
$ git diff --stat
 specs/SPEC_eigenstate_reachability.md | 18 +++++++++++++++++-
 1 file changed, 17 insertions(+), 1 deletion(-)
```

Plus this handoff file, `bd update chem-7cx --notes=...`, and `bd create` for chem-0m8 (not code
changes).

## Suggested commands for the human

```bash
git add specs/SPEC_eigenstate_reachability.md sandbox-handoffs/chem-3fd.md .beads/issues.jsonl
git commit  # see the actual commit(s) made this session for message text
git status
```

`bd close chem-3fd --reason=...` to be run immediately after this file is committed
(session-death-safety ordering), per the session's git policy: work is committed to
`sandbox/chem-3fd` here, host pushes / a human merges. Not pushed, merged, or rebased by this
session.
