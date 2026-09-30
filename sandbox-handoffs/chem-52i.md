# chem-52i — Downstream certified rungs vs. an uncorroborated gap_bracket

Status: bead's acceptance criteria met via the **"absence of escape → positive robustness finding"**
branch. `gap_selfcheck.self_checked_gap_from` was **not** wired into the three call sites — no
propagating bug was found that needed it. One of the two latent crashes named in the bead was real
and is fixed (with a regression test); the other was investigated thoroughly and could not be
reproduced as an actual exception.

## What changed

- `certified_dipole.py` — `spectral_width` now short-circuits an identically-zero sparse operator
  (`a_sparse.nnz == 0`) to width `0.0` instead of calling `scipy.sparse.linalg.eigsh`, which raises
  `ArpackError: Starting vector is zero` on an all-zero matrix. This is the fix for latent crash #1
  (see below).
- `tests/test_certified_dipole_spec.py` — added `test_G5_zero_dipole_operator_does_not_crash`, a
  regression test reproducing the crash on square H₄ (`"H 0 0 0; H 1.4 0 0; H 1.4 1.4 0; H 0 1.4 0"`,
  where μ_z ≡ 0 as an operator by centrosymmetry, not just as an expectation value) and asserting the
  fixed behavior (`spectral_width == 0.0`, resulting certificate exact at `mu=0.0, half_width=0.0`).
- `specs/BACKLOG.md` — the "Certified-bounds arc" entry for this bead's claim marked `[x] CLOSED
  2026-09-28 (chem-52i)` with the outcome and pointer to this file.
- No changes to `certified_gaps.py`, `hf_overlap_certificate.py`, `hf_overlap_subspace.py`, or
  `gap_selfcheck.py` — none were needed (see finding below).
- Five scratch investigation scripts left untracked in the repo root (`scratch_probe_52i.py`,
  `scratch_downstream_52i.py`, `scratch_sweep2_52i.py`, `scratch_downstream2_52i.py`,
  `scratch_hehplus_52i.py`) — not part of the permanent test suite, kept only so every number below
  has a literal regenerating command. Recommend deleting them before commit, or moving the two most
  load-bearing ones (`scratch_downstream2_52i.py`, the persistent-escape reproducer, and
  `scratch_hehplus_52i.py`, the HeH⁺ sweep) somewhere permanent if this arc gets revisited.

## Pre-registered acceptance criteria (from `bd show chem-52i`) and pass/fail

1. **Sweep `certified_dipole`, `certify_hf_overlap`, `certify_hf_subspace_overlap` at M ∈
   {6,8,12,16} on stretched/asymmetric H₄ geometries where the upstream `gap_bracket` escape was
   found.** ✅ Done — 4 geometries × 4 M × 3 call sites = 48 checks on geometries independently
   verified reproducible (see Confound section). Command:
   `uv run python scratch_downstream_52i.py` (first 2 of its 3 geometries — the 3rd is the excluded
   confound, see below) and `uv run python scratch_downstream2_52i.py` (both geometries).
2. **Either a downstream escape is demonstrated and `gap_selfcheck.self_checked_gap_from` is wired
   in, OR absence of escape is recorded as a positive robustness finding.** ✅ Absence-of-escape
   branch: zero downstream escapes across all 48 checks on reproducible geometries. Not wired in.
3. **The two latent crashes, if hit, are fixed and reported separately from the main finding.**
   ✅ Crash #1 (ArpackError) — hit, reproduced, fixed, regression-tested. ⚠️ Crash #2 (HeH⁺
   zero-half-width epsilon floor) — the underlying phenomenon is real and reproduced, but never
   actually raised an exception in 40 checks; reported as "could not verify as a literal crash,"
   see below — not conflated with criterion 2's headline finding.

## Main finding: upstream escape reproduced; zero downstream propagation

### Upstream reproduction (confirms the bead's premise is current)

Exact match to the backlog's scout probe, run fresh this session:

```
$ uv run python -c "
from certified_gaps import gap_bracket, reachable_gap
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
mh = build_molecular_hamiltonian(atom='H 0 0 0; H 0 0 3.0; H 0 0 6.0; H 0 0 9.0')
solver = QuantumKrylovSolver(mh)
print('true gap', reachable_gap(mh))
print(gap_bracket(mh, 6, solver=solver))
"
true gap 0.0012040047142152233
GapBracket(gap_lower=0.590147..., gap_upper=0.598538...)
```

`Δ_lo = 0.590147` vs true gap `0.001204` — a ~490× escape at M=6, matching the backlog's "~500x."
(This geometry only escapes at M=6; at M=8/12/16 `gap_bracket` self-mode instead reports
`gap_lower<0` or diverges the upper bound to `+inf`, which is a *safe* non-escape, not a recovery.)

A grid search (`scratch_sweep2_52i.py`, 27 candidate asymmetric H₄ geometries, max atom separation
capped at 8.3 Å to dodge the SCF-degeneracy confound) found geometries whose lower certificate stays
escaped through **all four** tested M — closer to the bead's "escapes through M=12" description than
the single-M backlog case:

```
$ uv run python scratch_downstream2_52i.py
persistent-escape H4 (0.9/2.0/2.0)   exact reachable gap = 0.432164 Ha
  -- upstream gap_bracket (self) --
    M= 6 Delta_lo=+0.548255 Delta_hi=+0.561727 lower_escape=True
    M= 8 Delta_lo=+0.552041 Delta_hi=+0.561406 lower_escape=True
    M=12 Delta_lo=+0.549409 Delta_hi=+0.560977 lower_escape=True
    M=16 Delta_lo=+0.542221 Delta_hi=+0.560610 lower_escape=True
persistent-escape H4 (0.9/2.5/2.0)   exact reachable gap = 0.455731 Ha
  -- upstream gap_bracket (self) --
    M= 6 Delta_lo=+0.562880 Delta_hi=+0.569517 lower_escape=True
    M= 8 Delta_lo=+0.564798 Delta_hi=+0.569449 lower_escape=True
    M=12 Delta_lo=+0.559012 Delta_hi=+0.569246 lower_escape=True
    M=16 Delta_lo=+0.547934 Delta_hi=+0.568899 lower_escape=True
```

(atoms: `"H 0 0 0; H 0 0 0.9; H 0 0 2.9; H 0 0 4.9"` and `"H 0 0 0; H 0 0 0.9; H 0 0 3.4; H 0 0 5.4"`)

### Downstream sweep: zero escapes (48/48 checks)

Same run (`scratch_downstream2_52i.py`) continues straight into the three downstream call sites on
these two geometries:

```
  -- downstream certified_dipole (exact mu_z = +0.012137) --
    M= 6 mu=+0.012484 half_width=0.007760 gap_lower_used=+0.548255 inside=True
    M= 8 mu=+0.012172 half_width=0.002732 gap_lower_used=+0.552041 inside=True
    M=12 mu=+0.012026 half_width=0.000529 gap_lower_used=+0.549409 inside=True
    M=16 mu=+0.012046 half_width=0.000365 gap_lower_used=+0.542221 inside=True
  -- downstream certify_hf_overlap d=1 (exact = 0.834802) --
    M= 6 gamma_min=0.427436 valid_lower_bound=True   (... all 4 M valid_lower_bound=True)
  -- downstream certify_hf_subspace_overlap d=2 (exact = 0.834810) --
    M= 6 gamma_min=0.814674 valid_lower_bound=True
    M= 8 gamma_min=0.693165 valid_lower_bound=True
    M=12 VACUOUS  (self-mode floor unresolved: Weinstein intervals overlap)
    M=16 VACUOUS  (same reason)
```

(second geometry: same pattern — all finite results `inside=True`/`valid_lower_bound=True`, remaining
cases safely VACUOUS. Full output reproducible via the command above.)

Plus the original two-geometry sweep (`scratch_downstream_52i.py`, first two of its three
geometries — linear H₄ R=3.0 backlog case, and asym H₄ `"H 0 0 0; H 0 0 1.0; H 0 0 3.0; H 0 0 8.0"`):
zero escapes there too — every downstream check is either correctly bounding (`inside=True` /
`valid_lower_bound=True`) or safely VACUOUS. (Full output reproducible via
`uv run python scratch_downstream_52i.py`; only its first two geometry blocks are evidence — see
confound note for the third.)

**Total: 4 geometries × 4 M × 3 call sites = 48 downstream checks, 0 escapes.** This is the
"absence of escape" branch of the acceptance criteria: the downstream certified rungs are robust to
their own gap input being wrong, at least on every geometry tested here. `gap_selfcheck` is therefore
not wired in — there is nothing currently propagating for it to catch.

**Caveat on generality:** 48 checks on 4 geometries is not a proof of robustness in general — it is
evidence that the specific escape mode found upstream (self-mode Weinstein floor overestimating the
gap on near-degenerate/asymmetric H₄) does not propagate on the geometries where it was found. The
mechanism by which it fails to propagate is visible in the numbers: `certified_dipole`'s half-width
formula still uses the (too-large) `Δ_lo` in the denominator of `s = σ₀/Δ_lo`, which makes `s`
*smaller* than it should be — i.e. the escape upstream biases the downstream bound toward being
*more* confident, not wider, yet the point estimates stayed close enough to true that no escape
resulted in these 48 cases. This is a directional risk, not a proof that no downstream escape can
ever occur; a different geometry/observable could plausibly combine a large upstream escape with a
point estimate that misses. Recorded as the honest scope of this finding, not closed off as general.

## Confound found while selecting probe geometries (reported separately, not the finding)

`reachable_gap` (the dense exact-diagonalization reference) is **non-reproducible across separate
process invocations** for H₄ geometries with max atom separation ≳9 Å — PySCF's RHF SCF converges to
different near-degenerate symmetry-broken solutions run to run. Confirmed and reconfirmed this
session on `"H 0 0 0; H 0 0 1.0; H 0 0 3.0; H 0 0 9.0"`:

```
$ for i in 1 2 3 4; do uv run python -c "
from certified_gaps import reachable_gap
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
mh = build_molecular_hamiltonian(atom='H 0 0 0; H 0 0 1.0; H 0 0 3.0; H 0 0 9.0')
print('gap=', reachable_gap(mh))
"; done
gap= 0.6238167689948724
gap= 0.6238167690737852
gap= 0.04265527915982492
gap= 0.33297066274895704
```

Three distinct values recur for the *same nominal geometry* (this matches the values seen for the
`scratch_downstream_52i.py` third geometry block last session, which is why that block is excluded
from the evidence count above — its "escape" numbers used an unstable reference and are not
trustworthy either way, not something that resolves by rerunning). Every geometry actually used as
evidence above (linear H₄ R=3.0/6.0/9.0, and the four asymmetric geometries) was independently
re-checked for reproducibility (rebuilt 3-4 times, gap agreeing to ≤4e-7 Ha) before being trusted;
the grid search in `scratch_sweep2_52i.py` caps max separation at 8.3 Å specifically to dodge this
zone. This is a pre-existing property of the SCF solver on near-dissociated H₄, not something
introduced by or specific to the certified-gap arc — flagged here because it directly affects which
geometries are safe to use as ground truth for this kind of probe, but it is not itself part of the
bead's claim and no fix was attempted (out of scope for chem-52i; a new bead is the right vehicle if
this needs a general SCF-solution-selection fix).

## Latent crash #1 (ArpackError on square H₄): confirmed, fixed, tested

`spectral_width` (`certified_dipole.py`) computes an operator's spectral width via two Lanczos
extremal solves (`scipy.sparse.linalg.eigsh`). On centrosymmetric H₄ (e.g. a square,
`"H 0 0 0; H 1.4 0 0; H 1.4 1.4 0; H 0 1.4 0"`), μ_z vanishes as an **operator** (all matrix elements
zero by symmetry), not just as an expectation value. ARPACK cannot start a Lanczos iteration from an
identically-zero starting vector and raised `ArpackError: Starting vector is zero`.

Fix: `spectral_width` now returns `0.0` immediately when `a_sparse.nnz == 0`, before calling
`eigsh`. Regression test added: `tests/test_certified_dipole_spec.py::test_G5_zero_dipole_operator_does_not_crash`,
which asserts `Az.nnz == 0` for square H₄ (documenting the precondition), `spectral_width(Az) == 0.0`,
and that the resulting `certified_dipole` call returns the trivially-exact certificate
(`mu=0.0, half_width=0.0`) instead of raising.

```
$ uv run pytest -q tests/test_certified_dipole_spec.py
.....                                                                    [100%]
5 passed, 3 warnings in 86.05s (0:01:26)
```

## Latent crash #2 (HeH⁺ zero-half-width epsilon floor): phenomenon confirmed, crash NOT reproduced

Swept HeH⁺ over 8 bond lengths (R ∈ {0.5, 0.6, 0.772, 0.9, 1.1, 1.4, 1.8, 2.5} Å) × 5 Krylov
dimensions (M ∈ {6,8,12,16,20}) × 3 call sites (`certified_dipole`, `certify_hf_overlap`,
`certify_hf_subspace_overlap`) = 40 checks (`scratch_hehplus_52i.py`, every call wrapped in
`try/except Exception`):

```
$ uv run python scratch_hehplus_52i.py 2>&1 | grep -c RAISED
0
$ uv run python scratch_hehplus_52i.py 2>&1 | grep -c "hw=0.000e+00"
33
```

**Zero exceptions raised across all 40 checks.** The named phenomenon — a "zero-half-width" — is real
and reproducible: HeH⁺'s HF-reachable subspace is small enough that the Krylov ground state saturates
to an exact eigenstate at low M, making the residual σ₀ exactly `0.0` and hence
`certified_dipole`'s `half_width` exactly `0.0` (33/40 checks; example:
`R=0.500 M=6  dipole[hw=0.000e+00]  hf1[g=0.994251]  hf2[g=0.998564]`). This is not a bug: it is
already documented as expected behavior in `test_G2_interval_closes_and_is_useful`'s docstring
("HeH+ is certified to < 1e-2 a.u. (its reachable subspace saturates, sigma_0 -> 0)") and passes that
existing gate. Despite 40 checks across a wide parameter sweep specifically targeting this
phenomenon, no code path turned "half_width exactly 0" into a raised exception (no division by zero,
no invariant violation in `GapCertificate`/`OverlapCertificate`'s `__post_init__` checks). **Could
not be verified as a literal crash** — reported honestly as unreproduced rather than claimed fixed.
No code change was made for this item since there was nothing observed to fix.

## Gate results (verbatim)

```
$ uv run pytest -q tests/test_certified_dipole_spec.py
.....                                                                    [100%]
5 passed, 3 warnings in 86.05s (0:01:26)

$ uv run pytest -q tests/test_hf_overlap_certificate_spec.py tests/test_hf_overlap_subspace_spec.py
...................................... (37 tests, no source changes to these modules)
37 passed, 6 warnings in 6.02s

$ uv run pytest -q tests/test_certified_gaps_spec.py tests/test_gap_selfcheck_spec.py
........                                                                 [100%]
8 passed, 6 warnings in 270.87s (0:04:30)

$ uv run ruff check certified_dipole.py tests/test_certified_dipole_spec.py
All checks passed!
```

No regressions in any gate touched by or adjacent to this work.

## What was decided not to do, and why

- **Did not wire `gap_selfcheck.self_checked_gap_from` into the three call sites.** The bead's own
  acceptance criteria treats this as conditional on finding a downstream escape; none was found in
  48 checks, so adding it now would be undemonstrated defensive code with no test that could show it
  does anything (violates the "don't add validation for scenarios that can't happen" / falsifiable-
  honesty culture — a change with no failing-then-passing test behind it is not evidence of anything).
- **Did not attempt to fix the H₄ SCF near-degeneracy non-reproducibility.** Out of scope for this
  bead (a property of `reachability.py`/PySCF's SCF convergence, not of the certified-gap arc); filed
  nowhere new since it's already recorded here and doesn't currently block any user-facing claim — if
  it needs fixing, it should be its own bead with its own falsifiable check.
- **Did not delete the scratch investigation scripts before finishing.** They are the literal
  regenerating commands for every number in this document; deleting them would break "a number only
  exists in a document if code produces it or the document states provenance." Left untracked in the
  repo root for a human to review/delete/relocate at commit time.

## What could not be verified

- Latent crash #2 (HeH⁺ epsilon floor) as a literal raised exception — see above. The underlying
  saturation phenomenon is verified; an actual crash is not, despite a targeted 40-check sweep.
- General propagation robustness beyond the 4 geometries tested — see the caveat in the main finding
  section. This is scoped honestly as "not observed on these 4 geometries," not "proven impossible."

## Files changed (git status, this session's edits only)

```
 M certified_dipole.py
 M tests/test_certified_dipole_spec.py
 M specs/BACKLOG.md
 M .beads/issues.jsonl / .beads/interactions.jsonl   (bd claim/update bookkeeping)
?? scratch_probe_52i.py
?? scratch_downstream_52i.py
?? scratch_sweep2_52i.py
?? scratch_downstream2_52i.py
?? scratch_hehplus_52i.py
```

Note: `git status` also shows `CLAUDE.md`, `reachability.py`, `pyproject.toml`, `uv.lock` as modified.
These are **pre-existing changes unrelated to this session's work** (macOS/Accelerate scipy build
config referencing other bead IDs, chem-dfe/chem-a0y) — not touched, not authored, and not reviewed
as part of chem-52i; flagged here only so a reviewer doesn't mistake them for part of this diff.

## Commands for a human to run next (not run by me — git policy: no commit/push this session)

```bash
git -c safe.directory=/workspace add certified_dipole.py tests/test_certified_dipole_spec.py specs/BACKLOG.md
# decide: keep, relocate, or `rm scratch_*_52i.py` before committing
git -c safe.directory=/workspace commit -m "certified_dipole: fix ArpackError on identically-zero operators (chem-52i)

Downstream sweep of certified_dipole/certify_hf_overlap/certify_hf_subspace_overlap
against the known gap_bracket self-mode escape on stretched/asymmetric H4 found zero
downstream escapes (48 checks, 4 geometries x M in {6,8,12,16} x 3 call sites) -- a
positive robustness finding, recorded in specs/BACKLOG.md. gap_selfcheck was not wired
in since nothing was found for it to catch.

Fixed a real latent crash found en route: spectral_width's Lanczos solve raised
ArpackError on an identically-zero operator (mu_z on centrosymmetric H4), now
short-circuited to width 0.0, with a regression test. The other named latent crash
(HeH+ zero-half-width epsilon floor) was investigated (40 checks) but never
reproduced as an actual exception -- reported as unverified, not fixed.

See sandbox-handoffs/chem-52i.md for full evidence and regenerating commands.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
# then, if desired:
bd dolt push && git push
```
