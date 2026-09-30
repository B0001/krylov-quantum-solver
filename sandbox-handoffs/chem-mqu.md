# Handoff: chem-mqu — the overlap certificate's noise sibling

## What I changed

New files (nothing else touched except the BACKLOG entry for this hypothesis):

- `hf_overlap_noise.py` — `exact_margin(mh, m, solver=None)` (zero-noise reference point: λ_u, r,
  δ, margin, exact overlap) and `hf_overlap_noise_coverage(mh, m, shots, z=2.0, trials=6000,
  seed=0, solver=None)` (Monte-Carlo coverage under shot noise that enters ONLY through β).
- `specs/SPEC_hf_overlap_noise.md` — the spec, with measured numbers in §5 (not `data/`, which is
  gitignored).
- `tests/test_hf_overlap_noise_spec.py` — gates G0–G7.
- `specs/BACKLOG.md` — closed the "overlap certificate has no noise sibling" entry (line ~250) with
  the finding, `[x] CLOSED (chem-mqu)`.

## Command that regenerates every number below

```
uv run python hf_overlap_noise.py
```
(full printed table, all 48 grid cells; ~7s wall time on this container)

```
uv run python -m pytest -q tests/test_hf_overlap_noise_spec.py -v
```
→ **10 passed** (G0–G7, some parametrized over shots). Ran standalone, not inside `make gates`;
this module never imports `pyscf`/`block2` in a way that needs process isolation (it's pure
PySCF+qiskit, same isolation class as the rest of `certified_*` — no DMRG).

```
uv run ruff check hf_overlap_noise.py tests/test_hf_overlap_noise_spec.py
```
→ `All checks passed!`

Regression check on the sibling module it depends on:
```
uv run python -m pytest -q tests/test_hf_overlap_certificate_spec.py
```
→ 23 passed (unchanged, confirms I didn't touch anything it relies on).

## Pre-registered criteria (from the bead) and outcome

| Criterion (bead's own words) | Outcome |
|---|---|
| ≥6000 seeded trials, shots∈{1e4,1e5,1e6}, z∈{0,1,2,3}, 4 systems spanning wide→vacuous margin | Done exactly as specified: M=8 self-mode, trials=6000, seed=0, full 4×3×4=48-cell grid. |
| Dies if z≤2 restores coverage≥0.9 **everywhere** | **Does not die.** H₄-moderate and H₄-thin top out at 0.8308 across the *entire* z≤2 sub-grid (all 3 shot counts) — not a near miss. See G3. |
| Determinant-only scope caveat stated in the writeup | Stated in §2/§8 R1 of the spec, in the BACKLOG close-out, and in `hf_overlap_noise.py`'s module docstring: none of this carries to the chained-Ritz bound (`krylov_refine.py`), whose v needs its own shot-noisy ⟨v\|H\|v⟩, ⟨v\|H²\|v⟩. |
| Coverage-vs-margin curve reported | §3 table (4 systems, margin from +1.419 to −0.255) + G4/G5 swing numbers at shots=1e4. |
| Bifurcation-vs-smooth-degradation question answered | Answered: not a hard binary cliff, but a *sharp sensitivity spike* — the z=0→3 coverage swing is 3–4× larger at thin margin than at wide/moderate margin (G4), and — the unpredicted part — **the swing changes sign** between regimes (G5): inflation helps at wide margin, hurts at moderate/thin margin. This is a sharper, more specific finding than the bead's original "smooth vs bifurcating" framing anticipated, and it directly explains *why* z≤2 can't just be cranked up everywhere the way it is in `certified_noise`/`gap_selfcheck_noise`. |

## The numbers (M=8, self-mode, trials=6000, seed=0)

Four systems, self-mode M=8 (same default as `hf_overlap_certificate`), spanning the margin
δ_exact − r as the bead specifies:

| system | margin | r/δ | exact overlap |
|---|---|---|---|
| H₂ eq (0.74 Å) | +1.419 (wide) | 0.113 | 0.9936 |
| H₄ chain (0.9/1.8/2.7 Å) | +0.151 (moderate) | 0.649 | 0.9769 |
| H₄ chain (1.0/2.0/3.0 Å) | +0.053 (thin) | 0.842 | 0.9677 |
| square H₄ a=1.2 Å | −0.255 (already vacuous, zero noise) | n/a | 0.0000 |

Coverage swing z: 0→3 at shots=1e4 (the sharpest cell): H₂ **+0.111**, H₄-moderate **−0.138**,
H₄-thin **−0.451**, square H₄ 0 (pinned). Full 48-cell table is reproduced by the command above and
is not re-pasted here — it's in the spec (§5) and the module's own `__main__` output.

## What I decided not to do, and why

- **Did not extend to the chained-Ritz bound.** The bead explicitly scopes this out ("a chained-Ritz
  state puts r back on the noise budget — state that scope boundary explicitly"). Doing it properly
  needs a two-noisy-input model (λ_v, r_v *and* β all noisy), which is a different, harder spec —
  filed nowhere yet because the bead didn't ask for it; if wanted, it's a natural chem-* follow-up
  bead, not squeezed into this one.
- **Did not sweep a continuum of margins.** The bead names exactly the endpoints (H₂ 0.74 Å →
  square H₄ 1.2 Å); I picked two named, already-existing intermediate geometries from
  `hf_overlap_certificate.py`'s and `hf_overlap_subspace.py`'s own `__main__` demos (H₄ chain
  0.9/1.8/2.7 Å and 1.0/2.0/3.0 Å) rather than inventing new ones, so every system here is one this
  repo already treats as a reference point elsewhere.
- **Did not test oracle-mode β.** Self mode is the production path the certificate exists for
  (`hf_overlap_certificate.py`'s own framing); oracle mode has no premise boundary to stress-test
  and wasn't asked for.
- **Did not try a symmetric (two-sided) pad on β as an alternative to the one-sided
  `gap_selfcheck_noise._pad` convention.** One-sided follows directly from β being a lower bound
  (same reasoning `gap_selfcheck_noise` already used for `gap_lower`); flagged as R3 in the spec as
  an alternative not tried, in case a future reader wants to check whether it changes the *magnitude*
  of the G5 sign reversal (it shouldn't change the *direction*, which comes from which failure mode
  dominates, but I did not verify that empirically).

## What I could not verify

- I did not independently re-derive the sign-reversal mechanism from the Davis–Kahan formula by
  hand beyond the qualitative argument in the spec/module docstrings (pushing β down shrinks δ,
  which can only ever decrease γ_min for fixed r — so more overclaim-fixing at wide margin and more
  vacuousness at thin margin follows from the same monotonicity, verified empirically by G1's
  frac_vacuous/frac_invalid monotonicity check, but not proved algebraically here beyond that).
- Single seed (seed=0) throughout, per the bead's own trial-count spec (≥6000 trials substitutes for
  multi-seed averaging here, as in every other noise spec in this repo) — I did not check a second
  seed for reproducibility of the exact swing ratios (3–4×), only that the *qualitative* pattern
  (G3–G5) is not a knife-edge: the margins between the pass/fail thresholds in the gates and the
  measured values are wide (e.g. G3's 0.83 vs the 0.9 bar; G4's 3.9×/3.3× vs the 2×/3× bars), so I'm
  confident this isn't seed-dependent, but that confidence is not itself a test artifact.

## Bead status

`bd close chem-mqu` — acceptance criteria met, gates green, evidence written up above and in the
spec before closing.

## Git state (left dirty, per instructions — do not commit/push)

```
 M specs/BACKLOG.md
?? hf_overlap_noise.py
?? specs/SPEC_hf_overlap_noise.md
?? tests/test_hf_overlap_noise_spec.py
```

Suggested commands for a human to run (not run here):
```
git add hf_overlap_noise.py specs/SPEC_hf_overlap_noise.md tests/test_hf_overlap_noise_spec.py specs/BACKLOG.md
git commit -m "hf_overlap_noise: the overlap certificate's noise sibling -- coverage bifurcates, inflation can reverse sign (chem-mqu)"
```
