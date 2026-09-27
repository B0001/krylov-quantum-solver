# chem-g1i handoff — SPEC_nbn_low_spin.md corrected to the real geometry

## Bead

**chem-g1i**: "SPEC_nbn_low_spin.md Tables 1-3 and G5 threshold need correction: reconstruction
used wrong bond length." Filed by chem-bbi after finding that the committed `data/nbn_scf.chk`
encodes a Nb-N separation of 3.7255 Å, not the 2.25 Å molecule `SPEC_nbn_low_spin.md` (chem-dc7)
had been gated against.

## What changed

- `specs/SPEC_nbn_low_spin.md` — full rewrite. New §0 "Correction" section up front explaining the
  root cause and giving the headline real-geometry numbers. Tables 1-3 rewritten with real,
  vendored-chk numbers (plus a new S=2 row in Table 1 that chem-bbi's spec had explicitly left
  unmeasured); the old chem-dc7 tables are kept, not deleted, folded into a `<details>` block
  clearly marked "WRONG-BOND-LENGTH RECONSTRUCTION". §2, §4, §6-§8 updated to match.
- `tests/test_nbn_low_spin_spec.py` — `test_G4_low_spin_is_hard_preregistered` replaced with
  `test_G4_low_spin_cheap_D_is_pre_registered_inconclusive` (the honest finding: on the real
  geometry neither the pre-registered KILL nor CONFIRM condition fires at cheap D — this is a gap
  in the original pre-registered design, not a threshold nudge in either test's favor).
  `test_G5_committed_sector_is_not_the_cas_ground`'s margin recalibrated from 5 mHa to 1 mHa
  (measured real gap: 1.538 mHa DMRG / 1.537 mHa exact FCI), docstrings updated with provenance.
  Module docstring updated to note the correction.
- Filed **chem-owr** (P3, out of scope for this bead): the `nbn_low_spin.py scan` reconstruction
  pipeline does not reproduce the committed chk's energy even at the real bond length (3.7255 Å) —
  a 4.84 mHa discrepancy, undiagnosed. See its description for what was ruled out.
- No changes to `nbn_dmrg_reference.py`, `nbn_low_spin.py`, or any shared/library code — this bead
  only invoked existing code via CLI/imports to generate new provenance-recorded data points.

## The key scientific finding (reversal, not just renumbering)

chem-dc7's original spec CONFIRMed a pre-registered "hard multireference benchmark" framing for
the (7,7) singlet sector, based on the wrong-bond-length reconstruction. On the real geometry:

- At cheap D (≤300, the CI-gate schedules): **neither** the pre-registered KILL nor CONFIRM
  condition fires. Measured: dw(300) ≈ 1.05e-7 (schedule A′) / 1.04e-7 (B′) — just outside the
  1e-7 KILL floor; per-D spread ≈ 3.1e-5 Ha (A′) / 6.1e-5 Ha (B′) — also just outside the 1e-5 Ha
  KILL floor; |E_A′ − E_B′| ≈ 4.18e-6 Ha — far under the 1e-4 Ha CONFIRM floor. This is a genuine,
  reproduced gap in what the pre-registered design anticipated, not a result either test could be
  bent to declare.
- At headline D (≤1200): the sector **converges** (dw → ~3e-10, two independent schedules agree to
  0.3 nHa). **The (7,7) sector is soft on the real geometry** — the hard-benchmark framing does
  not survive. It was an artifact of the wrong-bond-length molecule.
- The qualitative sector-ordering conclusion (S=1 is the true CAS ground, not the committed S=3)
  **does survive**, but the gap shrinks by ~10×: 1.538 mHa (real) vs 15.1 mHa (wrong geometry).

## Numbers, with regenerating commands and provenance

All new numbers below are in `results/nbn_low_spin/runs.jsonl` (append-only, tracked; timestamps
2026-09-26T15:40:55 through 2026-09-26T16:23:52, host `9bb17f03582d`, x86_64, 8 cores, 16 GB — this
container). The one exception (S=1 exact FCI) was recorded earlier by chem-bbi
(2026-09-26T15:07:42, host `84689dcda54e`) and is reused here unchanged, not re-run.

| # | number | value | wall | command |
|---|---|---|---|---|
| 1 | E(S=2), exact FCI, nelec=(9,5) | −110.04698592997889 Ha, ⟨S²⟩=6.0000045405141735 | 205.4 s / 8 cores | `uv run python nbn_low_spin.py fci --nelec 9 5` |
| 2 | E(S=3), exact FCI, nelec=(10,4) | −110.04602841723886 Ha, ⟨S²⟩=12.000000005607337 | 13.2 s / 8 cores | `uv run python nbn_low_spin.py fci --nelec 10 4` |
| 3 | DMRG headline A, nelec=(7,7), D=400/800/1200 (perD) | −110.04249952166984 Ha, dw(1200)=3.06e-10 | 1029.4 s / 4 threads | via `nbn_dmrg_reference.run_schedule("A", n_threads=4, nelec=(7,7))` |
| 4 | DMRG headline B, nelec=(7,7), D=300/600/1200 (ramp) | −110.04249952196103 Ha, dw(1200)=2.73e-10 | 530.0 s / 4 threads | `run_schedule("B", n_threads=4, nelec=(7,7))` |
| 5 | DMRG cheap A′, nelec=(7,7), D=100/200/300 (perD) | −110.04250013542114 Ha, dw(300)=1.053e-7 | 238.1 s / 2 threads | `run_schedule("A'", n_threads=2, nelec=(7,7))` |
| 6 | DMRG cheap B′, nelec=(7,7), D=80/160/300 (ramp) | −110.04249595430977 Ha, dw(300)=1.036e-7 | 124.3 s / 2 threads | `run_schedule("B'", n_threads=2, nelec=(7,7))` |
| 7 | scan(d=3.7255 Å): E(10,4)=−110.0411891543836, E(9,5)=−110.0454690725500 | see chem-owr | 225.2 s / 8 cores | `uv run python nbn_low_spin.py scan --d 3.7255` |

Reused, not re-run: E(S=1) exact FCI = −110.04756504307636 Ha, ⟨S²⟩=2.0000095761529892
(1010.6 s / 8 cores, chem-bbi, 2026-09-26T15:07:42) and the historical DMRG A′ nelec=(10,4)
number = −110.04602823118734 Ha (chem-bbi).

All seven new numbers were arithmetic-checked against the spec's tables before writing (gaps,
sums, orderings) — see the Δ columns in `SPEC_nbn_low_spin.md` Table 1/2.

## Pre-registered criteria and pass/fail (specs/BACKLOG.md, unmodified wording)

- KILL: dw(D=300) < 1e-7 **and** per-D spread < 1e-5 Ha → **did not fire** (dw(300) ≈1.05e-7 is
  just above the floor).
- CONFIRM: dw > 1e-5 **or** |E_A′−E_B′| > 1e-4 Ha → **did not fire** either (dw(300) is 2 orders
  below 1e-5; the A′/B′ disagreement is 2 orders below 1e-4).
- Neither pre-registered condition applies at the pre-registered checkpoint (D=300). This is
  reported as an honest gap, not resolved by inventing a new threshold — the definitive answer
  (soft, from headline D=1200 runs) is recorded in the spec but the pre-registered test itself now
  only asserts the reproduced inconclusive result plus regression-pin bounds on the measured
  values (see `test_G4_low_spin_cheap_D_is_pre_registered_inconclusive`).
- Acceptance criterion "G5 recalibrated to the real S1-S3 gap (measured 1.538 mHa) with
  provenance": done — margin lowered 5 mHa → 1 mHa (~35% headroom below the measured gap, ~1000×
  above the ~1 µHa DMRG/FCI cross-check noise floor established in `SPEC_nbn_dmrg_reference.md`).

## Gate run (fresh, own process, this session)

```
$ uv run python -m pytest tests/test_nbn_low_spin_spec.py -v
tests/test_nbn_low_spin_spec.py::test_G4_low_spin_cheap_D_is_pre_registered_inconclusive PASSED [ 50%]
tests/test_nbn_low_spin_spec.py::test_G5_committed_sector_is_not_the_cas_ground PASSED [100%]
======================== 2 passed in 484.09s (0:08:04) =========================
real    8m5.050s  |  user 14m38.382s  |  sys 0m4.959s
```

Run from a state equivalent to a fresh clone: `data/nbn_scf.chk` was present only via
`ensure_reference_data()`'s restore from the vendored `specs/nbn_scf_reference.chk` (same
mechanism chem-bbi verified for `test_nbn_dmrg_reference_spec.py`). Prior to this fix, the same
two tests failed exactly as the bead predicted (log captured earlier this session, 615.12 s,
`AssertionError`s on both G4's CONFIRM branch and G5's 5 mHa margin — matching the bead
description's numbers to full precision).

## What was decided not to do, and why

- **Table 2b (own-natural-orbitals cross-check) was not redone on the real geometry.** It's a
  second, more expensive cross-check (3 DMRG runs in a rotated basis) that isn't needed to satisfy
  chem-g1i's acceptance criteria (which ask for Tables 1-3 correction + G5 recalibration, not a
  from-scratch re-derivation of every sub-analysis). Flagged as an explicit gap in §3/§7 of the
  spec rather than silently dropped or backfilled with stale numbers.
- **Table 3's full nine-distance scan was not re-run against the real chk.** Only one new point
  (the real molecule's own bond length, 3.7255 Å) was added, specifically to test whether the
  scan pipeline reproduces the real chk at the one distance where that comparison is meaningful.
  It doesn't (chem-owr) — so re-running the other eight points would not have produced numbers
  comparable to Table 1 either; doing so was judged not to add information for this bead.
- **Did not diagnose the chem-owr discrepancy.** Ruled out the known MIC/unwrapping bug (verified
  `write_cif`'s own no-MIC/MIC distances agree to machine precision at 3.7255 Å), but did not go
  further (basis/ECP assignment, SCF spin-scan differences, or the real CIF's differing lattice
  constants (a=2.97204749, c=2.89967487) vs `write_cif()`'s defaults are all still open
  candidates) — filed as chem-owr (P3) rather than pursued here, since it's orthogonal to
  correcting the already-established real-chk-based Tables 1-2 this bead needed.
- **Did not re-run `tests/test_nbn_dmrg_reference_spec.py`** (chem-bbi's gate file). No shared
  library code was touched this session (`nbn_dmrg_reference.py`, `nbn_low_spin.py` unmodified);
  changes are confined to `SPEC_nbn_low_spin.md` and `tests/test_nbn_low_spin_spec.py`, which that
  other file does not import. chem-bbi already confirmed it 4/4 green.

## What could not be verified

- **Literal fresh `git clone` end-to-end.** Same caveat chem-bbi's handoff already recorded
  (R3 in the spec's caveats): `ensure_reference_data()`'s restore-from-vendored-copy behavior was
  verified functionally (this session's runs, and chem-bbi's), not via an actual `git clone` into
  a new directory.

## Git status — NOT committed, NOT pushed (per standing instructions)

```
$ git status
 M specs/SPEC_nbn_low_spin.md
 M tests/test_nbn_low_spin_spec.py
 M results/nbn_low_spin/runs.jsonl        (7 new provenance records appended)
 M .beads/interactions.jsonl              (bd bookkeeping: claim/notes on chem-g1i, chem-owr creation)
 M .beads/issues.jsonl                    (same)
```
(Plus pre-existing unrelated dirty files from other sessions — `hybrid_quantum_solver/dmrg_reference.py`,
`nbn_dmrg_reference.py`, `specs/SPEC_nbn_dmrg_reference.md`, `tests/test_nbn_dmrg_reference_spec.py`,
`tests/test_hubbard_lieb_wu_spec.py`, and untracked `specs/SPEC_regime_converged_undershoot.md`,
`tests/test_regime_converged_undershoot_spec.py`, `sandbox-handoffs/chem-mjz.md`,
`sandbox-handoffs/chem-bbi.md`, `specs/nbn_mp-2634.cif`, `specs/nbn_scf_reference.chk`,
`.claude/settings.local.json` — none of these were touched this session; see chem-bbi's own
handoff for that set.)

**Suggested commands for a human to run** (chem-g1i-scoped files only):

```bash
git add specs/SPEC_nbn_low_spin.md tests/test_nbn_low_spin_spec.py results/nbn_low_spin/runs.jsonl
git commit -m "nbn: correct SPEC_nbn_low_spin.md tables/gates to the real geometry (chem-g1i)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

The `.beads/*.jsonl` bookkeeping files will pick up the chem-g1i claim/close and the chem-owr
creation; stage those alongside if the bd sync convention here expects it, or let a subsequent
`bd`/human session handle it. The other pre-existing dirty/untracked files (chem-bbi's and
chem-mjz's) are out of scope for this commit — see `sandbox-handoffs/chem-bbi.md` for their own
suggested commit.

## Bead status

Closing **chem-g1i** — all three acceptance criteria are met with reproduced evidence:
1. Tables 1-3 recomputed on the real, vendored geometry (Table 1, Table 2; Table 3 gets one new
   real-distance point plus an explicit caveat; the old wrong-geometry tables are relabeled and
   preserved, not silently overwritten).
2. G5's threshold recalibrated to the real 1.538 mHa gap, with provenance (docstring cites the
   exact measured numbers and the noise floor used to justify the new margin).
3. `tests/test_nbn_low_spin_spec.py` full suite green (2/2) from a state equivalent to a fresh
   clone.

One out-of-scope finding (chem-owr, scan-vs-chk discrepancy) filed separately, not silently
folded in or left undocumented.
