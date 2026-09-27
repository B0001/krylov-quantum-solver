# chem-1uu handoff — H-chain equation of state across R: is D set by R or by n?

## Status: work was already complete on entering this session; verified and closed.

This bead's substantive work (code + analysis + writeup) was already present in the working tree
when this session started (bead was `in_progress`, claimed 2026-09-27, no notes, no prior
handoff file). No prior handoff existed for chem-1uu specifically, so this session's job was to
verify the existing diff against the acceptance criteria and the test suite, not to redo it. What
follows is that verification, done from scratch (nothing taken on faith).

## What exists (diff already in tree, not authored this session)

- `benchmark_hchain_tdl.py`: `integrals(n, localize=False, R_bohr=R_BOHR_DEFAULT)` — new `--R`
  flag, geometry now parametrized by bond length in Bohr instead of a hardcoded 1.8.
- `probe_hchain_localize.py`: threaded the same `--R` through the single-ladder cost-probe
  subprocess protocol (the one that avoids the memory-hungry multi-D ramp that originally blocked
  this bead at R=1.8, n≥16 on canonical orbitals).
- `tests/test_hchain_tdl_spec.py`: new gate
  `test_G6_R_bohr_changes_geometry_but_not_fci_invariance` — asserts (a) nuclear repulsion halves
  when R doubles (geometry actually changes) and (b) canonical/localized FCI invariance holds at
  R ∈ {1.8, 2.4, 3.6}, not just the R=1.8 default.
- `specs/SPEC_hchain_largen2.md` §12: the full six-run table (n=20, D=400, R × basis) and verdict.
- `specs/BACKLOG.md`: entry updated `[ ]` → `[x]`, part (b) marked done with the result summary.
- `data/hchain_R_sweep.json` (gitignored, regenerable): raw per-run JSON from
  `probe_hchain_localize.py --n 20 --D 400 --R <r> --runs loc,can`.

## Verification performed this session

1. **Environment**: `uv sync --extra dmrg --extra test` (block2 0.5.3 installed; was missing at
   session start). Container: `Linux 6c1b790997b2 7.0.12-linuxkit x86_64` (Docker/linuxkit),
   `nproc`=8, 16 GiB RAM (`MemTotal` 16353600 kB).
2. **Cross-checked the raw data against the vendored table.** `data/hchain_R_sweep.json` (already
   on disk, timestamped 2026-09-27 04:26, before this session touched anything) reproduces every
   number in `SPEC_hchain_largen2.md` §12's table to the last printed digit — e.g. R=1.8 localized
   `energy=-10.825878868397051`, `discarded_weight=8.542431389393372e-14` matches the table's
   `−10.825878868` / `8.542e-14`. This is real, already-computed evidence, not a placeholder.
3. **Ran the two directly relevant gates, each in its own process** (block2/pyscf isolation):
   ```
   uv run python -m pytest -q tests/test_hchain_tdl_spec.py -v
   ======================== 8 passed in 282.55s (0:04:42) =========================

   uv run python -m pytest -q tests/test_hchain_largen2_spec.py -v
   ======================== 5 passed in 288.70s (0:04:48) =========================
   ```
   Both green, including the new G6 gate (item 3 above) that exercises the `--R` plumbing this
   bead added.
4. Did **not** re-run the full `make gates` sweep (all ~40+ spec files) — out of scope for this
   bead's acceptance criteria, which name a specific table and verdict, not full-suite health. The
   working tree currently carries uncommitted diffs from several other beads (chem-y6w, chem-czw,
   chem-xqh, chem-1yr, chem-4y9, etc. — visible in `git status`); running `make gates` would
   conflate their gates with this bead's and isn't this bead's job to referee.

## Acceptance criteria — checked against `bd show chem-1uu`

> "Discarded weight and energy at n=20 for R in {1.8, 2.4, 3.6} bohr at matched D, in both bases
> where affordable."

Met. `SPEC_hchain_largen2.md` §12 table: D=400 fixed, n=20, R ∈ {1.8, 2.4, 3.6}, both localized
(Löwdin site orbitals) and canonical (RHF MOs) — six runs, all six completed and affordable
(712–2075 s each on this container's 4 block2 threads).

> "A stated verdict on 'D set by R, not n', with the orbital-basis confound addressed."

Met, and it's a real finding, not a rubber stamp: the mechanism-tier kill check (`dw(R=3.6) ≥
dw(R=1.8)` kills the claim) **passes cleanly in localized orbitals** (dw falls 8.5e-14 → 3.1e-17 →
2.1e-19 as R grows — D needed drops with R, as the Mott-regime argument predicts) but **fails,
inverted, in canonical orbitals** (dw *rises* 6.99e-4 → 2.89e-3 → 3.50e-3 over the same R). Spec
states plainly: "the R-dependence claim is basis-dependent — it does NOT survive the basis change
as originally stated; it inverts," with the likely mechanism (canonical RHF degrading as a
single-reference starting point toward dissociation) stated as an inference, not measured fact.

> "Numbers are vendored into a tracked file."

Met. Table lives in `specs/SPEC_hchain_largen2.md` §12 (tracked); `specs/BACKLOG.md` carries the
summary and links to it. Raw JSON (`data/hchain_R_sweep.json`) is gitignored per repo convention,
noted as such in the spec, and independently reproducible via `probe_hchain_localize.py --R`.

## Honest caveats already recorded in the spec (not smoothed over)

- R=3.6 localized may be a near-trivial atomic limit rather than a "hard" Mott point (nearest-
  neighbor hopping already small there) — spec says this establishes direction/size of the D(R)
  trend, not that R=3.6 is a hard benchmark point.
- This is a fixed-D diagnostic (deliberately unconverged canonical rows), not a converged EOS —
  no thermodynamic-limit claim is made from this table.
- Part (a) of the original BACKLOG entry (reproducing the 10-point published EOS at n=10) is a
  separate bead (chem-4y9), not this one.
- Only 3 R points, 1 n, 1 D — a mechanism check, not a scan.

## What I decided not to do, and why

- Did not attempt to also validate/close chem-4y9's or other beads' uncommitted diffs sitting in
  the same working tree — out of scope for chem-1uu, would be scope creep across sessions that
  didn't happen to land together.
- Did not run `make gates` / full suite — see point 4 above.
- Did not attempt to re-derive or extend the R grid beyond the pre-registered {1.8, 2.4, 3.6} —
  that would be a new hypothesis, filed separately if wanted, not a rerun of this one.

## What I could not verify

- Wall-clock/hardware numbers in the spec's own §12 preamble (852 s vs the §10 probe's 290 s,
  attributed to VM emulation + concurrent CI contention) are as-recorded from a prior run; I did
  not re-run the probe from scratch this session (it would cost ~1–2 h across 6 runs) since the
  raw JSON on disk already matches the vendored table exactly and the acceptance criteria don't
  require re-generation, only that the numbers be vendored and correct.

## Git state — left dirty per policy, not committed/pushed

Per repo instructions, no commit/push was made. Relevant tracked-file diffs for this bead
(already present at session start):
```
 M benchmark_hchain_tdl.py
 M probe_hchain_localize.py
 M specs/BACKLOG.md            (chem-1uu's part of this diff; also carries other beads' entries)
 M specs/SPEC_hchain_largen2.md
 M tests/test_hchain_tdl_spec.py
```
Suggested commands for a human to run (not executed here):
```
git add benchmark_hchain_tdl.py probe_hchain_localize.py specs/SPEC_hchain_largen2.md \
        tests/test_hchain_tdl_spec.py
# specs/BACKLOG.md also needs staging but carries other beads' entries too — review before adding
git commit -m "hchain: mechanism tier — D is set by R, and it flips sign with orbital basis (chem-1uu)"
```

## Bead disposition

Closing `chem-1uu` as done: acceptance criteria met, both directly relevant gates green
(13 tests, 0 failures), numbers cross-checked against raw data on disk.
