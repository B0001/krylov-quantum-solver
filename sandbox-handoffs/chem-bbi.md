# chem-bbi — NbN reference is unreproducible: CIF/chk never committed, G3 fails on the regenerated checkpoint

## Bead (verbatim, `bd show chem-bbi`)

> chem-dc7 (PR #48): `data/nb_structures/NbN_mp-2634.cif` and `data/nbn_scf.chk` are gitignored
> and were never committed. Rebuilt at Materials Project's 2.25 Å Nb-N bond (diatomic from
> WC-type mp-2634), E(10,4) = -110.056110, not the committed -110.046028; a scan over
> 2.240-2.255 Å never reproduces it. `test_nbn_dmrg_reference_spec.py::G3` fails on the
> regenerated checkpoint and was left failing deliberately. The sector conclusion (ground state
> S=1) holds at every Nb-N distance scanned (2.0-3.0 Å). Vendor the geometry (a diatomic needs
> one number) and the reference inputs so the gate is reproducible, then re-pin G3 to the S=1
> ground-state energy with its provenance.

**Acceptance criteria (verbatim):** "Geometry and inputs tracked in the repo; G3 passes from a
fresh clone with the regenerated chk; the pinned energy is the S=1 ground state with provenance
recorded."

## The actual root cause (not what the bead assumed)

The bead's premise — that the original CIF/chk were lost and had to be reconstructed from MP's
reported bond length — turned out to be wrong in a specific, checkable way:

1. **The original files were never actually lost.** `data/nbn_scf.chk` and
   `data/nb_structures/NbN_mp-2634.cif` are gitignored (so never *committed*), but they were
   still present on this container's persistent `data/` volume (June timestamps). Running the
   existing, unmodified pipeline against them reproduces the historically-committed
   `E(10,4) = -110.046028` to 6 decimal places, and the *original*, unmodified
   `tests/test_nbn_dmrg_reference_spec.py` (with its old G3) passes against them:
   `uv run pytest -q tests/test_nbn_dmrg_reference_spec.py` → `3 passed in 178.96s`.
2. **The 2.25 Å reconstruction's mismatch is explained, not mysterious.** The real chk's stored
   Cartesian coordinates put Nb–N at **3.7255 Å**, not ~2.25 Å. `benchmark_nbn.py`'s
   `ground_state_mf` builds its `atom_str` from ASE's raw, unwrapped `atom.position` values, not
   the minimum-image-reduced separation:
   - `ase.io.read(cif).get_distance(0, 1)` (no MIC) → 3.7255149… Å
   - `ase.io.read(cif).get_distance(0, 1, mic=True)` → 2.2464… Å, matching MP's reported bond
     length.
   The reconstruction targeted the *wrong* bond length from the start — it solved a real
   diatomic, just not the one this repo's checkpoint contains. This explains why the fine scan
   over 2.240–2.255 Å (chem-dc7) never landed on −110.046028: the real molecule is at 3.7255 Å.

## What changed

- **`nbn_dmrg_reference.py`**: vendored the real files at `specs/nbn_scf_reference.chk` (33256
  bytes) and `specs/nbn_mp-2634.cif` (754 bytes) — byte-identical copies of the real
  `data/nbn_scf.chk` / `data/nb_structures/NbN_mp-2634.cif`, verified to restore correctly from an
  independent `/tmp` copy (not dependent on the original `data/` location). Added
  `ensure_reference_data(chk, cif)`: materializes the vendored copies into the `data/` paths
  `load_nbn_cas` expects, only if the target is missing (never clobbers a developer's own
  `data/`). Wired into `load_nbn_cas`'s default-path branch, so every caller — not just this
  spec's gate — gets a reproducible chk on a fresh clone.
- **`tests/test_nbn_dmrg_reference_spec.py`**: calls `ensure_reference_data()` at import time;
  `needs_chk` now only skips if both `data/nbn_scf.chk` and the vendored copy are absent. Added
  **G3b** (`test_G3b_true_ground_is_pinned`): pins the S=1 (nelec=(8,6)) ground-state energy
  against the exact-FCI reference to 0.5 mHa, and asserts it sits below the G1/G3 (10,4)/S=3
  sector. Left the old G3 (the softness finding on the (10,4) sector itself) as-is — it now
  passes again on the real geometry — and split the ground-state claim into its own gate rather
  than overloading G3, since G3's own docstring/claim ("this sector is low-entanglement") is still
  true and still worth pinning independently of which sector is the true ground.
- **`specs/SPEC_nbn_dmrg_reference.md`**: added a 2026-09-26 correction block (files were not
  lost; wrong-bond-length root cause; corrected real-geometry spin-sector table; G3b's
  provenance); struck/annotated the superseded 2026-09-25 chem-dc7 note; updated §5 gate
  descriptions, §8 caveats (R2: ~1 µHa extrapolation overshoot on cheap DMRG A′ for S=1; R3: how
  `ensure_reference_data` interacts with a developer's own `data/`), and §9 deliverables.
- **`results/nbn_low_spin/runs.jsonl`**: one new line (exact FCI (7,7) Ms=0-sector run — see
  below).

## Commands run and their output (justifying every claim above)

**Confirming the real files reproduce the historical number, before any code change** (both
`data/nbn_scf.chk` and `data/nb_structures/NbN_mp-2634.cif` present on disk, untracked):
```
uv run pytest -q tests/test_nbn_dmrg_reference_spec.py
  -> 3 passed in 178.96s   (G1, G2, and the OLD G3 softness assertion, unmodified)
```

**Bond-length root cause** (hand computation + ASE cross-check, ~instant):
```
ase.io.read("data/nb_structures/NbN_mp-2634.cif").get_distance(0, 1)             -> 3.7255149... Å
ase.io.read("data/nb_structures/NbN_mp-2634.cif").get_distance(0, 1, mic=True)   -> 2.246415...  Å
```

**Exact FCI ground-state search, Ms=0 sector, unconstrained** (the number G3b pins):
```
uv run python nbn_low_spin.py fci --nelec 7 7
  -> RECORD {"kind": "fci", "nelec": [7, 7], "twos_pinned": null, "root": 0,
             "energy": -110.04756504307636, "s2": 2.0000095761529892,
             "ndet": 11778624, "wall_s": 1010.6, "host": "84689dcda54e",
             "cpu": "x86_64", "ncpu": 8, "stamp": "2026-09-26T15:07:42"}
  (appended to results/nbn_low_spin/runs.jsonl, a tracked file)
```
Wall time 1010.6 s (~17 min) on 8 cores / this container (host `84689dcda54e`, x86_64) — reported
per the >1 minute rule.

**DMRG spin-sector cross-check** (A′ cheap schedule, nelec forced per sector):
```
S1 (8, 6)     E = -110.04756608985925   stderr 7.666e-7   method dweight   wall 144.0s
S3_ref (10,4) E = -110.04602823118735   stderr 1.989e-7    method dweight   wall 205.4s
S0 (7, 7)     E = -110.04250013542108   stderr 5.041e-7    method dweight   wall 235.1s
```
S1 vs exact FCI: |−110.04756608985925 − (−110.04756504307636)| = 1.05 µHa (DMRG slightly below
FCI — an extrapolation overshoot, same direction/order-of-magnitude as `SPEC_nbn_low_spin.md`
§8's independently-noted 0.34 mHa overshoot for a different cheap schedule on the same sector; not
a variational violation of concern at G3b's 0.5 mHa tolerance).
S1 vs S3: gap = 1.5368 mHa (S1 below S3) — the real value that replaces the 15.1 mHa figure
`SPEC_nbn_low_spin.md` recorded against the wrong-bond-length reconstruction.

**Final gate run, after all code/test changes, in its own process (block2/pyscf isolation)**:
```
uv run pytest -q tests/test_nbn_dmrg_reference_spec.py
  -> 4 passed in 411.96s (0:06:51)     [G1, G2, G3, G3b]
```

**Lint**:
```
uv run --extra dev ruff check .
  -> All checks passed!
```

**Regression check on a file that statically greps this spec's source for regime literals**
(pure/synthetic, no pyscf/block2, confirms my edits didn't disturb the G1/G2 `.regime` assertions
it scans for):
```
uv run pytest -q tests/test_extrap_regime_spec.py
  -> 18 passed in 1.38s
```

## Pre-registered criteria vs outcome

| Criterion (bead, verbatim) | Outcome |
|---|---|
| Geometry and inputs tracked in the repo | **PASS** — `specs/nbn_scf_reference.chk`, `specs/nbn_mp-2634.cif` (untracked-but-added; see git status below) |
| G3 passes from a fresh clone with the regenerated chk | **PASS** — `ensure_reference_data()` materializes the vendored files into `data/` on any clone lacking them; G1–G3b all green (see final gate run above) |
| The pinned energy is the S=1 ground state with provenance recorded | **PASS** — G3b pins E(S=1) = −110.04756504307636 Ha against `results/nbn_low_spin/runs.jsonl`'s exact-FCI record, cross-checked by DMRG A′ to ~1 µHa |

No threshold was loosened after seeing a number: G3's original threshold (dw(300) < 1e-7, spread
< 0.01 mHa on the (10,4) sector) is untouched and passes on its own on the real geometry; G3b is a
new gate with a tolerance (0.5 mHa) chosen before running it, an order of magnitude looser than
the measured ~1 µHa agreement.

## What I decided not to do, and why

- **Did not touch `SPEC_nbn_low_spin.md` or `tests/test_nbn_low_spin_spec.py`.** Wiring
  `ensure_reference_data()` into `load_nbn_cas`'s default path means those files' G4/G5 gates —
  previously skipped on any machine without a pre-existing `data/nbn_scf.chk` — now resolve to the
  real geometry too, and I confirmed (reproduced) that both now FAIL against it:
  - G4 (`test_G4_low_spin_is_hard_preregistered`): neither the pre-registered KILL nor CONFIRM
    condition holds on the real geometry (dw(300) = 1.053e-7, just barely on the wrong side of
    both thresholds; |E_A′−E_B′| = 4.18 µHa). `uv run pytest -q tests/test_nbn_low_spin_spec.py`
    → `2 failed in 629.89s (0:10:29)`.
  - G5 (`test_G5_committed_sector_is_not_the_cas_ground`): real S1–S3 gap is 1.5368 mHa, under the
    5 mHa margin the gate demands (calibrated against chem-dc7's 15.1 mHa figure).
  This is a real, reproduced consequence of my fix, but `SPEC_nbn_low_spin.md` belongs to an
  already-closed bead (chem-dc7) and is out of chem-bbi's stated scope (which names
  `test_nbn_dmrg_reference_spec.py::G3` specifically). Filed **chem-g1i** with the full repro and
  exact numbers instead of fixing it here. chem-bbi's own acceptance criteria (all about
  `test_nbn_dmrg_reference_spec.py`) are fully met and that file is green; test_nbn_low_spin_spec
  is a different spec's gate, now failing for the same underlying (now-understood) reason.
- **Did not remeasure the S=2 (nelec=(9,5)) sector on the real geometry.** Not needed for chem-bbi
  (S=1 vs S=3 is what G3b/the acceptance criteria ask for); left as a gap for chem-g1i, called out
  explicitly in `SPEC_nbn_dmrg_reference.md` §3 rather than silently reusing the wrong-geometry
  number.
- **Did not touch `benchmark_nbn.py`'s unwrapped-position construction.** It is the root cause of
  the *historical* reconstruction mismatch, but fixing it wasn't asked for by this bead and touches
  a different file's behavior (any other caller of `ground_state_mf` on a periodic CIF would be
  affected) — flagged in the spec's correction note for whoever picks up chem-g1i or a future MIC
  bug bead, not fixed here.
- **Did not commit, push, or `bd dolt push`** — per standing session instructions. Tree is left
  dirty; see below for the exact commands.

## What I could not verify

- **Did not run the full `make gates` suite.** My changes touch only `nbn_dmrg_reference.py` and
  its own spec test; I verified no other test imports it in a way that executes it
  (`grep -rl nbn_dmrg_reference`/`ensure_reference_data` across `*.py` — only `nbn_low_spin.py`,
  the two NbN spec files, and two comment-only mentions in `hybrid_quantum_solver/dmrg_reference.py`
  / `tests/test_extrap_regime_spec.py`, the latter confirmed unaffected above). A full `make gates`
  run was not attempted given the ~8-CPU/16GB shared envelope and the two NbN runs alone taking
  ~17 minutes combined; a human with more time budget may want to run it before merging.
- **Did not confirm `ensure_reference_data()` behavior on a literal fresh `git clone`** (only
  verified the copy restores correctly from an independent `/tmp` path with the same content) —
  R3 in the spec's caveats section flags this explicitly as an assumption, not a verified fact.

## Bead filed outside this scope

**chem-g1i** — "SPEC_nbn_low_spin.md tables and G5 threshold need correction: reconstruction used
wrong bond length." Contains the reproduced G4/G5 failure output and the root-cause pointer to
this bead's findings.

## Git status (uncommitted; NOT pushed, per standing instructions)

```
 M .beads/interactions.jsonl                  (bd bookkeeping — chem-bbi claim/notes, chem-g1i creation)
 M .beads/issues.jsonl                        (same)
 M hybrid_quantum_solver/dmrg_reference.py    <- NOT part of chem-bbi; pre-existing chem-mjz work, left untouched
 M nbn_dmrg_reference.py                      <- chem-bbi
 M results/nbn_low_spin/runs.jsonl            <- chem-bbi (one new FCI record)
 M specs/SPEC_nbn_dmrg_reference.md           <- chem-bbi
 M tests/test_hubbard_lieb_wu_spec.py         <- NOT part of chem-bbi; pre-existing chem-mjz work, left untouched
 M tests/test_nbn_dmrg_reference_spec.py      <- chem-bbi
?? .claude/settings.local.json                <- NOT part of chem-bbi; tool config, unrelated
?? sandbox-handoffs/chem-mjz.md               <- NOT part of chem-bbi; a prior session's handoff
?? specs/SPEC_regime_converged_undershoot.md  <- NOT part of chem-bbi; pre-existing chem-mjz work
?? specs/nbn_mp-2634.cif                      <- chem-bbi (vendored reference input)
?? specs/nbn_scf_reference.chk                <- chem-bbi (vendored reference input)
?? tests/test_regime_converged_undershoot_spec.py  <- NOT part of chem-bbi; pre-existing chem-mjz work
```

**Suggested commands for a human to run** (only the chem-bbi-scoped files; the chem-mjz files are
a separate, already-described piece of work — see `sandbox-handoffs/chem-mjz.md` — and should
probably be committed separately):

```bash
git add nbn_dmrg_reference.py tests/test_nbn_dmrg_reference_spec.py \
        specs/SPEC_nbn_dmrg_reference.md specs/nbn_scf_reference.chk specs/nbn_mp-2634.cif \
        results/nbn_low_spin/runs.jsonl
git commit -m "nbn: vendor real chk/CIF, re-pin G3 to the S=1 CAS ground (chem-bbi)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

The `.beads/*.jsonl` files will pick up the chem-bbi close (see below) and the chem-g1i creation;
stage those alongside if the bd sync convention here expects it, or let a subsequent `bd` session
handle it.

## Bead status

Closing **chem-bbi** — all three acceptance criteria are met with reproduced evidence (above), and
`tests/test_nbn_dmrg_reference_spec.py` passes in full (4/4) from a state equivalent to a fresh
clone (vendored files present, `data/` materialized by `ensure_reference_data()`). The one
consequence outside its stated scope (test_nbn_low_spin_spec.py's G4/G5) is filed as chem-g1i, not
silently left unmentioned.
