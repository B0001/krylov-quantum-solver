# chem-czw handoff — NbN CAS(14,14): re-derive on the other UHF minimum, decide whether to switch reference

**Bead:** chem-czw (P3, task) — "Re-derive NbN CAS(14,14) reference numbers on the deterministic
lower-energy UHF minimum" (follow-up to chem-y6w).

**Resolution:** did the re-derivation for all four spin sectors on the deterministic lower-energy
minimum's orbitals, cross-checked against the existing (vendored-chk) numbers, and answered the
bead's question: **keep the vendored chk as the primary reference.** Not a cost-based deferral this
time (chem-y6w's reason) — a substantive one: the vendored chk's orbitals give a strictly lower
(more favorable) CASCI(14,14) energy than the alternative minimum's orbitals, in every sector
checked.

---

## 1. What changed

### 1a. Regenerated the lower-energy minimum's chk, confirmed it's the same minimum chem-owr/chem-y6w characterized

New one-off script `scripts/gen_nbn_lower_minimum_chk.py`: monkey-patches `benchmark_nbn.CHKFILE`,
calls `benchmark_nbn.ground_state_mf("specs/nbn_mp-2634.cif")` (the chem-y6w-fixed,
`lib.num_threads(1)`-pinned deterministic path), asserts `mf.e_tot` matches chem-owr/chem-y6w's
characterized lower minimum (`-110.02813374576796 Ha`) within `1e-9` Ha (loose enough to absorb
~5.7e-14 Ha cross-container float noise observed on this run, five orders of magnitude tighter
than the 0.548 mHa gap to the *other* minimum — still cleanly distinguishes "same minimum, float
noise" from "landed on a third solution"). The chk itself was **not** kept as a committed artifact
(see §5) — `results/` otherwise holds only text/CSV/log files, and the chk is deterministically
regenerable from tracked inputs (the script + `specs/nbn_mp-2634.cif`), unlike the vendored chk
which is frozen precisely *because* a fresh regen does not reproduce it.

### 1b. Added `--chk` to `nbn_low_spin.py`'s `fci` subcommand (previously `dmrg`-only)

```
$ git diff -- nbn_low_spin.py
```
Two small hunks: `cmd_fci` now passes `chk=args.chk` through to `load_nbn_cas` and records it in
the output JSON row; the `fci` argparse subparser gained `--chk` (default `"data/nbn_scf.chk"`,
backward-compatible) — mirrors the pattern already on `cmd_dmrg` exactly.

### 1c. Ran the full four-sector battery on the lower minimum's orbitals

| S | nelec | method | E (Ha) | wall | vs Table 1 (vendored chk, same sector) |
|---|---|---|---|---|---|
| 1 | (8,6) | exact FCI | −110.04577676472408 | 259.7s / 8 cores | vendored **1.788 mHa lower** |
| 2 | (9,5) | exact FCI | −110.04546990757419 | 192.2s / 8 cores | vendored **1.516 mHa lower** |
| 3 | (10,4) | exact FCI | −110.04119080464493 | 5.6s / 8 cores | vendored **4.838 mHa lower** (matches chem-owr's originally-cited "4.84 mHa" almost exactly) |
| 0 | (7,7) | DMRG A, perD 400/800/1200 | −110.0410570692685 | 901.6s / 4 threads, dw(1200)=8.2e-10 | vendored **1.440 mHa lower** |
| 0 | (7,7) | DMRG B, ramp 300/600/1200 | −110.04105709905623 | 444.0s / 4 threads, dw(1200)=7.7e-10 | agrees with A to 2.98e-8 Ha |

Commands (each `--chk results/nbn_lower_minimum/nbn_scf_lower.chk` while it existed; regenerate
via `scripts/gen_nbn_lower_minimum_chk.py` first — see §5):

```bash
uv run python nbn_low_spin.py fci --nelec 8 6 --chk <regen'd chk path>
uv run python nbn_low_spin.py fci --nelec 9 5 --chk <regen'd chk path>
uv run python nbn_low_spin.py fci --nelec 10 4 --chk <regen'd chk path>
uv run python nbn_low_spin.py dmrg --schedule A --nelec 7 7 --threads 4 --chk <regen'd chk path>
uv run python nbn_low_spin.py dmrg --schedule B --nelec 7 7 --threads 4 --chk <regen'd chk path>
```

Raw JSON records (with `"chk": "results/nbn_lower_minimum/nbn_scf_lower.chk"`) are in
`results/nbn_low_spin/runs.jsonl` (tracked, append-only — this session's 5 new lines are a pure
append, verified via `git diff --stat`, no existing lines touched).

### 1d. Findings

1. **Ordering is orbital-choice-independent.** S=1 < S=2 < S=3 < S=0 on the lower minimum's
   orbitals too, matching Table 1. The (7,7) sector's softness (two independent DMRG schedules
   agreeing to 2.98e-8 Ha) also survives. Nothing here undermines chem-dc7/chem-g1i's existing
   sector-ordering or hardness conclusions.
2. **The vendored chk's orbitals give a strictly lower (more favorable) CASCI(14,14) energy than
   the lower-SCF-energy minimum's orbitals, in every sector checked — by 1.44 to 4.84 mHa** —
   despite the vendored chk having a *higher* raw whole-molecule SCF energy (0.548 mHa). This is
   the opposite direction from the SCF-level comparison. Caveat (pre-registered as informal, not
   resolved): comparing CASCI totals across two different UHF references at fixed active-space
   size is a standard, useful empirical quality criterion but **not** a rigorous variational bound
   over orbital choice — a rigorous comparison needs CASSCF orbital relaxation, out of scope here
   (as it already was in `SPEC_nbn_dmrg_reference.md` §7).

**Decision: keep `specs/nbn_scf_reference.chk` / `data/nbn_scf.chk` as the primary reference.**
This upgrades chem-y6w's deferral (compute cost alone) to a substantive reason: the committed chk
is not merely cheaper to keep, it demonstrably gives a better-converged CAS(14,14) description in
every sector tested, by a comfortable margin (1.4-4.8 mHa, vs. the ~1e-7 to 1e-9 Ha noise floor of
these calculations).

### 1e. Specs updated (decision recorded in specs, not just this handoff)

- `specs/SPEC_nbn_low_spin.md` — new **§0b** with the full method, Table 4 (the cross-check table
  above), findings, and decision; top status line updated; **§6** gained **G6** (the new gate,
  below); **§7** "Out of scope" bullet updated from "declined" to "DONE (chem-czw)".
- `specs/SPEC_nbn_dmrg_reference.md` — **§8 R4** appended with the chem-czw resolution paragraph
  (the other minimum gives strictly higher CASCI energies than G1-G3b's own numbers, in every
  sector).
- `specs/SPEC_nbn_scf_determinism.md` — **§7** out-of-scope bullet updated to note chem-czw did
  the re-derivation and confirmed the decision.
- `specs/BACKLOG.md` — new `[x]` Done entry (top of `## Done`, right after chem-y6w's).

### 1f. New regression gate: `tests/test_nbn_czw_lower_minimum_spec.py` (own process)

Pyscf + qiskit only (imports `benchmark_nbn` → `hybrid_quantum_solver.molecular_hamiltonian` →
qiskit, for the fixture-gated live-regeneration sub-test). **No block2/DMRG import anywhere in
this file** (`nbn_low_spin.py` and `hybrid_quantum_solver.dmrg_reference` both import block2
lazily inside function bodies, not at module scope) — still kept as its own file per the
`test_*_spec.py` one-file-per-process convention, so it never risks sharing a process with a
block2 import regardless.

Five pinned-literal regression checks (Table 4's numbers, cheap, no compute) plus one fixture-gated
live check:

- Ordering (S=1 < S=2 < S=3 < S=0) holds on the lower minimum's orbitals.
- The two independent (7,7) DMRG schedules agree (< 1e-4 Ha).
- The vendored chk's CASCI energy is strictly lower than the lower minimum's, by at least 1 mHa,
  in each of the three exact-FCI sectors (parametrized) — **this is the decision-bearing check.**
- Same, for the DMRG S=0 sector.
- `test_lower_minimum_chk_regenerates_and_reproduces_S3_FCI` (skipped without
  `specs/nbn_mp-2634.cif`/`specs/nbn_scf_reference.chk`): live end-to-end regen of the *other*
  minimum (~10-20s: one tight SCF + one cheap (10,4) FCI) reproducing both the pinned SCF energy
  and the pinned (10,4) FCI energy within `1e-6` Ha — catches a regression in
  `benchmark_nbn._tight_scf`'s determinism fix or in `load_nbn_cas`/`spin_fci` without repeating
  the full ~25-minute four-sector re-derivation every cycle.

---

## 2. Gate-run evidence (verbatim)

```
$ uv run python -m pytest tests/test_nbn_czw_lower_minimum_spec.py -v
============================= test session starts ==============================
collected 7 items

tests/test_nbn_czw_lower_minimum_spec.py::test_lower_minimum_ordering_matches_vendored_S1_lowest PASSED [ 14%]
tests/test_nbn_czw_lower_minimum_spec.py::test_lower_minimum_dmrg_schedules_agree PASSED [ 28%]
tests/test_nbn_czw_lower_minimum_spec.py::test_vendored_chk_gives_lower_casci_energy_than_the_other_minimum[nelec0] PASSED [ 42%]
tests/test_nbn_czw_lower_minimum_spec.py::test_vendored_chk_gives_lower_casci_energy_than_the_other_minimum[nelec1] PASSED [ 57%]
tests/test_nbn_czw_lower_minimum_spec.py::test_vendored_chk_gives_lower_casci_energy_than_the_other_minimum[nelec2] PASSED [ 71%]
tests/test_nbn_czw_lower_minimum_spec.py::test_vendored_chk_gives_lower_dmrg_energy_for_S0_too PASSED [ 85%]
tests/test_nbn_czw_lower_minimum_spec.py::test_lower_minimum_chk_regenerates_and_reproduces_S3_FCI PASSED [100%]

============================== 7 passed in 11.59s ==============================
```

Lint on touched/new files:

```
$ uv sync --extra dev   # (re-synced --extra dmrg --extra test --extra dev together afterward —
                         # a bare `--extra dev` sync deselects dmrg/test and uninstalls block2/
                         # pytest; noted here so a future session doesn't repeat that surprise)
$ uv run ruff check tests/test_nbn_czw_lower_minimum_spec.py nbn_low_spin.py scripts/gen_nbn_lower_minimum_chk.py
All checks passed!
```

Existing gates re-run to confirm no regression (each its own process, per the block2/pyscf/
qiskit-aer isolation rule):

```
$ uv run python -m pytest tests/test_nbn_low_spin_spec.py -q
..                                                                       [100%]
2 passed in 505.32s (0:08:25)
```

Confirms G4/G5 (Table 1's own gates, still anchored to the vendored chk) are unaffected by this
bead's additive `--chk` edit to `nbn_low_spin.py::cmd_fci` — expected, since the default value is
unchanged and G4/G5 never pass `--chk` explicitly.

---

## 3. Every number, with regenerating command and vendor location

| Number | Value | Command | Vendor location |
|---|---|---|---|
| Lower minimum SCF energy (reproduced) | `-110.0281337457679 Ha` (diff 5.68e-14 Ha from chem-owr/chem-y6w's `-110.02813374576796`) | `uv run python scripts/gen_nbn_lower_minimum_chk.py` | not vendored — regenerated fresh, script is tracked |
| Lower minimum, S=1 (8,6) exact FCI | `-110.04577676472408 Ha` | `nbn_low_spin.py fci --nelec 8 6 --chk <regen'd>` | `results/nbn_low_spin/runs.jsonl` |
| Lower minimum, S=2 (9,5) exact FCI | `-110.04546990757419 Ha` | `nbn_low_spin.py fci --nelec 9 5 --chk <regen'd>` | `results/nbn_low_spin/runs.jsonl` |
| Lower minimum, S=3 (10,4) exact FCI | `-110.04119080464493 Ha` | `nbn_low_spin.py fci --nelec 10 4 --chk <regen'd>` | `results/nbn_low_spin/runs.jsonl`, also pinned as a literal in the new gate |
| Lower minimum, S=0 (7,7) DMRG A | `-110.0410570692685 Ha` | `nbn_low_spin.py dmrg --schedule A --nelec 7 7 --threads 4 --chk <regen'd>` | `results/nbn_low_spin/runs.jsonl` |
| Lower minimum, S=0 (7,7) DMRG B | `-110.04105709905623 Ha` | `nbn_low_spin.py dmrg --schedule B --nelec 7 7 --threads 4 --chk <regen'd>` | `results/nbn_low_spin/runs.jsonl` |
| Vendored-chk numbers (Table 1, unchanged, cross-checked against) | −110.04756504307636 / −110.04698592997889 / −110.04602841723886 / −110.04249952166984 Ha (S=1/2/3/0) | not recomputed this session | `specs/SPEC_nbn_low_spin.md` Table 1 (chem-g1i) |
| Per-sector shift (vendored minus lower minimum) | +1.788 / +1.516 / +4.838 / +1.440 mHa (S=1/2/3/0; positive = vendored lower/better) | arithmetic on the rows above | `specs/SPEC_nbn_low_spin.md` §0b Table 4 |

---

## 4. Pre-registered criteria and pass/fail

The bead asked to (1) decide whether the lower-energy minimum should replace the vendored chk, and
(2) if so, redo the FCI/DMRG table for every spin sector, cross-checking against the existing
numbers rather than silently overwriting them.

1. **Decide whether to switch reference** — **PASS, decided: no.** Recorded in
   `specs/SPEC_nbn_low_spin.md` §0b (not just this handoff), with the CASCI-energy-comparison
   finding as the substantive reason.
2. **If switching, redo the table for every sector, cross-checked, not overwritten** — **N/A** (the
   decision was not to switch) **but done anyway**, since answering (1) honestly required actually
   running all four sectors on the other minimum first. Table 4 in `SPEC_nbn_low_spin.md` §0b is
   explicitly a cross-check table, side-by-side with Table 1; nothing in Table 1,
   `SPEC_nbn_dmrg_reference.md` G1-G3b, or any other existing number was overwritten.
3. **Test suite passes** — **PASS.** New gate `tests/test_nbn_czw_lower_minimum_spec.py`: 7/7
   (§2). Pre-existing gate `tests/test_nbn_low_spin_spec.py` (whose numbers this bead cross-checks
   against, and whose host file `nbn_low_spin.py` this bead edited): 2/2, confirming the additive
   `--chk` option (default unchanged) did not regress G4/G5 (§2). Ruff clean on all touched/new
   files (§2). Full `make gates`/`make test` NOT run (would re-run every `test_*_spec.py` file in
   the repo, including many unrelated to this bead, each its own process — costly; scoped to the
   NbN-relevant files plus lint, matching chem-y6w's own precedent).

---

## 5. What was decided not to do, and why

- **Did not keep `results/nbn_lower_minimum/nbn_scf_lower.chk` as a committed artifact.** Every
  other file directly under `results/` is text/CSV/log (`find results -maxdepth 2 -type f`); this
  chk would have been the only binary. Unlike the vendored `specs/nbn_scf_reference.chk` — frozen
  precisely because a fresh regeneration does *not* reproduce it — this chk **is** deterministically
  regenerable from tracked inputs (`scripts/gen_nbn_lower_minimum_chk.py` + `specs/
  nbn_mp-2634.cif`), so committing the binary buys nothing a human couldn't get by running the
  script; the new gate's fixture-gated sub-test already proves the regeneration reproduces the
  pinned numbers on demand. Deleted it from the working tree before finishing (`rm -rf
  results/nbn_lower_minimum`) rather than leave a stray, ambiguous binary for a human reviewing
  `git status`.
- **Did not re-derive DMRG for the lower minimum's (10,4)/(9,5)/(8,6) sectors** — redundant, exact
  FCI already covers them at a fraction of the cost (same approach Table 1 itself used).
- **Did not attempt CASSCF orbital relaxation on either minimum** — would settle the "which
  reference is actually better" question rigorously (a true variational comparison), but is out of
  scope per `SPEC_nbn_dmrg_reference.md` §7 (unchanged by this bead).
- **Did not re-derive `SPEC_nbn_dmrg_reference.md`'s own G1-G3b gate *file*** on the lower minimum —
  Table 4's (10,4) number is the value those gates would need if the reference were switched, but
  since the decision is to keep the vendored chk, the gate file itself is left pointed at the
  vendored chk, unchanged.
- **Did not run full `make gates`/`make test`** — see §4 point 3.

---

## 6. What could not be verified

- **Whether the CASCI-energy-comparison finding would survive CASSCF orbital relaxation** — it is
  an informal, standard cross-reference criterion at fixed active-space size, not a rigorous bound;
  a rigorous answer needs CASSCF on both minima, out of scope here (as in `SPEC_nbn_dmrg_reference.md`
  §7, unchanged).
- **Whether either minimum is the *global* UHF minimum** (vs. some third, still-lower solution) —
  not investigated, same caveat chem-y6w left open.
- Nothing outstanding here — `tests/test_nbn_low_spin_spec.py` completed green (2/2, §2/§4) before
  this handoff was finalized.

---

## 7. Git status / suggested commands (NOT executed — human review required)

Per this session's git policy, no commit/push was made. Files touched by this bead's work:

```
 M nbn_low_spin.py
 M results/nbn_low_spin/runs.jsonl
 M specs/BACKLOG.md
 M specs/SPEC_nbn_dmrg_reference.md
 M specs/SPEC_nbn_low_spin.md
?? specs/SPEC_nbn_scf_determinism.md          (pre-existing untracked file from chem-y6w; one line edited by chem-czw in §7)
?? scripts/gen_nbn_lower_minimum_chk.py
?? tests/test_nbn_czw_lower_minimum_spec.py
?? sandbox-handoffs/chem-czw.md
```

(Every other modified/untracked file in `git status` — `be2_cbs.py`, `cross_check.py`,
`lambda_ladder.py`, `results/be2_avas_cbs/`, `specs/SPEC_be2_cbs.md`, etc. — predates this session:
other beads' uncommitted work already sitting in this shared container. Not touched here.)

Suggested commands for a human to review and run:

```bash
git add nbn_low_spin.py results/nbn_low_spin/runs.jsonl specs/BACKLOG.md \
        specs/SPEC_nbn_dmrg_reference.md specs/SPEC_nbn_low_spin.md \
        specs/SPEC_nbn_scf_determinism.md scripts/gen_nbn_lower_minimum_chk.py \
        tests/test_nbn_czw_lower_minimum_spec.py sandbox-handoffs/chem-czw.md
git commit -m "chem-czw: re-derive NbN CAS(14,14) on the other UHF minimum, keep vendored chk as reference

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
bd close chem-czw --reason="Re-derived all four spin sectors (exact FCI S=1/2/3, two independent DMRG schedules S=0) on the deterministic lower-energy UHF minimum chem-y6w characterized. Decision: keep the vendored chk as primary reference -- its orbitals give a strictly lower (more favorable) CASCI(14,14) energy than the alternative minimum's in every sector checked (1.44-4.84 mHa), despite a higher raw SCF energy. Cross-check table in SPEC_nbn_low_spin.md Section 0b; regression gate tests/test_nbn_czw_lower_minimum_spec.py green (7/7)."
```

All gates referenced above are confirmed green as of this handoff (§2/§4) — these commands are
ready to run as-is; no outstanding confirmation is needed.

---

## ADDENDUM 2026-09-27 (local macOS re-verification)

1. chem-y6w's `_tight_scf` rewrite had `mf.kernel(dm0=dm)` dedented out of the stability loop
   (UnboundLocalError when the first check is stable; no re-convergence otherwise). Fixed.
2. With the fix, macOS's deterministic regeneration lands on the VENDORED minimum
   (-110.02758538272664), not the lower one. Which minimum you get depends on the BLAS build.
   The regen gate now accepts either characterized minimum and checks its S=3 FCI (macOS
   reproduces Table 1's -110.04602841723886).
3. The lower-minimum numbers in §1c were produced by the buggy loop on Linux; they have not been
   re-run with the fixed loop (filed as a follow-up bead).
