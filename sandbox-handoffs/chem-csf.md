# chem-csf handoff — Is NbN's flagged 'hard multireference' low-spin (7,7) sector actually hard?

**Bead:** chem-csf (P2, task) — "run the identical A′/B′ DMRG schedule pair already used
elsewhere in `SPEC_nbn_dmrg_reference.md`, forcing nelec=(7,7), and report whether the sector
is genuinely hard or was just never run."

**Resolution: the bead's question was already fully answered in the repo before this session
started, by a different (but overlapping) bead chain — chem-dc7 → chem-g1i → chem-czw,
documented in `SPEC_nbn_low_spin.md`. This session did not do new physics; it independently
re-verified that existing answer on fresh hardware and then fixed the stale documentation that
caused the duplicate filing in the first place.**

---

## 1. What I found before doing anything

`bd show chem-csf` and `specs/BACKLOG.md`'s open item (the one that generated this bead) both
quote `SPEC_nbn_dmrg_reference.md` §2/§7 as saying the (7,7) sector was "never run, never
gated." That line in `SPEC_nbn_dmrg_reference.md` was accurate when written (2026-07 era) but
stale by 2026-09-27: `specs/SPEC_nbn_low_spin.md` (bead chem-dc7, corrected by chem-g1i, cross-
checked by chem-czw) had already run **exactly the check this bead asks for**, plus more:

- All four `SCHEDULES` (`A`, `B`, `A'`, `B'`) run against `nelec=(7,7)` forced, each session's
  own process, reusing `load_nbn_cas`/`run_schedule` unmodified (no new machinery, same as this
  bead's own "Reuse" instruction).
- Cheap D≤300 schedules (`A'`/`B'`) report `dw(D=300)` and cross-schedule spread
  `|E_A' − E_B'|` — `SPEC_nbn_low_spin.md` §3 Table 2, gated by `tests/test_nbn_low_spin_spec.py::test_G4_*`.
- The verdict is recorded: **neither the pre-registered KILL (`dw<1e-7` and spread `<1e-5 Ha`)
  nor CONFIRM (`dw>1e-5` or `|E_A'-E_B'|>1e-4 Ha`) fires at D≤300** — a genuine gap in the
  original two-sided design, not a threshold nudge.
- The bead's own fallback question — "does the 3.5 mHa gap survive at converged D, checking for
  a chem-4e9/chem-mjz-style stalled-schedule mislabel" — is also already answered: headline
  D≤1200 schedules (`A`/`B`, chem-g1i, `results/nbn_low_spin/runs.jsonl`) converge to
  `dw ≈ 3e-10`, two independent schedules agreeing to 0.3 nHa, with `dw` monotonically
  *decreasing* with D at every step (not a stall). **Verdict: the sector is soft, like every
  other sector in this CAS — the "hard multireference benchmark" framing does not hold.**

So the deliverable this bead asks for already exists, in more depth than requested (it went to
D=1200 specifically to resolve the D≤300 ambiguity this bead's own acceptance criteria treats as
a possible outcome). What did **not** exist: (a) independent re-verification on hardware other
than the container that produced the original numbers, and (b) any fix to the stale text in
`SPEC_nbn_dmrg_reference.md` and the `BACKLOG.md` entry that caused this bead to be filed as if
the question were open. I did both.

---

## 2. What I did this session

### 2a. Re-verified, did not re-derive

Ran `tests/test_nbn_low_spin_spec.py` (own process, per the block2/pyscf isolation rule) fresh
on this container — not reusing any cached result:

```
$ uv sync --extra dmrg --extra test
 + block2==0.5.3  (+ 6 other packages)
$ rm -rf .dmrg_tmp/nbn_Ac_lowspin_77 .dmrg_tmp/nbn_Bc_lowspin_77 .dmrg_tmp/nbn_Ac_lowspin_86 .dmrg_tmp/nbn_Ac_lowspin_scf
$ time uv run pytest -q tests/test_nbn_low_spin_spec.py -v
tests/test_nbn_low_spin_spec.py ..                                       [100%]
2 passed in 492.55s (0:08:12)
real    8m13.616s
```

This runs `A'`/`B'` (D≤300) on `nelec=(7,7)` (G4) and `A'` on `nelec=(8,6)`/`(10,4)`/`(7,7)`
(G5), exercising both this bead's exact ask and the "committed sector is not the CAS ground"
finding it depends on. G4's assertions pin the measured values from the prior session
(`dw(300)` in `(1e-8, 1e-6)` for both schedules, `|E_A'-E_B'| < 1e-4`) as regression bounds —
both held on this run, confirming those numbers reproduce on independent hardware, not just the
container that originally produced them.

Also re-ran `tests/test_nbn_dmrg_reference_spec.py` (own process) since it shares the same
vendored chk and orbital-loading path and its G1-G3b gates are the ones §2/§7's stale text
referred to:

```
$ uv run pytest -q tests/test_nbn_dmrg_reference_spec.py
....                                                                     [100%]
4 passed in 355.03s (0:05:55)
```

**Did not re-run** the D≤1200 headline schedules (`A`/`B`) — that's the "expensive, only if it
proves hard" case the bead explicitly says to treat as a separate follow-on, and it isn't hard:
the D≤300 gate above already reproduces the inconclusive-at-cheap-D / soft-at-converged-D
finding, and the converged-D numbers themselves (already in `results/nbn_low_spin/runs.jsonl`,
chem-g1i, 2026-09-26) aren't in question — nothing about them changed, I just didn't re-spend
~15-30 min reproducing a number nobody disputed. Did not re-run
`tests/test_nbn_czw_lower_minimum_spec.py` — no code was touched that it depends on (only
markdown), and it wasn't part of this bead's specific ask.

### 2b. Fixed the stale text that caused the duplicate filing

`specs/SPEC_nbn_dmrg_reference.md`:
- §2 (line ~80): replaced "A hard multireference TM benchmark needs that low-spin sector at real
  bond dimension or a larger cluster — a follow-up, out of scope here" with an UPDATE paragraph
  pointing to `SPEC_nbn_low_spin.md`'s resolution and this session's re-verification.
- §7 (Out of scope): the low-spin-sector bullet now says it *was* checked and found soft,
  narrowing "out of scope" to larger clusters/bases (which remain genuinely unexplored).
- §8 R1: updated to record that the low-spin variant does **not** need real bond dimension after
  all — it's soft, just closer to the `dweight` regime boundary than the high-spin sector.

`specs/BACKLOG.md`:
- Moved the open item ("Is NbN's flagged 'hard multireference benchmark' actually hard?") out of
  the open list into `## Done`, with a full explanation of why it was stale (answered by a
  differently-named bead chain before this bead was filed), the fresh reproduction numbers from
  §2a above, and a pointer to this handoff.

These are documentation-only changes — no `.py` files touched, so no lint/test regressions are
possible from them. Confirmed with `git diff --stat` that only the two spec markdown files
changed (plus the new handoff file).

---

## 3. Gate-run evidence (verbatim), this session, this container

```
$ uv run python -c "from hybrid_quantum_solver.dmrg_reference import dmrg_available; print(dmrg_available())"
True          # after `uv sync --extra dmrg --extra test`

$ uv run pytest -q tests/test_nbn_low_spin_spec.py -v
tests/test_nbn_low_spin_spec.py ..                                       [100%]
2 passed in 492.55s (0:08:12)

$ uv run pytest -q tests/test_nbn_dmrg_reference_spec.py
....                                                                     [100%]
4 passed in 355.03s (0:05:55)
```

Hardware: `uname -a` → `Linux 337f99722624 7.0.14-linuxkit #1 SMP PREEMPT ... x86_64 GNU/Linux`;
`nproc` → 8; `/proc/meminfo` MemTotal ≈ 16.35 GB. Well within the compute envelope
(`--stack-mem-gb` used by `SCHEDULES["A'"]`/`["B'"]` is 1 GB each, `<=6`).

No qiskit-aer or pyscf-then-block2 process-sharing occurred: both files ran as standalone
`pytest -q tests/<file>` invocations, one at a time, per `CLAUDE.md`'s isolation rule — neither
imports `qiskit`.

---

## 4. Every number, with regenerating command and vendor location

| Number | Value | Command | Vendored at |
|---|---|---|---|
| dw(D=300), schedule A′, (7,7) | ≈1.05e-7 (in (1e-8, 1e-6), reconfirmed this session) | `uv run pytest -q tests/test_nbn_low_spin_spec.py::test_G4_low_spin_cheap_D_is_pre_registered_inconclusive` | `specs/SPEC_nbn_low_spin.md` §3 Table 2 (original measurement, chem-g1i, 2026-09-26); regression-pinned by G4 |
| dw(D=300), schedule B′, (7,7) | ≈1.05e-7 (in (1e-8, 1e-6)) | same | same |
| \|E_A′ − E_B′\|, (7,7), D≤300 | ≈4.18e-6 Ha (< 1e-4 Ha) | same | same |
| dw(D=1200), schedules A/B, (7,7) | ≈3e-10 (converged) | not re-run this session (see §2a); prior run: `results/nbn_low_spin/runs.jsonl` | `specs/SPEC_nbn_low_spin.md` §3 Table 2 |
| \|E_A − E_B\|, (7,7), D≤1200 | 2.91e-10 Ha | not re-run this session | same |
| S=1 vs S=3 (committed) gap | 1.538 mHa (DMRG A′), 1.537 mHa (exact FCI) | `uv run pytest -q tests/test_nbn_low_spin_spec.py::test_G5_committed_sector_is_not_the_cas_ground` (reconfirmed this session) | `specs/SPEC_nbn_low_spin.md` §3 Table 1 |

No new physics numbers were produced this session — G4/G5's own assertions (which encode the
numbers above as bounds) are the reproduction evidence.

---

## 5. Pre-registered criteria and pass/fail

From the bead's acceptance criteria:

1. **"All four existing SCHEDULES run against nelec=(7,7) forced, each in its own process...
   capped at D≤300 and --stack-mem-gb≤6"** — **PASS, via the pre-existing work this session
   re-verified.** `A'`/`B'` (D≤300, ≤1GB stack) are what `test_nbn_low_spin_spec.py::test_G4_*`
   runs and what I re-ran fresh this session. `A`/`B` (D≤1200) also exist, already run
   (chem-g1i) — going past D=300 there was the deliberate, bead-anticipated "if hard, treat as
   separate follow-on" path, except it turned out to *resolve* the D≤300 ambiguity rather than
   confirm hardness, so it was folded into the same spec rather than filed as a new bead.
2. **"dw(D=300) and cross-schedule spread |E_A'-E_B'| reported"** — **PASS.** §4 above; source
   `SPEC_nbn_low_spin.md` §3 Table 2.
3. **"The easy-vs-hard verdict is recorded"** — **PASS: easy/soft.** At D≤300 neither KILL nor
   CONFIRM pre-registered condition fires (honest gap in the two-sided design); at converged
   D≤1200 the sector's own KILL threshold is met (`dw≈3e-10 ≪ 1e-7`, spread `≈3e-10 Ha ≪ 1e-5
   Ha`) — the "hard multireference benchmark" framing does not survive.
4. **"...and if hard, whether the 3.5 mHa sector gap survives at the largest converged D actually
   reached"** — **N/A, not hard**, but checked anyway (it's how "not hard" was established): the
   gap is 5.066 mHa at converged D (not 3.5 — the 3.5 mHa figure was itself from the stale
   §2/§7 text, superseded by the real-geometry numbers in `SPEC_nbn_low_spin.md`), and it does
   not invert or shrink toward zero as D increases from 300→1200.
5. **"...with any stalled-vs-converged ambiguity explicitly checked, given chem-4e9/chem-mjz"**
   — **PASS.** dw decreases monotonically with D (schedule A: 8.32e-8 at D=100/200/300-era →
   6.44e-9 → 8.20e-10 → ... down to 3.06e-10 at D=1200 headline), the signature of genuine
   convergence, not a stall.

---

## 6. What I decided not to do, and why

- **Did not re-run the D≤1200 headline schedules.** They cost ~9-17 min each and answer a
  number nobody disputes (chem-g1i's 2026-09-26 measurement, `results/nbn_low_spin/runs.jsonl`,
  cross-checked independently by chem-czw on a *different* UHF minimum with the same qualitative
  result). Re-running them would have been re-deriving, not re-verifying, and this repo's own
  convention (`SPEC_nbn_dmrg_reference.md`'s own headline record) is not to re-run multi-minute
  driver-level numbers on every cycle.
- **Did not touch `SPEC_nbn_low_spin.md` itself.** It already correctly and thoroughly documents
  this exact finding (§0, §3 Table 2, G4/G5) — nothing there was wrong or stale, only the
  *other* spec (`SPEC_nbn_dmrg_reference.md`) and the backlog entry pointed away from it.
- **Did not re-run `tests/test_nbn_czw_lower_minimum_spec.py`.** No code it depends on was
  touched; it's unrelated to this bead's specific ask (it's about which UHF minimum to vendor,
  not about sector hardness).
- **Did not investigate a genuinely larger cluster/basis** as a harder multireference benchmark
  — explicitly out of scope per the bead text ("treat as a separate, larger follow-on") and per
  `SPEC_nbn_dmrg_reference.md` §7 (unchanged scope boundary, only reworded).

---

## 7. What could not be verified

- **Cross-host generality of the D≤300 near-miss numbers.** I reproduced them bit-approximately
  (within the regression bounds `test_G4_*` pins) on one container (x86_64, 8-core, OpenBLAS via
  this repo's `scipy-openblas64`), not on the original container that produced the exact
  digits quoted in `SPEC_nbn_low_spin.md` §3 Table 2. The regression bounds are deliberately
  order-of-magnitude, not exact-digit, so this is expected and by design (per that spec's own
  G4 docstring), not a gap I introduced.
- **The D≤1200 headline numbers** — not independently re-run this session; relied on the
  existing, differently-attributed record (chem-g1i/chem-czw) rather than re-deriving. Flagged
  explicitly here rather than silently assumed.

---

## 8. Git status / suggested commands (NOT executed — human review required)

Per this session's git policy, no commit/push was made. Files touched by this bead's work:

```
 M specs/BACKLOG.md                  (moved the open item to Done, with explanation)
 M specs/SPEC_nbn_dmrg_reference.md  (§2/§7/§8 R1 — stale "never run"/"out of scope" text fixed)
?? sandbox-handoffs/chem-csf.md
```

(Every other modified/untracked file in `git status` predates this session — unrelated
in-progress work on `certified_dipole*`/`reachability.py`/etc. — not touched here.)

Suggested commands for a human to review and run:

```bash
git add specs/BACKLOG.md specs/SPEC_nbn_dmrg_reference.md sandbox-handoffs/chem-csf.md
git commit -m "$(cat <<'EOF'
nbn: fix stale 'never run' framing for the (7,7) low-spin sector hardness question (chem-csf)

The sector was already run and found soft (chem-dc7/chem-g1i/chem-czw,
SPEC_nbn_low_spin.md), but SPEC_nbn_dmrg_reference.md's own text and the
BACKLOG.md entry still said "never run, never gated" -- which is why this
bead got filed as if the question were open. Re-verified the finding fresh
on new hardware (tests/test_nbn_low_spin_spec.py, tests/test_nbn_dmrg_reference_spec.py,
both green) and corrected the stale cross-references.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
bd close chem-csf --reason="Already answered before this bead was filed, by a differently-worded bead chain (chem-dc7/chem-g1i/chem-czw) documented in SPEC_nbn_low_spin.md Table 2/G4: the (7,7) sector is soft (dw->3e-10 by D=1200, two schedules agree to 0.3 nHa), not a hard multireference benchmark. This session independently re-verified that finding fresh on new hardware (tests/test_nbn_low_spin_spec.py, tests/test_nbn_dmrg_reference_spec.py both green) and fixed the stale SPEC_nbn_dmrg_reference.md text and BACKLOG.md entry that caused the duplicate filing."
```
