# chem-ivu handoff — re-run chem-czw's lower-minimum NbN table on Linux with the fixed `_tight_scf` loop

**Bead:** chem-ivu (P3, task) — confirm (or re-derive) chem-czw's Linux lower-minimum SCF energy
and S=3 FCI now that `benchmark_nbn._tight_scf`'s R3 loop-body bug (`mf.kernel(dm0=dm)` dedented
out of the stability loop) is fixed.

**Resolution: confirmed, unchanged.** The fixed loop reproduces both pinned numbers to within
noise (SCF: 5.68e-14 Ha; S=3 FCI: bit-identical). No re-derivation needed.

---

## 1. Background

`sandbox-handoffs/chem-czw.md`'s ADDENDUM (local macOS re-verification, same session) found that
chem-y6w's rewrite of `_tight_scf` had a bug: the re-convergence call `mf.kernel(dm0=dm)` had been
dedented out of the stability-following loop, so it never ran after an unstable verdict (and threw
`UnboundLocalError` whenever the *first* stability check passed, which happens on every macOS run).
The bug was fixed in the working tree (`benchmark_nbn.py`, already uncommitted before this
session — see `git diff benchmark_nbn.py`). But chem-czw's Linux lower-minimum numbers — the SCF
energy `-110.02813374576796` pinned in `tests/test_nbn_scf_determinism_spec.py` G1/G2 and
`scripts/gen_nbn_lower_minimum_chk.py`, and the four-sector FCI/DMRG table in
`tests/test_nbn_czw_lower_minimum_spec.py` — had been produced by the **buggy** loop and were
flagged unverified against the fix (`specs/SPEC_nbn_scf_determinism.md` §8 R3, as written before
this session).

This bead's job: rerun `scripts/gen_nbn_lower_minimum_chk.py` with the fixed loop on the Linux
container, confirm it still lands on `-110.02813374576796` and reproduces the cheap S=3 FCI. If
not, re-derive.

## 2. What I ran

Host: `796434596efe` (from `uname -a`: `Linux 796434596efe 7.0.12-linuxkit #1 SMP PREEMPT ... x86_64`),
8 cores (`nproc`), 16 GB RAM (`/proc/meminfo`), no GPU.

```
$ uv sync --extra dmrg --extra test --extra dev   # combined, per chem-czw's noted gotcha
$ time uv run python scripts/gen_nbn_lower_minimum_chk.py /tmp/nbn_scf_lower_regen.chk
...
<class 'pyscf.soscf.newton_ah.SecondOrderUHF'> wavefunction is stable in the internal stability analysis
[CLASSICAL] tightened SCF: e_tot=-110.02813375 Ha  converged=True  <S^2>=12.0779
e_tot = np.float64(-110.0281337457679)
OK: matches the deterministic lower-energy minimum (diff 5.684e-14 Ha), written to /tmp/nbn_scf_lower_regen.chk
real    0m9.024s
```

```
$ time uv run python nbn_low_spin.py fci --nelec 10 4 --chk /tmp/nbn_scf_lower_regen.chk
RECORD {"kind": "fci", "nelec": [10, 4], "twos_pinned": null, "root": 0,
        "energy": -110.04119080464493, "s2": 12.000000000270802, "ndet": 1002001,
        "wall_s": 5.8, "chk": "/tmp/nbn_scf_lower_regen.chk", "host": "796434596efe",
        "cpu": "x86_64", "ncpu": 8, "stamp": "2026-09-27T22:59:33"}
real    0m7.790s
```

Both numbers match the pinned literals:
- SCF: `-110.0281337457679` vs pinned `-110.02813374576796` — diff `5.68e-14` Ha (well inside the
  script's own `1e-9` Ha tolerance, and 5 orders of magnitude under the 0.548 mHa gap to the other
  minimum).
- S=3 FCI: `-110.04119080464493` — **bit-identical** to
  `tests/test_nbn_czw_lower_minimum_spec.py::LOWER_MINIMUM_FCI[(10, 4)]`.

Total wall time for both: ~17s, one process (pyscf+qiskit only, no block2 — safe to combine).

## 3. Why the bug didn't corrupt these numbers

`mf.stability()` reports this minimum as internally stable on its **first** check (see the printed
`"wavefunction is stable in the internal stability analysis"` line above, reached after exactly one
prior `"has an internal instability"` cycle). The buggy re-convergence line
(`mf = mf.newton(); mf.kernel(dm0=dm)`) only executes in the branch taken *after* an unstable
verdict, to re-converge from the followed-instability density matrix. Since the Linux path only
needed that branch once (to get from the vendored-style saddle down to this minimum) and then hit
"stable" on the very next check — before the loop would exit anyway — the missing
re-convergence call was never on the execution path that determines the final `mf.e_tot`. The bug
was latent on this geometry/minimum, not triggered. (On macOS, by contrast, the very first
stability check already passes, so the loop body — including the buggy line, dedented to
unconditionally execute after the `try` block regardless of branch — hits the `UnboundLocalError`
immediately; that's a different failure mode than "silently wrong energy", which is why chem-czw's
addendum caught it as a crash, not a silent number change.)

## 4. Gate-run evidence (verbatim)

```
$ uv run python -m pytest tests/test_nbn_scf_determinism_spec.py tests/test_nbn_czw_lower_minimum_spec.py -v
collected 9 items
tests/test_nbn_scf_determinism_spec.py::test_G1_deterministic_across_ambient_thread_counts PASSED [ 11%]
tests/test_nbn_scf_determinism_spec.py::test_G2_deterministic_across_repeats PASSED [ 22%]
tests/test_nbn_czw_lower_minimum_spec.py::test_lower_minimum_ordering_matches_vendored_S1_lowest PASSED [ 33%]
tests/test_nbn_czw_lower_minimum_spec.py::test_lower_minimum_dmrg_schedules_agree PASSED [ 44%]
tests/test_nbn_czw_lower_minimum_spec.py::test_vendored_chk_gives_lower_casci_energy_than_the_other_minimum[nelec0] PASSED [ 55%]
tests/test_nbn_czw_lower_minimum_spec.py::test_vendored_chk_gives_lower_casci_energy_than_the_other_minimum[nelec1] PASSED [ 66%]
tests/test_nbn_czw_lower_minimum_spec.py::test_vendored_chk_gives_lower_casci_energy_than_the_other_minimum[nelec2] PASSED [ 77%]
tests/test_nbn_czw_lower_minimum_spec.py::test_vendored_chk_gives_lower_dmrg_energy_for_S0_too PASSED [ 88%]
tests/test_nbn_czw_lower_minimum_spec.py::test_lower_minimum_chk_regenerates_and_reproduces_S3_FCI PASSED [100%]
9 passed, 60 warnings in 32.43s
```

`test_lower_minimum_chk_regenerates_and_reproduces_S3_FCI` is itself a second, independent
regeneration (its own tmp chk, run inside pytest) that separately confirms the same SCF energy and
S=3 FCI — so this bead's manual run above and the gate's own fixture-gated sub-test agree with each
other as well as with the pre-existing pinned literals.

```
$ uv run ruff check scripts/gen_nbn_lower_minimum_chk.py benchmark_nbn.py \
    tests/test_nbn_scf_determinism_spec.py tests/test_nbn_czw_lower_minimum_spec.py
All checks passed!
```

Not re-run: `tests/test_nbn_low_spin_spec.py` (the full four-sector battery's own gate, ~8.5 min) —
out of scope. This bead's pre-registered check was the SCF energy plus the *cheap* S=3 FCI only
(chem-ivu's description explicitly names both); nothing in that check touches S=1/S=2/(7,7)-DMRG,
and confirming the SCF minimum plus one sector is sufficient to show the buggy loop's specific
defect (missing re-convergence after an unstable verdict) was never exercised on this path.

## 5. Numbers produced, with regenerating command and vendor location

| Number | Value | Command | Vendor location |
|---|---|---|---|
| Lower minimum SCF energy (reconfirmed, fixed loop) | `-110.0281337457679 Ha` (diff 5.68e-14 Ha from pinned) | `uv run python scripts/gen_nbn_lower_minimum_chk.py /tmp/out.chk` | not vendored (regenerated fresh; pinned literal already in `tests/test_nbn_scf_determinism_spec.py`/`scripts/gen_nbn_lower_minimum_chk.py`) |
| Lower minimum, S=3 (10,4) exact FCI (reconfirmed, fixed loop) | `-110.04119080464493 Ha` (bit-identical to pinned) | `nbn_low_spin.py fci --nelec 10 4 --chk /tmp/out.chk` | appended row in `results/nbn_low_spin/runs.jsonl` (this session); pinned literal already in `tests/test_nbn_czw_lower_minimum_spec.py::LOWER_MINIMUM_FCI[(10, 4)]` |

## 6. Pre-registered criteria and pass/fail

Bead's acceptance criterion: "Linux regen with the fixed loop reproduces LOWER_MINIMUM_E_SCF and
the S=3 FCI to 1e-6 Ha, or the table is re-derived and the gate constants updated."

- **SCF energy within 1e-6 Ha of `-110.02813374576796`** — **PASS.** Diff `5.68e-14` Ha.
- **S=3 FCI within 1e-6 Ha of `-110.04119080464493`** — **PASS.** Bit-identical (diff 0).
- **Table re-derivation** — **N/A, not triggered** (criterion is an "or"; the first branch held).
- **Test suite (scoped gates) passes** — **PASS.** 9/9 (§4). Ruff clean on touched/new files (§4).

## 7. What I changed

- `results/nbn_low_spin/runs.jsonl` — one appended line (this session's S=3 FCI confirmation run;
  pure append, no existing lines touched — same convention chem-czw used).
- `specs/SPEC_nbn_scf_determinism.md` §8 R3 — appended the confirmation (fixed loop reproduces both
  pinned numbers on this Linux container, plus the "why the bug was latent, not triggered, here"
  explanation from §3 above).
- `specs/BACKLOG.md` — new `[x]` Done entry (top of `## Done`) recording this bead's finding.
- `sandbox-handoffs/chem-ivu.md` — this file.

No code changes beyond what was already in the working tree before this session (the `_tight_scf`
fix in `benchmark_nbn.py` was pre-existing, uncommitted work from the prior session that produced
the macOS addendum — not touched further here).

## 8. What I decided not to do, and why

- **Did not re-run the full `tests/test_nbn_low_spin_spec.py` four-sector battery** (~8.5 min) —
  out of scope; the bead's own acceptance criterion names only the SCF energy and the cheap S=3
  FCI, and confirming those is sufficient to show the specific defect (missing re-convergence after
  an *unstable* stability verdict) was never on this geometry/minimum's execution path (§3).
- **Did not re-derive S=1/S=2/(7,7)-DMRG rows on a fresh fixed-loop chk** — same reason; would
  repeat ~25 minutes of compute the bead didn't ask for and the mechanism argument (§3) already
  covers.
- **Did not keep the regenerated `/tmp/nbn_scf_lower_regen.chk`** — `/tmp` is outside the repo and
  was deleted after use; deterministically regenerable from tracked inputs
  (`scripts/gen_nbn_lower_minimum_chk.py` + `specs/nbn_mp-2634.cif`), same convention chem-czw
  established for this artifact.
- **Did not touch any of the other pre-existing dirty files in the working tree**
  (`nbn_dmrg_reference.py`, `df_factorization.py`, `shift_both_sides.py`,
  `hybrid_quantum_solver/dmrg_reference.py`, their specs/tests, `.claude/settings.local.json`, the
  other `sandbox-handoffs/*.md` files, `specs/SPEC_nbn_scf_determinism.md`'s pre-existing content,
  `scripts/gen_nbn_lower_minimum_chk.py`, `tests/test_dmrg_scratch_isolation_spec.py`,
  `tests/test_nbn_scf_determinism_spec.py`) — all predate this session (other beads' uncommitted
  work in this shared container); confirmed via `git status` before and after.

## 9. What could not be verified

- Whether the mechanism argument in §3 (the buggy line never executes on this path because the
  loop only needed one unstable→stable transition) generalizes to *any* other geometry/minimum on
  this molecule, or is specific to this exact SCF trajectory. Not investigated further — out of
  scope for a confirmation bead about this specific pinned table.
- Cross-BLAS-backend / cross-host generality of "the bug was latent, not triggered" — verified only
  on this one Linux container, same scope limitation `SPEC_nbn_scf_determinism.md` R1 already
  states for the fix itself.

## 10. Git status / suggested commands (NOT executed — human review required)

Per this session's git policy, no commit/push was made.

```bash
git add results/nbn_low_spin/runs.jsonl specs/SPEC_nbn_scf_determinism.md specs/BACKLOG.md \
        sandbox-handoffs/chem-ivu.md
git commit -m "chem-ivu: confirm chem-czw's lower-minimum NbN SCF+S3-FCI survive the fixed _tight_scf loop

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
bd close chem-ivu --reason="Re-ran scripts/gen_nbn_lower_minimum_chk.py with the fixed _tight_scf loop on the Linux container: reproduces the pinned SCF energy (-110.0281337457679, diff 5.68e-14 Ha) and the pinned S=3 FCI (-110.04119080464493, bit-identical). The buggy loop's missing re-convergence line was never on this geometry's execution path (mf.stability() found it stable after a single unstable->stable transition), so it was latent here, not triggered -- no re-derivation needed. Gates green: test_nbn_scf_determinism_spec.py (2/2), test_nbn_czw_lower_minimum_spec.py (7/7). specs/SPEC_nbn_scf_determinism.md R3 and specs/BACKLOG.md updated with the confirmation."
```

(Note: all the *other* modified/untracked files `git status` shows predate this session and are
untouched by chem-ivu — see §8's list. A human `git add`/commit of this bead's work should only
stage the four files listed above.)

All gates referenced above are confirmed green as of this handoff (§4/§6).
