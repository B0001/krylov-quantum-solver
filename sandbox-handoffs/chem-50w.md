# chem-50w handoff — pin `pyscf.lib.num_threads(1)` around `reachability._driver_orbitals`'s own RHF solve

**Bead:** chem-50w (P2, bug) — "`reachability._driver_orbitals`'s unpinned `drv.run_pyscf()` shares
chem-7cx/chem-s6c's thread-nondeterminism mechanism."

**Resolution:** implemented and gated. `_driver_orbitals`'s one `drv.run_pyscf()` call is now pinned
to `pyscf.lib.num_threads(1)` (ambient count restored in a `finally`), the same pattern
chem-7cx/chem-s6c/chem-0yh used elsewhere for the identical mechanism. A new spec
(`specs/SPEC_driver_orbitals_scf_determinism.md`) and gate file
(`tests/test_driver_orbitals_scf_determinism_spec.py`, G1-G5) pin the claim and the explicit
before/after retry-count measurement the bead and chem-s6c's notes asked for.

---

## 1. What changed

### 1a. Root cause (confirmed directly, not inferred)

`reachability._driver_orbitals` re-runs an unsymmetrized (`symmetry=False`, PySCFDriver hardcodes
this) RHF solve to verify a rebuild reproduces `mh.qubit_hamiltonian`. PySCF's own multi-threaded
OpenMP integral evaluation is not bit-associative, so this independent solve can land on a
different locally-stable branch depending only on ambient thread count — the same mechanism
chem-7cx fixed in `build_molecular_hamiltonian` and chem-s6c fixed in
`scf_symmetry_status`/`hf_orbital_irreps`, but left untouched here because `_driver_orbitals`
already carries chem-g14's bounded retry + sign-gauge-search machinery for a related, distinct
mechanism.

Direct repro (square H4, sto-3g, `mh` built once and held FIXED under `num_threads(1)` to isolate
this from `build_molecular_hamiltonian`'s own separate, still-unpinned-on-this-branch
nondeterminism — see §1c below):

```
$ for t in 1 2 4 8; do uv run python /tmp/repro_counts3.py $t 1.19 1e-6; done   # before the fix
ambient_at_call=1 a=1.19 ct=1e-06 found=yes calls=1
ambient_at_call=2 a=1.19 ct=1e-06 found=NO  calls=5
ambient_at_call=4 a=1.19 ct=1e-06 found=yes calls=4
ambient_at_call=8 a=1.19 ct=1e-06 found=NO  calls=5
```

and, holding ambient threads fixed at 2 and calling the unpinned `_driver_orbitals(mh)` ten times in
a row (`a = 1.35`, `conv_tol = 1e-9`):

```
found:  yes yes yes yes  NO  NO yes yes yes yes
calls:    3   1   3   1   5   5   2   2   5   5
```

Both the classification (found/not-found) and the retry count vary run to run on an *identical*
input — exactly the mechanism the bead named. (Throwaway repro scripts, not committed; the spec's
§1/§9 reproduce the same table via the committed gate file's own subprocess machinery.)

### 1b. Fix

`reachability.py`, inside `_driver_orbitals`'s retry loop:

```python
ambient_threads = lib.num_threads()
lib.num_threads(1)
try:
    drv.run_pyscf()
finally:
    lib.num_threads(ambient_threads)
```

17 lines changed total (`git diff reachability.py`: +18/-2, one new top-level import
`from pyscf import lib`). The sign-gauge search itself makes no further SCF call (`rebuild_from`
only re-runs `drv.to_problem()` on already-computed `mo_coeff` — a fresh AO2MO transform, not a new
RHF solve), so nothing else in the function needed touching.

### 1c. New spec + gate: `specs/SPEC_driver_orbitals_scf_determinism.md`,
`tests/test_driver_orbitals_scf_determinism_spec.py`

Every gate builds its own `mh` under `pyscf.lib.num_threads(1)` **before** varying ambient threads
only around the call under test. This is a deliberate choice, explained in the spec (§2/§8 R1): it
isolates "does `_driver_orbitals`'s own RHF call add nondeterminism on top of a fixed comparison
target" (this bead's actual question) from "is the comparison target itself reproducible"
(`build_molecular_hamiltonian`'s own free RHF call — chem-7cx's question, and chem-7cx's fix lives
on a sibling sandbox branch **not merged into `sandbox/chem-50w`**, which was cut from the same
point on `main` chem-7cx/chem-s6c/chem-0yh were). Building `mh` under `num_threads(1)` is itself
deterministic — confirmed directly (8 repeated builds of the same witness gave a bit-identical
`qubit_hamiltonian`/`hf_energy`) — so this is a legitimate, honest test fixture, not a tautology.

Gates, pyscf + qiskit-nature only (imports `reachability` and
`hybrid_quantum_solver.molecular_hamiltonian`; **not** block2 — own file, per the `test_*_spec.py`
process-isolation convention):

- **G1 — cross-ambient-thread determinism, fixed `mh` (DEFINITION OF DONE part 1).** For each of
  three witnesses (`a=1.19,ct=1e-6`; `a=1.10,ct=1e-6`; `a=1.35,ct=1e-9`), `_driver_orbitals(mh)`
  called once per ambient thread count in `{1,2,4,8}` (each its own subprocess) gives the same
  found/MO-coefficient outcome across all four.
- **G2 — repeat-call determinism, one fixed ambient setting.** Ten fresh subprocess calls, ambient
  fixed at 2, are pairwise bit-identical.
- **G3 — same-process repeat-call determinism + thread-count restored.** Ten fresh in-process calls
  are pairwise bit-identical and leave `lib.num_threads()` unchanged after every call.
- **G4 — retry loop stops firing from thread noise alone (DEFINITION OF DONE part 2 — the explicit
  before/after retry-count measurement the bead asked for).** For the fixed `mh`, across all four
  ambient thread counts, `_driver_orbitals` makes exactly **one** `PySCFDriver.run_pyscf` call —
  down from the 1-5 observed before the fix (§1b's table).
- **G5 — no regression on `SPEC_rebuild_gauge_search.md`'s own mechanics.** Re-runs
  `tests/test_rebuild_gauge_search_spec.py` (G3-G5: cheap-when-unneeded call counting, energy
  pre-filter retry-not-false-accept, oversized-active-space refusal) as a subprocess check.

## 2. Gate-run lines (verbatim)

```
$ uv run pytest -q tests/test_driver_orbitals_scf_determinism_spec.py
.....
5 passed, 1 warning in 76.57s (0:01:16)

(rerun for stability:)
$ uv run pytest -q tests/test_driver_orbitals_scf_determinism_spec.py
.....
5 passed, 1 warning in 76.80s (0:01:16)
```

```
$ uv run pytest -q tests/test_rebuild_gauge_search_spec.py
.....
5 passed, 5 warnings in 2.27s
```

```
$ uv run ruff check .
All checks passed!
```

## 3. Regression check — pre-existing flakiness unaffected (not caused, not cured)

These three files call through `reachability.py` but their known flakiness is independently driven
by `build_molecular_hamiltonian`'s own still-unpinned free RHF call (chem-7cx's territory, not
merged on this branch) — confirmed by running each **before** (`git stash` on just `reachability.py`)
and **after** this fix and finding the same failure pattern/magnitude both times:

```
tests/test_symmetry_reachability_spec.py:   5 failed, 9 passed   -- IDENTICAL set of 5, before AND after
tests/test_reachability_tolerance_spec.py:  5 failed, 7 passed   -- before AND after (different specific
                                                                     IDs each run even within one state --
                                                                     pre-existing flakiness, per
                                                                     SPEC_rebuild_gauge_search.md's own
                                                                     prior finding on this file)
tests/test_eigenstate_reachability_spec.py (49 cases, flaky; 3x baseline + 2x after-fix runs):
  before: 14 failed/35 passed, 12 failed/37 passed, 15 failed/34 passed
  after:  17 failed/32 passed, 16 failed/33 passed
```

The after-fix counts sit in the same range as baseline's (12-15 before vs 16-17 after, both well
within this file's own documented baseline range of 14-18/45 pre-chem-g14 and 14-15/45
post-chem-g14 — see `SPEC_rebuild_gauge_search.md` §5 G2). This fix does not clear
`test_G2_decision_equals_symmetry_adapted_fci[1.19]` (still fails both before and after — confirmed
directly, isolated single-case rerun) because that specific flake is a *combination* of this
mechanism and `build_molecular_hamiltonian`'s own (chem-7cx's, unmerged here); fixing only one side
of a two-sided nondeterminism does not make the pair deterministic. This is the explicit, measured
scope boundary the spec's §2/§7/§8-R1/§9-R2 record — not a surprise found after the fact.

## 4. Every number produced, and where it is vendored

| Number | Command | Vendored in |
|---|---|---|
| Before-fix 10-repeat table (found/calls) | throwaway `/tmp` script, not committed; reproduced by the committed gate's own machinery | `specs/SPEC_driver_orbitals_scf_determinism.md` §1, §9 |
| Before-fix cross-thread table (a=1.19) | same | §1, §9 |
| After-fix: G1/G2/G4 all-pass, 1 call every time | `tests/test_driver_orbitals_scf_determinism_spec.py` (committed, reruns the claim itself) | §9, and this handoff §2 |
| Unisolated (R2) residual gap | throwaway `/tmp` script; the gap itself is NOT this bead's to close | §8 R1, §9 R2 |
| Regression counts (§3 above) | `uv run pytest -q <file>`, run 2-3x each, before/after via `git stash` | this handoff §3 only (not re-vendored into a spec — these are pre-existing numbers that belong to chem-7cx's/chem-g14's specs, not a new claim of this bead's) |

## 5. Pre-registered criteria vs. outcome

| Criterion (spec §5) | Outcome |
|---|---|
| G1 — cross-ambient-thread determinism, fixed `mh`, 3 witnesses | **PASS** |
| G2 — repeat-call determinism, one ambient setting | **PASS** |
| G3 — same-process repeat + thread-count restored | **PASS** |
| G4 — retry loop settles to 1 call across all ambient settings | **PASS** (down from 1-5 before) |
| G5 — no regression on `SPEC_rebuild_gauge_search.md` mechanics | **PASS** |

G1 and G4 were pre-registered as the definition of "done" — both pass.

## 6. What I decided not to do, and why

- **Did not touch `build_molecular_hamiltonian`'s own unpinned RHF call.** That is chem-7cx's bead,
  on a sibling sandbox branch not merged into `sandbox/chem-50w`. Pulling it in here would silently
  widen this bead's scope and make this spec's claim ambiguous (determinism from *which* fix?). The
  spec's §7/§8 R1/§9 R2 instead measure and name the residual gap explicitly.
- **Did not attempt to clear `test_eigenstate_reachability_spec.py::test_G2[1.19]` or any other
  pre-existing failing case.** Confirmed (§3) that doing so requires chem-7cx's fix too; declaring
  victory on a two-sided nondeterminism by fixing only one side would be exactly the kind of
  unearned confidence this repo's culture (CLAUDE.md) rules out.
- **Did not raise `_MAX_SIGN_SEARCH_ORBITALS`/`_MAX_REBUILD_RETRIES`** or otherwise touch
  `SPEC_rebuild_gauge_search.md`'s own mechanics — out of scope, and G5 confirms no regression there.
- **Did not run the full `make gates`/`scripts/run_gates.sh` sweep** (block2/DMRG gates are
  unrelated to this change and expensive — ~8 CPUs/16GB shared envelope). Ran the specific files
  that import/exercise `reachability.py` instead (§2/§3). If a reviewer wants the full sweep, it is
  `scripts/run_gates.sh` or `make gates`; nothing in this diff touches DMRG/block2 code paths.

## 7. What I could not verify

- **Cross-host/cross-pyscf-version reproduction.** Verified only on this container (Linux x86_64, 8
  cores, OpenBLAS, pyscf 2.13.1, `uname -a`:
  `Linux 5404a8aa3661 7.0.14-linuxkit #1 SMP PREEMPT Fri Sep 18 10:19:32 UTC 2026 x86_64 GNU/Linux`).
  Not tested on macOS or a different pyscf version — same caveat every prior spec in this family
  (`SPEC_reachability_scf_determinism.md` R3, `SPEC_scf_branch_determinism.md` R1) records.
- **Whether the now-deterministic rebuild branch is "physically correct."** Explicitly out of scope
  (spec §2) — this fix makes the result stable, not necessarily right; chem-g14's own
  cross-platform/version branch-correctness divergence is untouched.

## 8. Timing / hardware

Full gate file (`tests/test_driver_orbitals_scf_determinism_spec.py`, G1-G5, includes 2 full
`pytest` subprocess re-invocations of `test_rebuild_gauge_search_spec.py` inside G5 plus ~40
subprocess `PySCFDriver` builds across G1/G2/G4): **~77s** wall time, two independent runs (76.57s,
76.80s). `test_rebuild_gauge_search_spec.py` alone: 2.3s. `test_eigenstate_reachability_spec.py`
(49 cases, unrelated regression check): ~115-127s per run, 5 runs total (~10 min). All well inside
the single-worker 8-CPU/16GB envelope; no DMRG/block2 code touched, no `--stack-mem-gb` relevant.

## 9. Git status at handoff

```
$ git status
On branch sandbox/chem-50w
Changes not staged for commit:
	modified:   .beads/issues.jsonl       (bd claim)
	modified:   reachability.py           (the fix, +18/-2)
Untracked files:
	.claude/settings.local.json           (not mine to touch; pre-existing sandbox artifact)
	specs/SPEC_driver_orbitals_scf_determinism.md
	tests/test_driver_orbitals_scf_determinism_spec.py
	sandbox-handoffs/chem-50w.md
```

Per this session's git policy, the tree is left **uncommitted**. Suggested commit (one commit; `git
add reachability.py specs/SPEC_driver_orbitals_scf_determinism.md
tests/test_driver_orbitals_scf_determinism_spec.py sandbox-handoffs/chem-50w.md .beads/issues.jsonl`
then commit) — left for the human to run, per the git-on-this-run instructions in the task prompt.

## 10. Beads

- **chem-50w**: close, as DONE — the candidate fix in the bead description is implemented, gated,
  and measured exactly as the bead's "Candidate fix" paragraph asked (before/after retry-count
  measurement included).
- No new follow-up bead filed: the residual (unisolated) gap is entirely chem-7cx's (already filed,
  already has its own spec in progress on a sibling branch), not new scope this bead's investigation
  uncovered.
