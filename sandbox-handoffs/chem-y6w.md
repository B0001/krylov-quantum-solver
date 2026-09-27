# chem-y6w handoff — NbN SCF reference: metastability + determinism

**Bead:** chem-y6w (P2, bug) — "Committed NbN SCF reference (data/nbn_scf.chk, S=3 sector) sits at
a threading-nondeterministic, marginally-unstable UHF solution."

**Resolution:** acceptance-criteria path **(b)** — explicit spec-recorded decision to keep the
current metastable chk as the reference, with caveats on every dependent gate/table, plus a
determinism-forcing fix so future regenerations are reproducible rather than a coin flip.

---

## 1. What changed

### 1a. Root cause (re-confirmed, not new — established by chem-owr)

The committed `specs/nbn_scf_reference.chk` (10,4)/S=3 UHF solution sits essentially exactly on a
pyscf internal-stability-Hessian zero crossing. Two genuinely, separately-locally-stable UHF minima
exist: the committed one, and another 0.548 mHa lower in SCF energy (4.84 mHa lower after
CASCI(14,14)/FCI(10,4) — about 3× the pinned S1-S3 gap of 1.538 mHa). Which one a *fresh*
`benchmark_nbn._tight_scf` regeneration finds used to depend on ambient multi-threaded-BLAS
floating-point non-associativity inside pyscf's Davidson-based `mf.stability()` (chem-owr measured
3/4 default-thread trials finding the lower minimum, 1/4 reproducing the committed one).

### 1b. Fix: pin pyscf's own thread count for the duration of `_tight_scf`

`benchmark_nbn.py::_tight_scf` now does:

```python
ambient_threads = lib.num_threads()
lib.num_threads(1)
try:
    ... (unchanged SCF + Newton/SOSCF + stability-following loop) ...
finally:
    lib.num_threads(ambient_threads)
return mf
```

`pyscf.lib.num_threads(n)` is pyscf's own runtime-settable thread control for its compiled OpenMP
extensions — its own docstring explicitly prefers it over `os.environ['OMP_NUM_THREADS']`, which
pyscf reads once at import and cannot be changed afterward (confirmed by testing: setting the env
var mid-process has no effect). Full diff:

```
$ git diff -- benchmark_nbn.py
```
(34 insertions / 15 deletions — wraps the existing kernel/newton/stability loop in a thread pin,
restored in `finally`; no other logic changed. Callers `ground_state_mf` and `nbn_low_spin.py`'s
`cif`/`scan` subcommands are unaffected in interface, only in the determinism of what they return.)

### 1c. New spec + gate: `specs/SPEC_nbn_scf_determinism.md`, `tests/test_nbn_scf_determinism_spec.py`

Two gates, pyscf-only (imports `benchmark_nbn` → `qiskit`/`qiskit_nature`; **not** block2 — own
file per the `test_*_spec.py` process-isolation convention):

- **G1 — cross-ambient-thread determinism.** Two fresh `ground_state_mf` runs on
  `specs/nbn_mp-2634.cif` (distinct chkfiles), `lib.num_threads()` = 8 then 2, must give
  bit-identical `e_tot` (`==`).
- **G2 — repeat-call determinism.** Three fresh runs at ambient thread count must be pairwise
  bit-identical.

**RED before fix** (verified by `git stash push -- benchmark_nbn.py`, rerun, `git stash pop`):
2 failed. **GREEN after fix** (restored):

```
$ UV_PROJECT_ENVIRONMENT=/tmp/venv uv run pytest -q tests/test_nbn_scf_determinism_spec.py
..                                                                       [100%]
2 passed, 50 warnings in 22.86s
```

Measured energy across ambient thread counts 1/2/4/8 and 5 repeated trials at ambient=8 (previously
the flaky setting): bit-identical `-110.02813374576796 Ha` (agrees with chem-owr's own
single-threaded measurement) every time. This is the *lower* minimum, not the vendored chk's
`-110.0275853827266 Ha` — see §3.

### 1d. Decision recorded in specs (not just this handoff)

- `specs/SPEC_nbn_low_spin.md` §0 — new paragraph: "chem-y6w's decision (2026-09-27) — option
  (b): keep the vendored chk, caveat every dependent, fix the regeneration path's determinism,"
  with reasoning (compute cost of full re-derivation vs. CLAUDE.md's own long-running-job
  deferral guidance) and a pointer to the declined-for-now follow-up, bead **chem-czw**.
- `specs/SPEC_nbn_low_spin.md` §6 — caveat markers added to gates **G4** and **G5** ("every
  number this gate pins comes from orbitals anchored to the vendored, metastable (10,4)/S=3 UHF
  solution").
- `specs/SPEC_nbn_low_spin.md` §7 — "Out of scope" updated with the decision + chem-czw reference
  and its ~50-70 min wall-clock cost estimate.
- `specs/SPEC_nbn_dmrg_reference.md` §8 — new **R4** caveat covering gates **G1, G2, G3, G3b**
  (all four load CAS integrals from the same vendored, metastable chk), with the same decision and
  chem-czw pointer.
- `specs/BACKLOG.md` — new `[x]` Done entry (top of `## Done`) documenting the hypothesis, the
  kill of the naive `os.environ['OMP_NUM_THREADS']` fix, the working fix, and the explicit
  "what this does NOT do" caveat.
- Follow-up bead **chem-czw** (P3, open) created for the declined full re-derivation (option a),
  so it isn't lost as an informal note.

---

## 2. Gate-run evidence (verbatim)

All three NbN-related spec gates, each run in its own process (per the block2/pyscf+qiskit-aer
OpenMP isolation rule in CLAUDE.md):

```
$ UV_PROJECT_ENVIRONMENT=/tmp/venv uv run pytest -q tests/test_nbn_scf_determinism_spec.py
..                                                                       [100%]
2 passed, 50 warnings in 22.86s

$ UV_PROJECT_ENVIRONMENT=/tmp/venv uv run pytest -q tests/test_nbn_dmrg_reference_spec.py
....                                                                     [100%]
4 passed in 341.48s (0:05:41)

$ UV_PROJECT_ENVIRONMENT=/tmp/venv uv run pytest -q tests/test_nbn_low_spin_spec.py
..                                                                       [100%]
2 passed in 484.05s (0:08:04)
```

`test_nbn_dmrg_reference_spec.py` and `test_nbn_low_spin_spec.py` load the vendored chk directly
via `nbn_dmrg_reference.load_nbn_cas` — they never call `benchmark_nbn._tight_scf`, so the fix in
§1b cannot regress them; both still pass green with the vendored (metastable) chk untouched, as
expected under decision (b). These two files' own `M` status in `git status` predates this session
(already modified in the git-status snapshot taken before I started; not touched by chem-y6w work).

Lint on touched files:

```
$ UV_PROJECT_ENVIRONMENT=/tmp/venv uv run --extra dev ruff check benchmark_nbn.py tests/test_nbn_scf_determinism_spec.py
All checks passed!
```

---

## 3. Every number, with regenerating command and vendor location

| Number | Value | Command | Vendor location |
|---|---|---|---|
| Committed (metastable) SCF energy | `-110.0275853827266 Ha` | loaded from chk, not recomputed here | `specs/nbn_scf_reference.chk` (vendored) |
| Deterministic fresh-regeneration SCF energy (the *other*, lower minimum) | `-110.02813374576796 Ha` | `UV_PROJECT_ENVIRONMENT=/tmp/venv uv run pytest -q tests/test_nbn_scf_determinism_spec.py` (G1/G2 assert this bit-identically; chem-owr's handoff has the original single-threaded measurement) | not vendored — reproduced fresh each gate run into `tmp_path`, never written to `specs/` |
| SCF-level gap between the two minima | `0.548 mHa` | difference of the two rows above | chem-owr's handoff (`sandbox-handoffs/chem-owr.md`) |
| CASCI(14,14)/FCI(10,4)-level gap between the two minima | `4.84 mHa` | chem-owr's handoff; not recomputed this session (would require the declined full re-derivation, chem-czw) | chem-owr's handoff |
| Pinned S1-S3 gap (unaffected by this bead) | `1.538 mHa` | `specs/SPEC_nbn_low_spin.md` (chem-g1i) | `specs/SPEC_nbn_low_spin.md` |

No new physics numbers are claimed by this bead's fix itself — G1/G2 are a reproducibility claim
about code (same input → same output), which per `specs/SPEC_nbn_scf_determinism.md` §3 does not
need an FCI/DMRG reference (the "reference" is reproducibility itself, checked directly by
comparing repeated runs bit-for-bit).

---

## 4. Pre-registered criteria and pass/fail

From the bead's acceptance criteria, path (b):

1. **"An explicit, spec-recorded decision to keep the current metastable chk as the reference"**
   — **PASS.** Recorded in `specs/SPEC_nbn_low_spin.md` §0 and `specs/SPEC_nbn_dmrg_reference.md`
   §8 R4 (not just this handoff).
2. **"...with a caveat on every table/gate that depends on it"** — **PASS.**
   `SPEC_nbn_low_spin.md` G4/G5 caveated; `SPEC_nbn_dmrg_reference.md` G1/G2/G3/G3b caveated
   (all four share the same vendored orbitals, so one R4 note covers all of them explicitly by
   gate name).
3. **"...plus a determinism-forcing fix ... so future regenerations are reproducible rather than
   a coin flip"** — **PASS.** `pyscf.lib.num_threads(1)` pin in `_tight_scf`; G1 (cross-thread)
   and G2 (repeat-call) both green, RED confirmed without the fix.
4. **Test suite passes** — **PASS** for every gate exercised (§2 above). Full `make gates`/
   `make test` were NOT run in full this session (would re-run all `test_*_spec.py` gates,
   including unrelated ones, each in its own process — costly; scoped instead to the three
   NbN-relevant files plus lint, since nothing else was touched). This is a partial-suite claim,
   noted explicitly rather than implied as "full suite green."

---

## 5. What was decided not to do, and why

**Declined:** re-deriving the CAS(14,14) reference numbers (all spin sectors S=0,1,2,3) on the
deterministic lower-energy UHF minimum, i.e. acceptance-criteria path (a).

**Why:** `nbn_dmrg_reference.load_nbn_cas` loads the UHF `mo_coeff`/`mo_occ` from the chk **once**
and builds every spin sector's CASCI integrals from those same orbitals (`get_h1eff`/`get_h2eff`).
Switching the reference orbitals is therefore not a one-number patch — it means re-running FCI/DMRG
for all four sectors and re-deriving `SPEC_nbn_dmrg_reference.md`'s G1-G3b and
`SPEC_nbn_low_spin.md`'s Tables 1-2 / S1-S3 gap / G4-G5, estimated at similar wall-clock cost to
chem-bbi + chem-g1i combined (~50-70 minutes of FCI/DMRG runs), which CLAUDE.md's own
"Long-running jobs" section says should be deferred to the user's own terminal rather than run in
this container. Filed as **chem-czw** (P3, open) instead of attempted here, with the exact
re-derivation scope and cost estimate recorded in its description so a future session (or the user)
can pick it up without re-diagnosing.

---

## 6. What could not be verified

- **Cross-host / cross-BLAS-backend generalization of the fix.** Verified only on this container
  (x86_64, 8 cores, OpenBLAS via `scipy-openblas64`). `specs/SPEC_nbn_scf_determinism.md` §8 R1
  flags this explicitly: a pyscf build linked against a different BLAS, or a stability-solver path
  that bottlenecks on raw numpy/BLAS calls outside pyscf's own OMP control, could still be
  non-deterministic elsewhere.
- **Whether the deterministic lower-energy minimum is itself the "true" global UHF minimum** (vs.
  a third, still-lower one) — not investigated; out of scope for a threading-determinism bug fix.
- **Full `make gates`/`make test`** was not run this session (see §4 point 4) — only the three
  NbN-relevant gate files plus lint on touched files.

---

## 7. Git status / suggested commands (NOT executed — human review required)

Per this session's git policy, no commit/push was made. Files touched by this bead's work:

```
 M benchmark_nbn.py
 M specs/BACKLOG.md
 M specs/SPEC_nbn_dmrg_reference.md
 M specs/SPEC_nbn_low_spin.md
?? specs/SPEC_nbn_scf_determinism.md
?? tests/test_nbn_scf_determinism_spec.py
?? sandbox-handoffs/chem-y6w.md
```

(Other modified/untracked files in `git status` predate this session — unrelated prior work, not
touched here.)

Suggested commands for a human to review and run:

```bash
git add benchmark_nbn.py specs/BACKLOG.md specs/SPEC_nbn_dmrg_reference.md \
        specs/SPEC_nbn_low_spin.md specs/SPEC_nbn_scf_determinism.md \
        tests/test_nbn_scf_determinism_spec.py sandbox-handoffs/chem-y6w.md
git commit -m "chem-y6w: pin pyscf thread count for deterministic NbN SCF regen; keep metastable chk as reference with caveats

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
bd close chem-y6w --reason="Decision (b): kept vendored metastable chk as reference with explicit caveats on every dependent gate (SPEC_nbn_low_spin.md G4/G5, SPEC_nbn_dmrg_reference.md G1-G3b); fixed regeneration-path determinism via pyscf.lib.num_threads(1) pin in benchmark_nbn._tight_scf (SPEC_nbn_scf_determinism.md, G1/G2 green, RED confirmed pre-fix). Full re-derivation on the alternate minimum declined as out of scope, filed as chem-czw."
```

---

## 8. Verification addendum (second session, same day)

The bead was found already `IN_PROGRESS` with the above work done and this handoff file already
written, but not yet closed. Before closing, independently re-verified rather than trusting the
prior session's notes at face value:

- `git diff -- benchmark_nbn.py`: confirmed the `lib.num_threads(1)`/`finally: lib.num_threads(ambient_threads)`
  pin is present exactly as described, wrapping the existing kernel/newton/stability loop, `lib`
  already imported (`from pyscf import gto, scf, mcscf, ao2mo, lib`).
- Re-ran the fast gate fresh (not reusing the recorded output):
  ```
  $ uv run pytest -q tests/test_nbn_scf_determinism_spec.py
  2 passed, 50 warnings in 32.86s
  ```
- Re-ran lint on touched files:
  ```
  $ uv run --extra dev ruff check benchmark_nbn.py tests/test_nbn_scf_determinism_spec.py
  All checks passed!
  ```
- Grepped `specs/SPEC_nbn_low_spin.md`, `specs/SPEC_nbn_dmrg_reference.md`, `specs/BACKLOG.md` for
  `chem-y6w`: confirmed the §0 decision paragraph, the R4 caveat covering G1-G3b, and the BACKLOG
  Done entry all exist as claimed (not just asserted in this handoff).
- Confirmed follow-up bead `chem-czw` exists, open, P3, with the re-derivation scope recorded.
- Did **not** re-run `test_nbn_dmrg_reference_spec.py` (341s) or `test_nbn_low_spin_spec.py` (484s)
  this session — neither touches `benchmark_nbn.py` (they load the vendored chk directly via
  `nbn_dmrg_reference.load_nbn_cas`), so the fix cannot regress them, and their prior session's
  green run (§2 above) plus the unchanged file contents (`git diff` shows no changes to
  `nbn_dmrg_reference.py` or the chk since that run) stand as evidence rather than being redundantly
  re-run. This is an explicit reliance decision, not an oversight.

Closed `chem-y6w` after this verification. No code changes made in this second session — the
first session's work was correct and complete as documented.
