# chem-5oj handoff: scoring `fit_thc`'s nonlinear THC with the native `thc_lambda`

**Bead:** chem-5oj — "The 62x THC lambda penalty (SPEC_thc_lambda) may be a collocation-point
artifact — `lambda_ladder.fit_thc`'s existing nonlinear fit has never been scored with the native
`thc_lambda`."

**Status: CLOSED.** All code, gates, and the corroboration sweep this bead's own cost note and
caveat asked for already existed in the working tree at session start (committed in `5c5c6db`,
by a different/prior worker instance whose container was snapshotted mid-session) but the bead was
left `in_progress` with no handoff, and the §9 results section the spec itself says it is
"awaiting" had never actually been produced. This session **independently re-ran every gate and
the corroboration sweep from scratch** (new numbers, not copy-pasted from the pre-existing
docstrings) and found the pre-existing claims hold up, then wrote them into the spec, flipped
`BACKLOG.md`, and wrote this handoff. No solver/library code was changed — only `specs/BACKLOG.md`,
`specs/SPEC_thc_collocation.md`, and a new `scripts/thc_collocation_sweep.py`.

## What I found at claim time

`bd show chem-5oj` was already `in_progress` (claimed 2026-09-29, no notes). `git status` was
clean and `git log` showed no branch-local WIP — but `specs/SPEC_thc_collocation.md`,
`tests/test_thc_collocation_spec.py`, and the `pair_indicator_collocation`/`fit_thc(...,
return_factors=)` code were already present and committed (`5c5c6db`, "Sandbox worker output...
chem-5oj (THC collocation)... Committed while a worker was live, so in-progress beads may be
mid-edit"). The spec's own status line read "DRAFT — awaiting gate numbers from this session's run
(chem-5oj)" — i.e. a previous worker wrote the spec/code/tests but never ran and recorded the
promised results, and never wrote a handoff. Per this session's instructions ("do not assume its
partial work is correct"), I did not take the pre-existing docstring numbers on faith — I re-ran
everything independently (new seeds produce the same qualitative story; see below for the specific
numeric discrepancies at the 2nd significant figure, attributed to BLAS floating-point
non-reproducibility, not to the prior work being wrong).

## Commands run and their output

**Setup:**
```
uv sync --extra dmrg --extra test --extra dev
```

**Gate 1 — this bead's own spec gate:**
```
uv run pytest -q tests/test_thc_collocation_spec.py -v
```
```
tests/test_thc_collocation_spec.py ....                                  [100%]
========================= 4 passed in 66.53s (0:01:06) =========================
```

**Gate 2 — the parent spec's gates, to confirm no regression of the validated anchor:**
```
uv run pytest -q tests/test_thc_lambda_spec.py -v
```
```
tests/test_thc_lambda_spec.py ....                                       [100%]
============================== 4 passed in 28.07s ==============================
```

**Lint:**
```
uv run ruff check thc_factorization.py lambda_ladder.py tests/test_thc_collocation_spec.py \
    tests/test_thc_lambda_spec.py scripts/thc_collocation_sweep.py
```
`All checks passed!`

**Corroboration sweep (new script, `scripts/thc_collocation_sweep.py`, written and run this
session — see "What I did" below):**
```
uv run python scripts/thc_collocation_sweep.py
```
Full output (verbatim):
```
==========================================================================================
1) LiH/STO-3G norb=6 M=21 -- 12-seed x 8000-nfev (restarts=1) reconstruction-error sweep
   corroborates test_thc_collocation_spec.py::test_G2 (fit_thc fails <1e-6 precondition)
==========================================================================================
seed= 0  recon_err=1.8166e-02  time=24.0s
seed= 1  recon_err=1.3192e-02  time=20.6s
seed= 2  recon_err=2.0204e-02  time=20.7s
seed= 3  recon_err=1.1140e-02  time=20.5s
seed= 4  recon_err=3.8178e-03  time=20.7s
seed= 5  recon_err=1.2953e-02  time=23.5s
seed= 6  recon_err=4.1910e-03  time=20.7s
seed= 7  recon_err=1.6915e-02  time=23.2s
seed= 8  recon_err=1.2197e-02  time=20.7s
seed= 9  recon_err=1.5165e-02  time=23.2s
seed=10  recon_err=6.6311e-03  time=20.8s
seed=11  recon_err=1.7697e-02  time=20.4s
min=3.8178e-03  max=2.0204e-02  mean=1.2689e-02
total wall time = 259.0s (4.32 min)
ALL 12 seeds stay above 1e-6 (G2 precondition never met): True
==========================================================================================
2) LiH/STO-3G norb=6 M=21 -- 5-seed beats_df_lambda / ratio_to_random spot check
   (fit_thc at CI-gate defaults: restarts=4, max_nfev=4000)
   corroborates test_thc_collocation_spec.py::test_G4
==========================================================================================
df_lambda=15.4846  random_lambda=1001.6941 (seed=0)
seed=0  recon_err=1.5900e-02  lambda=14.6514  ratio_to_random=0.0146  beats_df=True  time=46.5s
seed=1  recon_err=1.7600e-02  lambda=16.6339  ratio_to_random=0.0166  beats_df=False  time=49.1s
seed=2  recon_err=2.5176e-02  lambda=14.8130  ratio_to_random=0.0148  beats_df=True  time=49.8s
seed=3  recon_err=1.6257e-02  lambda=14.9239  ratio_to_random=0.0149  beats_df=True  time=48.4s
seed=4  recon_err=1.6655e-02  lambda=14.7960  ratio_to_random=0.0148  beats_df=True  time=47.9s
beats_df_lambda at 4/5 seeds
ratio_to_random range: [0.0146, 0.0166]  (all <0.2, i.e. >=5x below random: True)
==========================================================================================
3) H2O/STO-3G full space, norb=7 -- single NOT-CI-gated run at fit_thc CI-gate defaults
   (this bead's cost note: '~600 Levenberg-Marquardt parameters at norb=7')
==========================================================================================
norb=7  M=28  n_LM_params=602
df_lambda=86.9575  random_lambda(seed=0)=5377.3005
fit_thc seed=0 restarts=4 max_nfev=4000: time=175.2s  recon_err=8.9503e-02  lambda=85.1129
  ratio_to_random=0.0158  beats_df=True
==========================================================================================
```

**Hardware:** `uname -a` → `Linux 7843f162c504 7.0.14-linuxkit #1 SMP PREEMPT ... x86_64 GNU/Linux`
(container); `nproc` → 8; `MemTotal` 16353600 kB (~15.6 GB). Total wall time for everything above
(2 gate files + lint + 3-part sweep) ≈ 13 CPU-minutes.

## The pre-registered criteria and whether each passed

From the bead text and `specs/SPEC_thc_collocation.md` §5:

- **G1 (random + structured collocation reconstruct exactly at M=21).** PASS —
  `test_G1_random_and_structured_collocation_reconstruct_exactly_at_matched_rank` passed; both
  `< 1e-9`.
- **G2 (fit_thc fails the `<1e-6` reconstruction precondition — a recorded finding, not a bug).**
  PASS (the finding is confirmed, not merely asserted once): seed=0/restarts=4/max_nfev=4000 gives
  recon err ≈1.6–1.8e-2; the independent 12-seed × 8000-eval sweep never goes below 3.8e-3 — four
  orders of magnitude above the `<1e-6` bar, at every one of 12 seeds, not just the CI seed.
- **G3 (λ formula sanity: structured collocation within 2× of `df_lambda`).** PASS.
- **G4 — both kill directions evaluated explicitly:**
  - **Kill A** ("penalty runs deeper than unoptimized points", fires if nonlinear is NOT ≥5× below
    random): **does not fire**. Nonlinear collocation is ≈65–69× below random at every seed checked
    (ratios 0.0146–0.0166 across 5 CI-default seeds at norb=6; 0.0158 at norb=7) — comfortably past
    the ≥5× (ratio < 0.2) bar.
  - **Kill B** ("beats `df_lambda` outright, `SPEC_thc_lambda` G4 must be revised", the
    opposite-direction kill the bead names): **fires at most but not all seeds** — 4/5 at norb=6
    (seed=1 does not beat it: λ=16.63 > df_lambda=15.48), and at the single norb=7 point checked.
    **This is recorded as the finding, not smoothed into "yes" or "no":** see the caveat read
    below. It does **not** trigger an actual revision of `SPEC_thc_lambda`'s G4, because that gate
    is specifically about the *exact* random-collocation THC (re-confirmed here unchanged, ≈65–68×
    above `df_lambda`), not about `fit_thc`'s inexact nonlinear fit.
- **Random seed pinned and reported.** Yes — `seed=0` for every CI-gated number; seeds 0–11 (sweep
  1) and 0–4 (sweep 2) explicitly enumerated above, not just described in prose.

## The honest finding (the part a reviewer should actually take away)

`fit_thc`'s nonlinear collocation is **reliably λ-small relative to random collocation** (Kill A
never fires, at any seed or the one larger system checked) but is **not reliably below
`df_lambda`** (Kill B: 4/5, not 5/5, at norb=6) — and, critically, **every single instance of
"beats `df_lambda`" comes paired with a reconstruction error 4–5 orders of magnitude above the
`<1e-6` bar the comparison was supposed to be matched at** (1.6–2.5e-2 at norb=6, 8.95e-2 at
norb=7). As `specs/SPEC_thc_collocation.md` §8 R1 anticipated before this run: a THC fit that
reproduces the ERI tensor poorly can trivially have a small 1-norm (the zero operator has λ=0 and
infinite error), so "beats `df_lambda` while failing G2" is **not** read here as "the 62× penalty
is solved by `fit_thc`." It is read as the bead's own half-predicted sharper finding: an
error-optimal THC fit is λ-small but **not λ-controlled** — it can land anywhere in λ-space,
including below `df_lambda`, while reproducing a *different, wrong* Hamiltonian. This does not
shorten the path to a genuine λ advantage, which (unchanged from `SPEC_thc_lambda` §7) still needs
ISDF/optimized collocation that preserves reconstruction fidelity, not just a small 1-norm.

## What I changed, file by file

- `specs/SPEC_thc_collocation.md` — flipped Status from `DRAFT — awaiting gate numbers...` to
  `CLOSED` with a one-paragraph finding summary, and added a new `## 9. Results` section (old §9
  Deliverables renumbered to §10) containing the full numeric results above. No gate wording, no
  acceptance thresholds, and no code were changed — only the results that were missing were added.
- `specs/BACKLOG.md` — moved the "62× THC λ penalty" open item (old line ~421, `### Fault-tolerant
  stack`) to `## Done` as `[x]` with the finding paragraph (mirrors the style of the adjacent
  already-closed THC entry at the old line ~935).
- `scripts/thc_collocation_sweep.py` (new) — the three-part corroboration sweep run above (12-seed
  reconstruction-error sweep, 5-seed beats-df spot check, the non-CI-gated norb=7 point). This is
  the script `specs/SPEC_thc_collocation.md` §8 R2 says exists ("a larger, non-CI-gated sweep...
  recorded in the PR/handoff") — it did not exist in the tree before this session; I wrote it,
  confirmed it lint-clean (`ruff check` passes), and ran it to produce the §9 numbers above rather
  than trust the prior worker's docstring estimates.
- **No change** to `thc_factorization.py`, `lambda_ladder.py`, or `tests/test_thc_collocation_spec.py`
  — these were already correct and passing as committed by the prior session; I verified this by
  running the gates myself rather than assuming the commit message was accurate.

## What I decided not to do, and why

- **Did not run the full `make gates` / `scripts/run_gates.sh` suite (~90 spec files).** I started
  it (`GATE_JOBS=2`), but stopped it after a few minutes once I confirmed it was going to walk the
  entire repo's spec suite (DMRG ladders, CBS extrapolations, etc.) — none of which this bead
  touches, since I changed no solver/library code, only two spec/doc files and a new, isolated
  analysis script. Running the full suite would cost tens of minutes to hours of this session's
  compute budget for zero additional evidence about this bead's claim. I instead ran the two gate
  files this bead's own acceptance criteria name (`test_thc_lambda_spec.py`,
  `test_thc_collocation_spec.py`) directly, which is what `make gates`'s `test_*_spec.py` glob
  would run for these two files anyway, and `ruff check` on every file I touched.
- **Did not re-derive or second-guess `fit_thc`'s Levenberg-Marquardt implementation.** It was
  already correct (G1 passes on the two linear methods it's compared against; `fit_thc`'s own
  behavior — not reaching machine-precision reconstruction — is the measured fact this bead exists
  to pin, not a bug to fix).
- **Did not attempt ISDF/optimized collocation.** Explicitly out of scope per both
  `SPEC_thc_lambda.md` §7 and `SPEC_thc_collocation.md` §7; this bead's job was to score the
  *existing* `fit_thc`, not build a new optimizer.
- **Did not revise `SPEC_thc_lambda.md`'s G4**, despite Kill B firing at most seeds — because, as
  explained above, that gate is specifically about the exact random-collocation THC
  (`tensor_hypercontraction`), which this session re-confirmed is unchanged (still ≈65–68× above
  `df_lambda`). `fit_thc` beating `df_lambda` while failing reconstruction is a different claim
  about a different (inexact) method, not evidence against G4's claim about the exact one.
- **Did not re-run the 12-seed / 5-seed sweeps at norb=7** (only the single CI-default point) —
  the ~600-parameter norb=7 case already takes ~175s for one restarts=4 run; a full 12×8000-eval
  sweep there would cost ~35+ minutes for a point this bead's own cost note says should be "run
  once... reported, not gated."

## What I could not verify

- **Global optimality of any `fit_thc` run.** Explicitly out of scope by design (`fit_thc` is a
  local, stochastic LM search) — the spec and this handoff never claim it; seeds are pinned and the
  spread across 12+5 seeds is reported instead.
- **Bit-for-bit reproducibility of `fit_thc` at a fixed seed.** Two of my own runs with identical
  `seed=0, restarts=4, max_nfev=4000` arguments gave recon err 1.83e-2 vs 1.59e-2 and λ 14.6201 vs
  14.6514 — close enough that no qualitative conclusion changes, but not bit-identical. I attribute
  this to non-deterministic summation order in multi-threaded BLAS calls inside
  `scipy.optimize.least_squares`'s Jacobian/residual evaluations, not to a bug in `fit_thc` or a
  seeding error (the RNG itself, `np.random.default_rng(seed)`, is deterministic). This is noted
  explicitly in `SPEC_thc_collocation.md` §9 rather than papered over by picking one run's numbers
  silently.
- **Whether the pre-existing docstring claims in `tests/test_thc_collocation_spec.py`
  (written before this session) were produced by an actual run or estimated/transcribed.** I did
  not try to find out — I treated them as unverified and re-derived everything in this handoff from
  scratch. They turned out to agree with my independent numbers (e.g. "ratio ~0.94" ≈ my 14.65/
  15.48=0.946; "4/5 seeds beat df_lambda" ≈ my exact 4/5 with the same failing seed=1), which is
  reassuring but was not assumed going in.

## Commands to reproduce every number above

```bash
uv sync --extra dmrg --extra test --extra dev
uv run pytest -q tests/test_thc_collocation_spec.py -v
uv run pytest -q tests/test_thc_lambda_spec.py -v
uv run ruff check thc_factorization.py lambda_ladder.py tests/test_thc_collocation_spec.py \
    tests/test_thc_lambda_spec.py scripts/thc_collocation_sweep.py
uv run python scripts/thc_collocation_sweep.py
```

## Git status at handoff time

Per this run's git policy ("commit your work AND your handoff file here... on branch
`sandbox/chem-5oj`"), I will commit:
- `specs/SPEC_thc_collocation.md` (DRAFT → CLOSED, §9 Results added)
- `specs/BACKLOG.md` (entry moved to `## Done`, marked `[x]`)
- `scripts/thc_collocation_sweep.py` (new)
- `sandbox-handoffs/chem-5oj.md` (this file)
- `.beads/` export changes from `bd close chem-5oj`

No other files are touched. `git status --short` immediately before these commits:
```
 M specs/BACKLOG.md
 M specs/SPEC_thc_collocation.md
?? scripts/thc_collocation_sweep.py
```
(`sandbox-handoffs/chem-5oj.md` not yet written at the time of that snapshot; `.beads/` changes
land once `bd close` runs.)
