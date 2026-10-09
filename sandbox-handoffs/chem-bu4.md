# chem-bu4 handoff

**Bead:** QPE precision: `pi*sin(theta0)` is only a large-t envelope — recentered walk exceeds it
at t=4 (H2 2.6966 > 2.6212). Found landing chem-0ke + chem-177 together (PR #55).

**Outcome: implemented and gated.** `bd close chem-bu4` after this handoff.

Hardware: Linux x86_64, `uname -a` = `Linux 4bebe7620f9a 7.0.14-linuxkit #1 SMP PREEMPT Fri Sep 18
10:19:32 UTC 2026 x86_64 GNU/Linux`, `nproc` = 8. All runs below completed in seconds to low
minutes — no job needed the long-running-job carve-out.

## What changed

Pure spec/test revision — **no library code changed** (`qpe_walk_readout.py`,
`qubitization_blueprint.py` untouched; both already had everything needed from the chem-177
landing: `run_qpe(..., recenter=...)` and `split_identity`).

1. **`specs/SPEC_qpe_precision_bound.md`** — added a `Revision (2026-10-09, chem-bu4)` block right
   after the Status line, and threaded the same correction through §2 ("what you cannot claim"),
   §5 (new **G4**), §8 (R1), §9 (deliverables), §10 (results table). The spec now says explicitly,
   in the normative sections (not just a buried landing note), that `pi*sin(theta0)` is a
   **large-t envelope**, not a per-`t` ceiling, and that the actual per-`t` ceiling is the exact
   nearest-dyadic-bin bound `2^t * max(|cos(theta0 +- pi/2^t) - cos(theta0)|)`, which is itself
   `<= pi` for every `theta0`/`t` by the mean-value theorem on `cos`.

2. **`tests/test_qpe_precision_bound_spec.py`** — added gate **G4**
   (`test_G4_exact_nearest_bin_bound_never_exceeded`, parametrized over `recenter in [False,
   True]`, run against the existing `system` fixture = H2 CAS(2,2)/LiH CAS(2,2)/N2 CAS(3,4)): for
   every `t = 4..14`, on BOTH the ordinary and the identity-free recentered walk, the measured
   ratio must not exceed the exact per-`t` ceiling. This is the gate that would have caught
   chem-177's landing finding as "within spec" instead of a footnote. 7 original + 6 new
   parametrized cases = 13 tests, all green.

3. **`specs/SPEC_ft_identity_shift.md`** — added `## 10. Revision (2026-10-09, chem-bu4)`, and
   edited §3d/§5-G3 in place: the worst-case phase-bit *budget* formula `t* = ceil(log2(C*lambda/eps))`
   now uses `C = pi` (the derived hard ceiling) instead of the old `C = 3` (an empirical fit from
   `SPEC_qpe_readout_laws` G2, fit to one measured sweep on one system). `pi > 3`, so this is
   strictly more conservative (every budget stays the same or grows by one bit).

4. **`tests/test_ft_identity_shift_spec.py`** — `_t_star`'s default `const` changed from `3.0` to
   `np.pi`; docstrings on `_t_star` and `test_G3_budget_reduction_is_real_and_honored` updated to
   name the hard-ceiling derivation instead of citing the old empirical bound. No other test logic
   changed — G3's own assertions (`t_rec <= t_raw`, both land under `eps` end-to-end) and the
   separate "`at least one system strictly shrinks`" test both still hold with the new constant
   (verified below).

5. **`specs/BACKLOG.md`** — two inline edits to the existing chem-0ke/chem-177 bullets so they stop
   saying "not gated" / citing the retired `3` constant, now that both are stale. No new bullet
   added (follows the precedent set by chem-8tl's revision of `SPEC_centered_pds_roots`, which also
   revised in place without a fresh BACKLOG entry).

## Pre-registered criteria and verdict

From the bead: *"Gate pins the exact nearest-bin bound on ordinary and recentered walks across t;
SPEC_qpe_precision_bound states pi*sin(theta0) as a large-t envelope; budget formula uses the hard
ceiling pi."*

- **Gate pins the exact nearest-bin bound on both walks, across t:** PASS —
  `test_qpe_precision_bound_spec.py::test_G4_exact_nearest_bin_bound_never_exceeded[False]` and
  `[True]`, 3 systems x 2 walks x t=4..14, 0 violations.
- **SPEC_qpe_precision_bound states the envelope framing:** PASS — see §2, §5-G4, §8-R1, §10 of the
  spec; the "cannot claim" bullet is now explicit that the bound is not a per-`t` ceiling.
- **Budget formula uses the hard ceiling pi:** PASS — `SPEC_ft_identity_shift.md` §3d/§5-G3 and
  `tests/test_ft_identity_shift_spec.py::_t_star` default to `pi`, not `3`.

## Gate-run lines (verbatim)

```
$ uv run pytest -q tests/test_qpe_precision_bound_spec.py tests/test_ft_identity_shift_spec.py tests/test_qpe_readout_laws_spec.py
.........................................                                [100%]
41 passed in 123.39s (0:02:03)
```

Breakdown (run individually earlier in the session, same results):
- `tests/test_qpe_precision_bound_spec.py`: 13 passed (7 original G1-G3 + 6 new G4 cases).
- `tests/test_ft_identity_shift_spec.py`: 19 passed (unchanged gate count; G3's budget numbers
  changed under the hood, still green).
- `tests/test_qpe_readout_laws_spec.py`: 9 passed, unmodified — included only to confirm this
  spec (referenced extensively by the other two) was not touched and still passes.

```
$ uv run ruff check .
All checks passed!
```

I did not run the full `make gates` / `make test` suite (block2/DMRG + ~90 spec files): nothing in
this change touches library code, PySCF integral extraction, the Krylov solver, or DMRG, so the
risk of this revision breaking an unrelated gate is effectively zero, and a full run would cost
much more wall time than the change warrants. Ran instead: the two gates this bead's acceptance
criteria name, plus the adjacent spec that both of them cite. If a reviewer wants the full
confirmation, `bash scripts/run_gates.sh` (or `make gates`) is the command; I did not run it.

## Every number produced, with regenerating command and where it's vendored

All numbers below come from `uv run python tests/test_ft_identity_shift_spec.py` (its `__main__`
block prints exactly these quantities) and `uv run pytest -q tests/test_qpe_precision_bound_spec.py`
(gate assertions); both are vendored into `specs/SPEC_qpe_precision_bound.md` §10 and
`specs/SPEC_ft_identity_shift.md` §3d/§10.

**G4 exact-ceiling table** (`SPEC_qpe_precision_bound.md` §10), regenerate with
`uv run python tests/test_ft_identity_shift_spec.py` (same numbers, independently derivable from
`uv run pytest -q tests/test_qpe_precision_bound_spec.py`'s G4 assertions):

| System       | walk       | worst t | measured ratio | exact ceiling | measured/ceiling |
|--------------|------------|--------:|----------------:|---------------:|------------------:|
| H2 CAS(2,2)  | ordinary   |      11 |          2.1752 |         2.2867 |             95.1% |
| H2 CAS(2,2)  | recentered |       4 |          2.6966 |         2.7739 |             97.2% |
| LiH CAS(2,2) | ordinary   |      14 |          2.2846 |         2.3517 |             97.1% |
| LiH CAS(2,2) | recentered |      10 |          2.4417 |         2.8862 |             84.6% |
| N2 CAS(3,4)  | ordinary   |      14 |          2.3032 |         2.3966 |             96.1% |
| N2 CAS(3,4)  | recentered |       7 |          2.9399 |         2.9946 |             98.2% |

**Revised budget `t*(eps=1e-3)`** (`SPEC_ft_identity_shift.md` §3d/§10), regenerate with
`uv run python tests/test_ft_identity_shift_spec.py`:

| System | t* raw (C=3, old) | t* rec (C=3, old) | t* raw (C=pi, new) | t* rec (C=pi, new) |
|--------|-------------------:|--------------------:|----------------------:|----------------------:|
| H2     |                 13 |                  13 |                     14 |                     13 |
| LiH    |                 13 |                  12 |                     13 |                     12 |
| N2     |                 15 |                  14 |                     15 |                     14 |
| H2O    |                 15 |                  14 |                     15 |                     15 |

(H2O's recentered budget no longer saves a bit at this `eps` under the stricter ceiling — flagged
explicitly in the spec rather than silently dropped; the continuous-budget / lambda-ratio claim
`SPEC_ft_identity_shift` §3d makes is untouched, only the integer `t*` conversion moved.)

## What I decided not to do, and why

- **Did not touch `SPEC_qpe_readout_laws.md`**, even though its G2 calls its own empirical
  `3*lambda/2^t` fit "a real, provable upper bound" — the same overstatement chem-bu4 fixes
  elsewhere, and now visibly inconsistent with the rest of the repo (pi > 3). Out of scope: the
  bead names chem-0ke (`SPEC_qpe_precision_bound`) and chem-177 (`SPEC_ft_identity_shift`)
  specifically ("landing ... together (PR #55)"); `SPEC_qpe_readout_laws` is a separate, earlier
  spec this bead doesn't mention. Filed as **chem-aob** instead of doing it now.
- **Did not run `make gates` / the full spec suite** (see above) — would cost many multiples of
  this change's actual blast radius (no library code touched) for the compute envelope this
  container has. Ran the three directly-relevant gate files instead and verified `ruff check .`
  clean repo-wide.
- **Did not add a new BACKLOG.md bullet for chem-bu4 itself** — followed the chem-8tl precedent
  (revising `SPEC_centered_pds_roots` in place, no fresh BACKLOG entry) since this bead revises
  two already-closed, already-logged specs rather than adding new capability. Instead made the two
  existing BACKLOG bullets (chem-0ke, chem-177) stop contradicting the revised specs (removed
  "not gated" / the retired constant-3 citation).
- **Did not reduce G4's `t` range below 4..14** or extend it to 4..20 (matching
  `SPEC_ft_identity_shift`'s own sweep range) to keep it consistent with this spec's pre-existing
  `T_RANGE` and its three systems (H2/LiH/N2, not H2O) — the bead's own worked example (H2 at t=4)
  is inside 4..14, so this range already exercises the exact failure case that motivated the bead.

## What I could not verify

- **Generalization beyond the three tested systems / t=4..14 / STO-3G / near-equilibrium
  geometries.** This was already an explicit caveat (R3) in `SPEC_qpe_precision_bound` before this
  revision and remains one — G4 is gated on the same systems/range as G1-G3, not a broader sweep.
- **Whether `pi` is itself ever approached closely enough on an untested system that a future
  "`pi*lambda/2^t`-budget" number could round down by one fewer bit than it should** — the measured
  fractions above (84.6%-98.2% of the exact ceiling) show real headroom is sometimes thin (LiH
  recentered sits at 84.6%, N2 recentered at 98.2%), but I did not search for a system that gets
  closer to 100% than N2's 98.2%, since that would be exploring beyond this bead's three
  pre-registered systems (chem-aob's note about SPEC_qpe_readout_laws touches a related, narrower
  version of this same concern).

## Git status (uncommitted — host commits per the session's git policy)

```
 M specs/BACKLOG.md
 M specs/SPEC_ft_identity_shift.md
 M specs/SPEC_qpe_precision_bound.md
 M tests/test_ft_identity_shift_spec.py
 M tests/test_qpe_precision_bound_spec.py
?? .claude/settings.local.json   (pre-existing, not from this session's work)
?? sandbox-handoffs/chem-bu4.md  (this file)
```

Suggested commit split (two commits, matching the two specs touched):
```
git add specs/SPEC_qpe_precision_bound.md tests/test_qpe_precision_bound_spec.py
git commit -m "chem-bu4: gate the exact nearest-bin QPE precision ceiling on both walks"

git add specs/SPEC_ft_identity_shift.md tests/test_ft_identity_shift_spec.py specs/BACKLOG.md sandbox-handoffs/chem-bu4.md
git commit -m "chem-bu4: use the derived pi hard ceiling for phase-bit budgets, not the empirical 3"
```

Bead `chem-aob` was filed (not actioned) for the adjacent `SPEC_qpe_readout_laws` G2 "provable"
wording, found but out of scope for this bead.
