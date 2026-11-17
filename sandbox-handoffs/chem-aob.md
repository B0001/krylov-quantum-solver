# chem-aob handoff

## Bead

`SPEC_qpe_readout_laws` G2 called the `3*lambda/2^t` constant "a real, provable upper bound." It
wasn't: it was a bare fit to one measured sweep (max ratio 2.175 on H2 CAS(2,2)), never derived, and
`3 < pi` meant it wasn't even guaranteed to be a valid ceiling on an untested system/walk variant.

## What I changed

Took the bead's suggested fix (b): swapped the constant itself to a value that actually *is*
provable, rather than just softening the word "provable" in the prose.

**The derivation** (now written into the spec and the test docstring): `run_qpe`'s argmax decode
picks the dyadic phase bin nearest the true phase `theta_0`, so its phase error is at most half a
bin, `pi/2^t`, out of a full `2*pi` circle divided into `2^t` bins. `E = lambda*cos(theta)` has
`|dE/dtheta| <= lambda` everywhere (mean value theorem, `|cos'| <= 1`), so
`err(t) <= lambda * pi/2^t` unconditionally — true for any system, not fit to this one. This is the
same "exact nearest-bin" fact already recorded (but not gated here) in `SPEC_ft_identity_shift.md`
§3d (`2^t*max|cos(theta0 +- pi/2^t) - cos(theta0)| <= pi`), which is closed/landed on this branch
already — no need for chem-bu4's (separate, unmerged) branch to make this claim.

Files touched:

- `specs/SPEC_qpe_readout_laws.md` — G2 now states and derives `err(t) <= pi*lambda/2^t`, not the
  unexplained `3`. Background, "Out of scope," and R1 caveat sections updated to match: the `pi`
  ceiling is now correctly called provable/universal; the *tight* realized ratio (2.175) is still
  explicitly flagged as a measurement, not derived — that tightening to a per-system
  `pi*sin(theta_0)` envelope remains `SPEC_qpe_precision_bound` (chem-0ke)'s job, out of scope here.
  Added an explicit revision note citing chem-aob and why `3` was wrong (empirical fit, `3 < pi`).
- `tests/test_qpe_readout_laws_spec.py` — `test_G2_precision_bounded_by_pi_lambda_over_2t` (renamed
  from `..._constant_times_lambda_over_2t`) now asserts `err <= np.pi * lam / 2**t` instead of
  `3.0 * lam / 2**t`. Docstring states the derivation inline.
- `specs/BACKLOG.md` — updated the two spots that described the old `3x`/"provable" framing to
  match: the `SPEC_qpe_readout_laws` closed-entry bullet (now says "PROVABLE ceiling `pi*lambda/2^t`"
  and records the chem-aob correction), and the `chem-177`/`SPEC_ft_identity_shift` bullet's
  cross-reference, which previously attributed its own `3lambda/2^t` constant to
  "`SPEC_qpe_readout_laws` G2" — that attribution is now stale since G2 itself uses `pi`, so I
  re-pointed it at `test_ft_identity_shift_spec.py::_t_star`'s own (unchanged) local constant.

## What I deliberately did NOT change (and why)

- **`tests/test_ft_identity_shift_spec.py::_t_star`** still hardcodes `const=3.0`. This bead's own
  text scopes that file to chem-bu4 ("Out of scope for chem-bu4 itself... Out of scope for this
  bead too by the same logic — the bead names only `SPEC_qpe_readout_laws.md` G2 as its target).
  Note for whoever picks up the adjacent work: there **is** a commit doing exactly that swap,
  `e512c25` "chem-bu4: use the derived pi hard ceiling for phase-bit budgets, not the empirical 3",
  but it lives on branch `sandbox/chem-bu4`, which is **not** merged into `main` and is not an
  ancestor of this branch (`sandbox/chem-aob`, cut from `main`). `main` and this branch both still
  have `SPEC_ft_identity_shift`'s `_t_star` at `const=3.0`, and `SPEC_qpe_precision_bound.md` in
  this tree has only G1-G3 (no G4) — the chem-bu4 G4 gate for the exact nearest-bin ceiling across
  both walks does not exist here either. If chem-bu4's branch lands later, re-check this spec's new
  wording against whatever final form it settles on (I did not depend on anything from that branch
  — the `pi` derivation here stands on its own via MVT + half-bin dyadic rounding, and is also
  independently recorded in this tree's already-closed `SPEC_ft_identity_shift.md` §3d).
- Did not touch `SPEC_qpe_precision_bound.md` — its own numbers (2.175/2.285/etc.) and framing were
  already honest (explicitly "not derived analytically," "not universal") and don't claim
  "provable"; nothing there needed correcting for this bead.
- Did not bump the spec's top `Status: DRAFT` line. It's inconsistent with the fact that the spec is
  implemented and gated (9/9 tests green, has been for a while per BACKLOG's `[x]` entry), but that
  mismatch predates this bead and isn't what chem-aob asked me to fix — flagging it here rather than
  silently changing scope.

## Gate-run lines (verbatim-equivalent)

```
$ uv run pytest -q tests/test_qpe_readout_laws_spec.py tests/test_qpe_precision_bound_spec.py
................                                                         [100%]
16 passed in 23.95s
```

Split out:
- `tests/test_qpe_readout_laws_spec.py` — 9 passed (G1-G4, including the new `pi`-based G2).
- `tests/test_qpe_precision_bound_spec.py` — 7 passed (unaffected; checked for regression since it
  directly depends on/references this spec's constant).

```
$ uv run ruff check .
All checks passed!
```

Hardware: Linux `6f0eb2fe1bf4` (Docker VM), x86_64, 8 cores (`nproc`). Both pytest runs together
took 23.95s — far under the "report wall time if >1 min" threshold, so no resource notes needed.

## Pre-registered criteria (from the bead) and outcome

The bead offered two options, (a) retitle to "empirical" or (b) swap the constant to `pi`. I took
(b), since `pi` is a genuinely provable, universal ceiling (not system-specific), so G2 can keep
calling itself "provable" truthfully instead of being downgraded to "empirical." Both options were
explicitly sanctioned by the bead; (b) keeps the gate's informational content (a hard, checkable
ceiling) rather than weakening it to a label change.

- Criterion: G2's constant must be actually derivable, not a measured fit — **PASS** (MVT + half-bin
  dyadic argument, independent of system; verified numerically: 2.175 < pi on H2 CAS(2,2)).
- Criterion: the spec's prose must not call anything "provable" that isn't — **PASS**: re-read the
  whole file after editing; the only remaining "provable"/"derived" claims are G1 (monotonicity,
  unchanged, already a direct consequence of dyadic bin refinement) and the new G2.
- Criterion: don't silently lose the honest "one-system-only" caveat — **PASS**, R1 rewritten to
  separate what's now universal (the `pi` ceiling) from what's still only measured on H2 CAS(2,2)
  (the realized 2.175 ratio).

## What I decided not to do, and why

- Did not pull in or depend on chem-bu4's unmerged `SPEC_qpe_precision_bound` G4 / pi-ceiling gate,
  since that branch isn't merged here and this bead's fix doesn't need it (see above).
- Did not change `test_ft_identity_shift_spec.py` — explicitly out of scope per the bead text, and
  chem-bu4 already addresses it on its own branch.
- Did not file a new bead for the stale `Status: DRAFT` header or for reconciling branch
  `sandbox/chem-bu4` into `main` — the latter is presumably the host's/a human's merge-queue
  decision, not a code-review finding; happy to file either if wanted, but didn't want to invent
  scope unasked.

## What I could not verify

- Whether `sandbox/chem-bu4`'s PR (if any) is still pending review, abandoned, or something else —
  I only confirmed via `git log --all` / `git branch --contains` that it exists as a local branch
  object and is not an ancestor of `main` or of this branch. Did not check GitHub for an open PR
  (no `gh` lookup performed — out of scope for this bead, flagging as unverified rather than
  guessing).

## Suggested commands for the human

```
git add specs/SPEC_qpe_readout_laws.md specs/BACKLOG.md tests/test_qpe_readout_laws_spec.py sandbox-handoffs/chem-aob.md
git commit -m "chem-aob: SPEC_qpe_readout_laws G2 now proves the pi*lambda/2^t ceiling instead of asserting the empirical 3x as 'provable'"
```

(Left uncommitted per the sandbox git policy — pre-commit hook refuses commits with
`BEADS_ACTOR=sandbox` on this branch regardless.)
