# Handoff: chem-loe — classical shadows for ⟨H⟩/⟨H²⟩ vs the λ² model, Temple coverage

## Starting state

This bead was already `in_progress` (claimed 2026-09-29, no notes) when I picked it up. A previous
worker had committed a near-complete implementation in `5c5c6db` ("Sandbox worker output... chem-loe
(shadow_temple.py)... in-progress"): `shadow_temple.py`, `tests/test_shadow_temple_spec.py`, and
`specs/SPEC_shadow_temple.md` (status DRAFT) all existed, composing the pre-existing closed specs
`SPEC_classical_shadows` and `SPEC_certified_noise`. I did **not** assume that work was correct —
I ran the gate first, cold, before touching anything.

## What I found: the gate did not pass as left

```
uv run pytest -q tests/test_shadow_temple_spec.py -v
```
→ **2 failed, 3 passed** (G1, G4, G4-sanity passed; G2 and G3 failed). Both failures were real,
substantive, not flaky — I reproduced each at 10 independent seeds before touching the test.

### G2: a real units bug, not a scientific finding

`test_G2` compared `shadow_norm` (a variance-scale HKP bound, same units as λ²) against the bare
amplitude ratio `λ_H2/λ_H` instead of the λ²-scale ratio `(λ_H2/λ_H)²`. Measured: λ-ratio = 7.21,
λ²-ratio = 52.0, shadow ratio = 29.7 — inside the squared bracket (29.7 < 52.0, confirming "shrinks
not removed"), outside the unsquared one (29.7 > 7.21, hence the original failure). The backlog's
own cited "λ-based ratio of 35.6" is itself a λ²-scale number: `(62.9/10.3)² ≈ 37.3` at the
scout-probe's pre-`SPEC_lambda_h2_bridge` λ values — confirming the squared reading was the one
originally intended, just implemented with the wrong units in the gate. Fixed in
`tests/test_shadow_temple_spec.py` and the matching print in `shadow_temple.py`'s `__main__`;
amended `specs/SPEC_shadow_temple.md` §5/§1 to state the ratio correctly and record why.

### G3: not a bug — a real, reproducible, unpredicted finding

`test_G3` asserted the empirical single-shot variance stays within 5% of `shadow_norm` (the HKP
bound) for both H and H². It does for H (ratio measured/bound 0.69–0.95 over 10 seeds, 16k
snapshots each) but **does not** for H²: ratio 1.16–1.49× across the same 10 seeds — consistently
above 1, never below, which rules out a one-off sampling fluctuation. I checked this directly rather
than trust the first failing seed:

```
uv run python -c "... loop over seed in range(5,15), print empirical_var/shadow_norm for H and H2 ..."
```
→ H: 0.688–0.950 (bound holds, every seed). H²: 1.160–1.486 (bound violated, every seed).

Mechanism (not fully derived, stated honestly as a hypothesis in the spec): the HKP-style additive
formula `shadow_norm = Σ|c_k|²3^{w_k}` treats each Pauli term's contribution to variance as
independent, ignoring cross-term correlations. `SPEC_classical_shadows` validated this formula only
on H (185 terms, H₂'s small system) and it held comfortably there, and continues to hold on H₄'s
larger H (185 terms). H²'s 1775-term, heavily support-overlapping expansion (O(N⁸) blowup) is the
first operator in this repo's history where that independence approximation has actually been
stress-tested at this term count — and it breaks. This does not kill the headline efficiency claim
(empirical variance for H² is still well under λ²: 1384–1774 vs λ²=3050.85, so shadows remain
cheaper than the repo's assumed noise model), but it does mean `shadow_norm` itself is not a safe
shot-budget bound for composite operators like H². I revised G3 to assert the measured direction
honestly (bound holds for H, bound violated for H² by >1.1×) rather than loosen the tolerance to
hide it — per this repo's "a gate that proves unsatisfiable is the finding" convention.

## What I changed

- `tests/test_shadow_temple_spec.py` — fixed G2's units bug (squared λ-ratio); rewrote G3 to assert
  the bound holds for H and is reproducibly violated for H² (pinning the violation itself, not
  loosening a tolerance).
- `shadow_temple.py` — fixed the matching units bug in the `__main__` summary print; added a
  bound-holds/violated annotation per operator; expanded the module docstring to state the G3
  finding.
- `specs/SPEC_shadow_temple.md` — status DRAFT → CLOSED; corrected the G2/G3 acceptance-criteria
  text (squared ratio; the H vs H² bound-violation split); updated the stale pre-`lambda_h2_bridge`
  λ numbers in §1 to the currently-measured values (λ_H=7.66, λ_H2=55.23 vs the old 62.9/10.3,
  noting the correction and that the qualitative finding is unchanged); filled in §9 with the real,
  regenerable numbers below; added an R1 follow-up note that G3's direct measurement is what caught
  the bound violation a bound-vs-bound-only gate would have missed.
- `specs/BACKLOG.md` — closed the "Method rungs" entry (was line ~322) `[x] CLOSED (chem-loe)` with
  the three-part finding (efficiency claim confirmed, coverage claim confirmed, shadow_norm-bound
  claim falsified for H²) and the units-bug note.
- No other files touched. `.beads/issues.jsonl` changed only via `bd update --claim`.

**Not touched, pre-existing in the working tree, unrelated to this bead:** `FINDINGS_chem-3z8.md`,
`SUMMARY_chem-3z8.md`, `specs/SPEC_hchain_vsdp_bracket.md`, `vsdp_hchain_lower_bound.py` — these
belong to `chem-3z8`, a different in-progress bead owned by a different assignee (`B0001`) sharing
this filesystem. Left untouched and uncommitted; the human reviewer should not attribute them to
this handoff. `.claude/settings.local.json` was also already untracked at session start — not mine,
left alone.

## Commands that regenerate every number in this report

```
uv run pytest -q tests/test_shadow_temple_spec.py -v
```
→ **5 passed** (G1, G2, G3, G4, G4-sanity), 118.6 s wall time.

```
uv run python shadow_temple.py
```
→ full printed summary (reproduced below), ~120 s wall time:
```
H4/STO-3G, M=12: H has 185 Pauli terms, H^2 has 1775 (O(N^8) blowup, recorded not hidden)
lambda_H=7.66  lambda_H2=55.23  (lambda^2 ratio 52.04)
H: shadow_norm=   40.18  lambda^2=   58.63  empirical_var=   29.75  (empirical/shadow_norm=0.740, bound HOLDS)  (empirical/lambda^2=0.507)
H2: shadow_norm= 1193.38  lambda^2= 3050.85  empirical_var= 1458.91  (empirical/shadow_norm=1.222, bound VIOLATED)  (empirical/lambda^2=0.478)
shadow_norm ratio H2/H = 29.70  vs  lambda^2 ratio = 52.04  (asymmetry shrinks, does not vanish)
------------------------------------------------------------------------------
shadow-estimated-moment coverage (200 trials, 4000 shots each, 113.9s): cov_raw=0.395  cov_upper=0.505
lambda^2-model coverage (4000 trials, matched shot budget): cov_raw=0.385  cov_upper=0.487
```

```
uv run python -m pytest -q tests/test_classical_shadows_spec.py tests/test_certified_noise_spec.py
```
→ **8 passed** — regression check on the two closed specs this one composes; confirms nothing in
their shared surface (`shadow_norm`, `hamiltonian_one_norms`, `collect_classical_shadow`) regressed.

```
uv run ruff check shadow_temple.py tests/test_shadow_temple_spec.py
```
→ `All checks passed!` (repo-wide `ruff check .` has 18 pre-existing errors, all in files this bead
never touched — `vsdp_hchain_lower_bound.py`, `scripts/spec_pm3_subspace_eta_bound.py`, etc., from
other in-progress sandbox beads).

Hardware: `Linux ... 7.0.14-linuxkit x86_64`, 8 vCPUs (`nproc`), single worker, no GPU. No run here
exceeded ~2 minutes.

## Pre-registered criteria vs outcome

| Criterion | Outcome |
|---|---|
| G1 unbiased on Krylov Ritz state (4·stderr, 16k snapshots) | **PASS** |
| G2 shadow_norm < λ² for both H, H²; shadow ratio < λ²-ratio; shadow ratio > 1 | **PASS** (after fixing the units bug — see above) |
| G3 empirical variance < λ² for both (the efficiency claim) | **PASS** |
| G3 (as literally specified) empirical variance ≤ shadow_norm for both | **FALSE for H²** — recorded as the finding, not hidden |
| G4 200-seed Temple coverage: shadow coverage not materially (>0.15) above λ²-model coverage | **PASS** — Δ = 0.395 − 0.385 = 0.01 |

The bead's headline question — "do shadows buy Temple-bracket coverage?" — is answered **no**: the
variational knife-edge (ρ₀ → E₀) is estimator-independent, exactly as hypothesized. The bead's
secondary claim — "shadows are cheaper than the λ² model" — survives for both moments. The one part
of the original hypothesis that does **not** survive as stated is that the HKP `shadow_norm` bound
itself is trustworthy for H²; it measurably is not, by 16–49% across seeds.

## What I decided not to do

- Did not re-derive or fix the HKP shadow-norm formula to account for Pauli cross-correlations —
  out of scope for this bead (it would be new estimator theory/code, not a bridge between two
  existing closed specs). Filed as a follow-up candidate below rather than attempted here.
- Did not re-run `SPEC_classical_shadows`'s own gates with H² added as a new case (would touch a
  closed spec outside this bead's scope) — the regression run above only confirms the existing H-only
  gates still pass unchanged.
- Did not investigate whether grouped/derandomized shadows would close the H² bound gap — explicitly
  out of scope per `SPEC_shadow_temple.md` §7 (inherited from `SPEC_classical_shadows`).

## What I could not verify

- The mechanism behind the H² bound violation (cross-term correlations between overlapping Pauli
  strings) is my best explanation given the measured pattern (H holds, H² doesn't, reproducibly), not
  a derivation I verified analytically. Flagged as a hypothesis, not a proven cause, in
  `tests/test_shadow_temple_spec.py` and `specs/SPEC_shadow_temple.md` §5/§8.
- Single molecule (H₄/STO-3G), single Krylov depth (M=12) — generalizing the H²-bound-violation
  finding to other systems/term-counts is unverified and explicitly out of scope (§7).

## Follow-up worth a new bead (not filed — out of this bead's scope, flagging for a human)

Whether the `shadow_norm` HKP-bound violation found here for H² is specific to Jordan-Wigner
molecular H² or a general property of random-Pauli shadows on any operator with many
overlapping-support terms would be worth its own scout probe (e.g. a synthetic high-term-count
operator, or H² on a second molecule) before it's trusted as a general statement about the shadow
protocol rather than an H₄-specific artifact.

## Suggested commands for the human

```bash
git add shadow_temple.py tests/test_shadow_temple_spec.py specs/SPEC_shadow_temple.md \
        specs/BACKLOG.md .beads/issues.jsonl sandbox-handoffs/chem-loe.md
git commit -m "shadow_temple: fix G2 units bug, record G3's real shadow_norm bound violation on H^2 (chem-loe)"
```

(`.beads/issues.jsonl` only reflects the `bd update chem-loe --claim` / `bd close chem-loe` status
changes.) Do not push — per sandbox git policy, a human pushes after reviewing the diff.
