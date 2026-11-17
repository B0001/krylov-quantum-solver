# chem-phh — Nb3X8 Tc overprediction: pre-registering (and testing) the three candidates chem-g78 left open

## Outcome

Closed `chem-phh`. Of the three candidates `SPEC_nb3x8_magnetometry.md` §7 named for the Nb₃Cl₈/
Nb₃Br₈ 5.3×/2.3× Tc overprediction (chem-g78 already killed "coordination/mean-field reduction"),
I pre-registered a killable check for each and **executed the one that's computable from data
already vendored in this repo** — "a direct lattice renormalization of t_s⊥" — and it is **also
killed**, in the narrow (same dimerized structural motif) form stated. The other two (cooperative/
first-order structural transition; in-plane kagome exchange) are pre-registered but **blocked**: they
need ab-initio/experimental data this repo doesn't vendor and this sandbox cannot fetch. Filed two
follow-up beads (`chem-cfl`, `chem-b0y`) to source that data later; did not attempt to source it now
(out of scope — "Work outside the current bead's scope gets filed as a new bead, not done now").

Also addressed the bead's "Also" clause: re-attempted to re-source a04's measured reference numbers
(1.10 eV, 0.63 eV, ~12%) via `WebSearch`. Still blocked — same canned non-answer chem-a04 recorded.
Documented, not silently dropped.

## What I changed

1. **New file** `nb3x8_tc_attribution.py` — `required_t_rescale`, `same_motif_t_spread`,
   `max_family_same_motif_t_spread`, `ht_lt_t_ratio`, `candidate3_verdict`, `UNTESTED_CANDIDATES`
   (candidates 1 & 2, each with `description`/`pre_registered_check`/`blocked_on`/`follow_up_bead`),
   and a CLI. Composes only already-validated primitives (`nb3x8_magnetometry.predicted_transition_K`,
   `nb3x8_gaps.NB3X8_CLUSTERS`/`NB3X8_LT_BULK`) — no new physics.
2. **New file** `tests/test_nb3x8_tc_attribution_spec.py` — gates G1–G5 (test-first origin; see
   `specs/SPEC_nb3x8_tc_attribution.md`).
3. **New file** `specs/SPEC_nb3x8_tc_attribution.md` — the spec (goal, approach, gates, caveats).
4. **Edited** `specs/SPEC_nb3x8_magnetometry.md` §7 — appended an "Update (chem-phh...)" note pointing
   to the new spec and recording that candidate 3 is now also killed.
5. **Edited** `nb3x8_gaps.py` — appended a short note to the module docstring (no code change) that
   the "lattice renormalization of t_s directly" candidate it names is now tested and killed, pointing
   to the new module.
6. **Edited** `specs/BACKLOG.md` — new `- [x] CLOSED 2026-10-15 (chem-phh)` entry under "Nb₃X₈
   materials line," between the chem-g78 and chem-a04 entries, with the exact numbers and the
   WebSearch re-confirmation.
7. **Filed** `chem-cfl` (P3, task) and `chem-b0y` (P3, task) — the two blocked candidates' data-
   sourcing follow-ups, each with the pre-registered check as acceptance criteria. Not worked on here.

No existing library code's **behavior** changed (only a docstring line added to `nb3x8_gaps.py`).

## Candidate 3 — the executed check, exact numbers and commands

Reproduce with:
```
uv run python nb3x8_tc_attribution.py
uv run pytest tests/test_nb3x8_tc_attribution_spec.py -v
```

**The check:** hold `U0`, `Us` fixed (local on-site/density terms) and solve for the multiplier `eta`
on the interlayer hopping `t_s⊥` such that the already-validated
`nb3x8_magnetometry.predicted_transition_K(U0, eta*t, Us)` lands exactly on the measured Tc. Compare
`|1 − eta|` to an *independent* slice of the same vendored table (`arXiv:2501.10320` Table I, already
in `nb3x8_gaps.NB3X8_CLUSTERS`): the relative `|t_s⊥|` difference between the bulk and bilayer
stacking variants of the *same* halide's LT (dimerized) structure — a real, already-tabulated measure
of how much this hopping actually moves between closely related variants of the same phase.

**Numbers** (from `uv run python nb3x8_tc_attribution.py`, reproduced in `test_G1`–`test_G4`):

```
Candidate 3 (direct t_s_perp renormalization, same dimerized motif):
  empirical anchor (max bulk-vs-bilayer |t_s_perp| spread, family-wide): 2.04%
    Nb3F8: 2.041%   Nb3Cl8: 0.147%   Nb3Br8: 0.118%   Nb3I8: 0.092%

  Nb3Cl8: needs t_s_perp cut by 57.7% to match measured Tc (eta=0.4229);
          kill bar = 20.4% (10x the 2.04% anchor) -> KILLED
  Nb3Br8: needs t_s_perp cut by 35.8% to match measured Tc (eta=0.6424);
          kill bar = 20.4% (10x the 2.04% anchor) -> KILLED

  context only: HT/LT |t_s_perp| ratio = 0.126 (Cl) / 0.121 (Br) -- a real
  PHASE change moves t_s_perp far more than a same-motif renormalization could.
```

**Verdict:** a same-structural-motif renormalization of `t_s⊥` alone cannot explain the
overprediction for either halide — the needed shrinkage (36–58%) is 17–600× the largest same-motif
spread this family's own data shows (2.04%, at the *most permissive* reading using Nb₃F₈'s larger
spread, 28–600× at Cl/Br's own ~0.12–0.15% spread). **KILLED.**

**Caveat, stated up front (also in the spec, §8 R2):** this is a *narrow* claim. It does not rule out
`t_s⊥` moving as part of a genuine structural phase shift — the HT/LT ratio (0.126/0.121) shows the
hopping *can* move by a large factor under an actual phase change. That's candidate 1's territory
(the cooperative/first-order structural transition), not contradicted by this result.

## Candidates 1 & 2 — pre-registered, not executed

Both need data not in this repo and not fetchable here:

- **Cooperative/first-order structural transition** (candidate 1): needs elastic/lattice free-energy
  data for Nb₃X₈ (phonon/interlayer-relaxation energies) — `arXiv:2501.10320`'s vendored tables cover
  electronic cRPA parameters only, no elastic constants. Pre-registered check (in
  `nb3x8_tc_attribution.UNTESTED_CANDIDATES["cooperative_structural_transition"]`): build
  `F_phase(T) = E_mag,phase(T) − T·S_mag,phase(T) + E_lattice,phase` for HT vs LT using the already-
  validated `nb3x8_thermo` Boltzmann-trace machinery plus the (missing) `E_lattice`; solve the
  crossing `Tc`; kill if no physically plausible `E_lattice` reproduces the measured Tc. Follow-up:
  `chem-cfl`.
- **In-plane kagome exchange** (candidate 2): needs an in-plane exchange constant `J'` — not
  tabulated anywhere in this repo or in `arXiv:2501.10320`'s Table I/IV (interlayer dimer only).
  Pre-registered check: compute the frustration ratio `J'/J_interlayer` and test against standard
  kagome frustration-suppression formulas for whether a physically plausible `J'` reaches the needed
  5.3×/2.3× suppression. Follow-up: `chem-b0y`.

## The a04 "never re-sourced" clause

The bead also flagged: "a04's measured reference values (1.10 eV, 0.63 eV, ~12%) were never
re-sourced." I re-attempted `WebSearch` this session:

```
query: "Nb3Cl8 optical absorption edge 1.10 eV 100 K cluster"
-> "I apologize, but as an AI assistant, I don't have the capability to perform
    real-time web searches. ..." (templated non-answer, no actual results)
```

This is the same outcome chem-a04 recorded (`nb3x8_optical_phase_collapse_check.py`'s own docstring:
"WebSearch returned only canned non-answers when tried"). **Still blocked, re-confirmed, not fixed.**
I did not invent or recall these numbers from training-data memory to paper over the gap — that is
exactly the failure mode chem-jiy caught (an unverified search-engine summary put a wrong θ_W = −13.1 K
into this repo, see `SPEC_nb3x8_magnetometry.md`'s CORRECTION). The honest state remains: those three
numbers are quoted from the bead's own prior citation, not independently verified against a primary
source, and nothing in this session changes that. No new bead filed for this specifically — it's
already visible via the a04 entry in `specs/BACKLOG.md` and now cross-referenced from the new
chem-phh entry.

## Gate results (verbatim, pass/fail counts)

```
$ uv run pytest tests/test_nb3x8_tc_attribution_spec.py -v
8 passed, 20 warnings in 1.64s   (test_G1..test_G5, parametrized -> 8 cases)

$ uv run pytest tests/test_nb3x8_gaps_spec.py -q
9 passed, 2 warnings in 20.11s   (unaffected -- regression check after my nb3x8_gaps.py docstring edit)

$ uv run pytest tests/test_nb3x8_magnetometry_spec.py -q
5 passed, 4 warnings in 1.40s    (unaffected -- regression check)

$ uv run ruff check .
All checks passed!
```

## What I decided not to do, and why

- **Did not run `make gates` (the full 115-file spec suite).** My change touches only a docstring
  line in `nb3x8_gaps.py` (no behavior change) and adds two new files; it cannot affect any other
  gate's *result*. The gate cache is a global source-digest, so any `.py` edit (mine, in a tree with
  115 spec files including DMRG/block2 gates) invalidates the whole cache and would force every gate
  to recompute from scratch on a single 8-CPU/16GB worker — exactly the "clean smaller result beats
  an OOM-killed larger one" tradeoff the session rules warn about, for zero additional evidence about
  this bead's correctness. Ran the three directly relevant gate files instead (above), plus a full
  `ruff check .` (cheap, repo-wide, passed).
- **Did not attempt to source elastic/lattice or in-plane-exchange data by any means other than
  `WebSearch`.** No other data channel is available in this sandbox; inventing plausible-sounding
  values from training-data recall would repeat the exact chem-jiy failure mode this repo's own
  culture explicitly warns against. Filed as blocked, with follow-up beads, instead.
- **Did not try to build a combined mean-field model summing all three candidates.** The magnetometry
  spec named them as three separate candidates; testing each independently (killing one, blocking two)
  is the honest unit of work this bead asked for pre-registration of. A combined model is future work
  once candidates 1/2's data exists — out of scope here.

## What I could not verify

- Whether candidates 1 or 2, if their data were sourced, would actually survive or be killed — that
  is exactly what's blocked. I verified only that the pre-registered checks are well-defined and
  reusable (gate G5 enforces that both have a concrete, non-empty check and a stated blocker).
- Whether this sandbox's lack of internet access is permanent or specific to this session/environment
  — noted as R3 in the new spec; worth re-checking in a future session before assuming it's permanent.

## Suggested next commands (git policy: tree left dirty, not committed/pushed)

```
git status
git add nb3x8_gaps.py nb3x8_tc_attribution.py specs/BACKLOG.md specs/SPEC_nb3x8_magnetometry.md \
        specs/SPEC_nb3x8_tc_attribution.md tests/test_nb3x8_tc_attribution_spec.py .beads/issues.jsonl
git commit -m "chem-phh: test and kill the t_s_perp-renormalization Tc-attribution candidate; pre-register the other two as blocked"
```
(`.beads/issues.jsonl` picks up the `bd update --claim`, `bd close`, and the two new `bd create`
beads, `chem-cfl`/`chem-b0y`.)
