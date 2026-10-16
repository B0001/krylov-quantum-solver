# SPEC: Nb3X8 Tc-overprediction attribution — pre-registering the three open candidates

**Status:** IMPLEMENTED — gates G1–G5 green. Candidate 3 is tested (killed, in the narrow form
stated); candidates 1 and 2 are pre-registered but BLOCKED on data this sandbox cannot source.

---

## 1. Goal

`SPEC_nb3x8_magnetometry.md` §G3 found that the isolated Nb₃X₈ interlayer dimer overpredicts the
measured singlet-formation Tc by 5.3× (Nb₃Cl₈) / 2.3× (Nb₃Br₈). Its §7 follow-up named "coordination
/ mean-field reduction" as the likely fix; chem-g78 (`SPEC_nb3x8_gaps.md` G8) tested that and found
it **false** — the spin-channel J_eff does not fall with cluster size for either halide. That left
the attribution open, naming three untested candidates. This spec pre-registers a killable check for
each, and executes the one of the three that existing vendored data supports: **a direct lattice
renormalization of t_s⊥** — specifically, whether shrinking the interlayer hopping alone, within the
*same* dimerized structural motif, can plausibly account for the miss. The claim is falsifiable: if
the needed hopping change is no bigger than what this family's own data shows between closely related
dimerized stacking variants, the candidate survives; if it requires a renormalization far outside that
range, it is killed.

## 2. Background and honest framing

- **What we can claim if gates pass:** a quantitative kill test for candidate 3, using an anchor
  (bulk-vs-bilayer |t_s⊥| spread) that is already in the vendored ab-initio table — not a new fit, not
  a plausibility judgment call. And an honest, explicit pre-registration of candidates 1 and 2's kill
  conditions, for whenever the needed data becomes available.
- **What we cannot claim:** candidates 1 (cooperative/first-order structural transition) and 2
  (in-plane kagome exchange) are **not tested here**. Testing them needs data this repo does not
  have and this sandbox cannot fetch (no elastic/lattice free energies; no in-plane exchange
  constant anywhere in `arXiv:2501.10320`'s vendored tables; no working internet access in this
  sandbox — re-confirmed this session, matching chem-a04's prior finding in
  `nb3x8_optical_phase_collapse_check.py`). Killing candidate 3 in its narrow (same-motif) form does
  **not** kill the broader idea that lattice relaxation affects t_s⊥ — the HT/LT phase ratio shows
  |t_s⊥| *can* move by a large factor under an actual phase change, which is candidate 1's territory,
  not candidate 3's.

## 3. Approach

**Reference = the family's own vendored ab-initio table** (`arXiv:2501.10320` Table I/IV, already in
`nb3x8_gaps.NB3X8_CLUSTERS`/`NB3X8_LT_BULK`), used two ways in the same computation:

1. Solve for the multiplier `eta` on t_s⊥ (U0, Us held fixed) such that the already-validated
   `nb3x8_magnetometry.predicted_transition_K` lands on the measured Tc (`EXPERIMENT[name]["Tc_K"]`).
2. Compare `|1 - eta|` to an *independent* slice of the same table: the relative |t_s⊥| difference
   between the bulk and bilayer stacking variants of the same halide's LT (dimerized) structure — an
   empirical measurement of how much this hopping actually varies between closely related variants of
   the same phase, already present in the table for reasons unrelated to this check.

Candidates 1 and 2 are pre-registered (description, the check that would kill/support them, and why
they are blocked) in `nb3x8_tc_attribution.UNTESTED_CANDIDATES` — not executed.

## 4. Public interface

```
nb3x8_tc_attribution.required_t_rescale(name) -> float          # eta s.t. predicted Tc = measured Tc
nb3x8_tc_attribution.same_motif_t_spread(name) -> float          # bulk-vs-bilayer |t_s_perp| spread
nb3x8_tc_attribution.max_family_same_motif_t_spread() -> float   # family-wide anchor
nb3x8_tc_attribution.ht_lt_t_ratio(name) -> float                # context: HT/LT t ratio (Cl/Br only)
nb3x8_tc_attribution.candidate3_verdict(name, margin=10.0) -> dict
nb3x8_tc_attribution.UNTESTED_CANDIDATES                         # candidates 1 & 2, pre-registered
nb3x8_tc_attribution (CLI)                                       # full report
```

Reuses `nb3x8_magnetometry.predicted_transition_K`/`EXPERIMENT` and `nb3x8_gaps.NB3X8_CLUSTERS`/
`NB3X8_LT_BULK` verbatim; no new physics primitives.

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_tc_attribution_spec.py`.

- **G1 — machinery anchor.** At `eta=1` (no rescaling), the module's own Tc-ratio reproduces the
  already-gated overprediction factors (5.3× Cl / 2.3× Br, matching
  `nb3x8_magnetometry.overprediction_factor` to <1e-9) — if this drifts, nothing below is trustworthy.
- **G2 — self-consistency.** `required_t_rescale` round-trips: plugging `eta` back into
  `predicted_transition_K` reproduces the measured Tc to <1e-6 relative, for both halides, with
  `0 < eta < 1` (a genuine reduction, matching the overprediction's sign).
- **G3 — anchor regression pin.** The bulk-vs-bilayer |t_s⊥| spreads (2.041% F, 0.147% Cl, 0.118% Br,
  0.092% I) are pinned against the vendored table.
- **G4 — candidate 3 (DEFINITION OF DONE).** The rescaling needed to match measured Tc (57.7% Cl,
  35.8% Br) exceeds 10× the family-wide same-motif anchor (2.04%, i.e. a 20.4% bound) for **both**
  halides — candidate 3 is **KILLED** in the narrow (same-structural-motif) form stated. Also checks
  that the HT/LT phase ratio moves |t_s⊥| by more than the needed same-motif rescaling — i.e. killing
  candidate 3 does not imply anything about candidate 1.
- **G5 — honest boundary.** `UNTESTED_CANDIDATES` names exactly candidates 1 and 2, each with a
  non-empty description, a pre-registered check containing an explicit kill condition, and a stated
  reason it is blocked — so a future edit cannot silently mark either "tested" without the check
  actually running.

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_tc_attribution_spec.py` encoding G1–G5 (initially failing — module missing).
2. `nb3x8_tc_attribution.py` (composition of validated primitives + the anchor comparison +
   `UNTESTED_CANDIDATES`).
3. `make gates` (own process; no block2).

## 7. Out of scope

- Actually running candidates 1 or 2 — they need data (elastic/lattice free energies; in-plane kagome
  exchange) not vendored in this repo and not fetchable from this sandbox.
- Re-sourcing `arXiv:2501.10320` or primary literature for the needed data (filed as follow-up beads,
  not attempted here).
- A quantitative mean-field model combining all three candidates — this spec tests them
  independently, as the magnetometry spec's §7 named them.

## 8. Caveats and risks

- **R1 — the margin (10×) is a generous kill bar, not a tight one.** The needed rescaling (36–58%)
  is 17–28× the 2.04% anchor even at this family's *largest* same-motif spread (Nb₃F₈, a different
  halide from the ones being tested) — the conclusion is robust to the exact margin chosen; a 2× or
  20× margin gives the same verdict.
- **R2 — "killed" is narrow.** This only rules out a same-structural-motif renormalization of t_s⊥
  alone. It does not address t_s⊥ changing as part of a genuine structural phase shift (candidate 1)
  or in combination with a competing in-plane coupling (candidate 2) — seeing both killed or surviving
  needs the blocked data.
- **R3 — the no-internet-access finding is sandbox-specific**, not a permanent property of the task;
  re-attempted this session (`WebSearch`) and got the same canned non-answer chem-a04 recorded — worth
  re-checking in a future session with different sandbox network access before assuming candidates 1/2
  are permanently blocked.

## 9. Deliverables

- `nb3x8_tc_attribution.py` — `required_t_rescale`, `same_motif_t_spread`,
  `max_family_same_motif_t_spread`, `ht_lt_t_ratio`, `candidate3_verdict`, `UNTESTED_CANDIDATES`, CLI.
- `tests/test_nb3x8_tc_attribution_spec.py` — gates G1–G5.
- `specs/SPEC_nb3x8_tc_attribution.md` — this spec.
- `specs/BACKLOG.md` — Nb₃X₈ materials line entry recording the verdict and the a04 re-sourcing
  re-confirmation.
- Follow-up beads for candidates 1 and 2, filed as out-of-scope work (data sourcing), not done here:
  chem-cfl (structural-transition elastic/lattice data) and chem-b0y (in-plane kagome exchange data).
