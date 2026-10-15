# SPEC: Nb3X8 Tc-overprediction attribution — pre-registering the three open candidates

**Status:** IMPLEMENTED — gates G1–G7 green. Candidates 2 and 3 are tested (both killed, in the
forms stated); candidate 1 is pre-registered but BLOCKED on data this sandbox cannot source.

---

## 1. Goal

`SPEC_nb3x8_magnetometry.md` §G3 found that the isolated Nb₃X₈ interlayer dimer overpredicts the
measured singlet-formation Tc by 5.3× (Nb₃Cl₈) / 2.3× (Nb₃Br₈). Its §7 follow-up named "coordination
/ mean-field reduction" as the likely fix; chem-g78 (`SPEC_nb3x8_gaps.md` G8) tested that and found
it **false** — the spin-channel J_eff does not fall with cluster size for either halide. That left
the attribution open, naming three untested candidates. This spec pre-registers a killable check for
each. Two are now executed:

- **Candidate 3 (chem-phh): a direct lattice renormalization of t_s⊥** — whether shrinking the
  interlayer hopping alone, within the *same* dimerized structural motif, can plausibly account for
  the miss. Falsifiable: if the needed hopping change is no bigger than what this family's own data
  shows between closely related dimerized stacking variants, the candidate survives; if it requires a
  renormalization far outside that range, it is killed.
- **Candidate 2 (chem-b0y): in-plane kagome exchange** — whether a frustrated in-plane coupling
  between Nb₃ trimers, now sourced from a primary ab-initio reference (`arXiv:2508.17352` Table 1),
  can plausibly suppress Tc by the needed 5.3×/2.3× via standard frustration/mean-field physics.
  Falsifiable: if the effective in-plane coupling needed to do the suppression job (via the standard
  RPA self-consistency closure) is unphysically large relative to the sourced value, the candidate is
  killed; if a modest fraction of the sourced coupling suffices, it survives.

## 2. Background and honest framing

- **What we can claim if gates pass:** quantitative kill tests for candidates 2 and 3, each using an
  anchor from real data (bulk-vs-bilayer |t_s⊥| spread for candidate 3; a primary-source in-plane
  exchange constant, verified against the downloaded PDF text per the chem-jiy methodology, for
  candidate 2) — not a new fit, not a plausibility judgment call. And an honest, explicit
  pre-registration of candidate 1's kill condition, for whenever the needed elastic/lattice data
  becomes available.
- **What we cannot claim:** candidate 1 (cooperative/first-order structural transition) is **not
  tested here**. Testing it needs elastic/lattice free-energy data this repo does not have and this
  sandbox cannot fetch (not a search-engine-access problem this time — `arXiv:2501.10320`, the source
  for all cRPA parameters this repo vendors, simply does not report phonon/elastic constants; no other
  attempted source does either). Killing candidate 3 in its narrow (same-motif) form does **not** kill
  the broader idea that lattice relaxation affects t_s⊥ — the HT/LT phase ratio shows |t_s⊥| *can*
  move by a large factor under an actual phase change, which is candidate 1's territory, not candidate
  3's. Killing candidate 2 does not address whether the in-plane coupling, together with a lattice
  relaxation (candidate 1), could jointly explain the miss — only that it cannot do so **alone** at the
  sourced scale.
- **Process note (corrects prior sessions' finding):** chem-phh and chem-a04 both reported "no working
  internet access" and treated candidate 2 as permanently blocked on that basis. That conclusion was
  only ever tested against the `WebSearch` tool, which genuinely is non-functional in this sandbox
  (returns generic non-answers). Direct bash `curl` to `export.arxiv.org` (and following its
  redirect) works and returns real content; the arXiv API plus direct PDF download + text extraction
  (`pypdf`) found and verified the primary source used here. The sandbox's internet access was never
  actually absent — only the AI-level search tool was broken. Candidate 1 remains blocked for a
  different reason (the data genuinely isn't reported anywhere found), not an access problem.

## 3. Approach

**Candidate 3 reference = the family's own vendored ab-initio table** (`arXiv:2501.10320` Table
I/IV, already in `nb3x8_gaps.NB3X8_CLUSTERS`/`NB3X8_LT_BULK`), used two ways in the same computation:

1. Solve for the multiplier `eta` on t_s⊥ (U0, Us held fixed) such that the already-validated
   `nb3x8_magnetometry.predicted_transition_K` lands on the measured Tc (`EXPERIMENT[name]["Tc_K"]`).
2. Compare `|1 - eta|` to an *independent* slice of the same table: the relative |t_s⊥| difference
   between the bulk and bilayer stacking variants of the same halide's LT (dimerized) structure — an
   empirical measurement of how much this hopping actually varies between closely related variants of
   the same phase, already present in the table for reasons unrelated to this check.

**Candidate 2 reference = a sourced primary-source in-plane exchange constant.** Zhang, Dai, Wang,
Cheng & Ji, "Robust Mottness and tunable interlayer magnetism in Nb3X8 (X = F, Cl, Br, I) bilayers"
(`arXiv:2508.17352`), Table 1, gives DFT+U total-energy-mapped in-plane NN (J1) / NNN (J2) exchange
between Nb₃ trimers for monolayer Nb₃Cl₈/Nb₃Br₈ (meV per trimer, "positive = AFM"), verified against
the downloaded PDF text (chem-jiy methodology — never trust a search-engine/AI summary). Two
findings correct the conjecture's premise before any suppression arithmetic: (a) the paper's own
in-plane lattice of Nb₃ trimers is **triangular**, not kagome; (b) the coupling is not robustly
antiferromagnetic — Nb₃Br₈'s J1 is FM by the paper's sign convention, and Nb₃Cl₈'s J1 is a tiny
near-degenerate AFM/FM difference (the paper: "Nb3Cl8 has nearly degenerate FM, ZZ-AFM and
Stripe-AFM configurations"). The quantitative check is still run, steelmanning the conjecture past
this premise problem:

1. Compute the effective in-plane coupling `J'_eff` that the standard RPA/mean-field
   self-consistency closure for weakly-coupled quantum dimers, `J'_eff * chi_dimer(Tc) = 1` (same
   structure that produces the Curie-Weiss law), would require to pull the isolated interlayer
   dimer's susceptibility crossover down to the measured Tc. `chi_dimer` reuses the already-validated
   `nb3x8_susceptibility.susceptibility`.
2. Compare `J'_eff` to the most generous possible reading of the sourced J1/J2: full
   triangular-lattice coordination (z1=z2=6) on the absolute values, summed with no cancellation
   (`sourced_inplane_bound`) — deliberately overstating the physically expected in-plane scale to
   give the candidate its strongest possible case.

Candidate 1 is pre-registered (description, the check that would kill/support it, and why it is
blocked) in `nb3x8_tc_attribution.UNTESTED_CANDIDATES` — not executed.

## 4. Public interface

```
nb3x8_tc_attribution.required_t_rescale(name) -> float          # eta s.t. predicted Tc = measured Tc
nb3x8_tc_attribution.same_motif_t_spread(name) -> float          # bulk-vs-bilayer |t_s_perp| spread
nb3x8_tc_attribution.max_family_same_motif_t_spread() -> float   # family-wide anchor
nb3x8_tc_attribution.ht_lt_t_ratio(name) -> float                # context: HT/LT t ratio (Cl/Br only)
nb3x8_tc_attribution.candidate3_verdict(name, margin=10.0) -> dict
nb3x8_tc_attribution.NB3X8_INPLANE_EXCHANGE                      # sourced J1/J2, arXiv:2508.17352 Table 1
nb3x8_tc_attribution.required_inplane_coupling(name) -> float    # RPA-required J'_eff (meV)
nb3x8_tc_attribution.sourced_inplane_bound(name, z1=6, z2=6) -> float  # most generous sourced bound
nb3x8_tc_attribution.candidate2_verdict(name) -> dict
nb3x8_tc_attribution.UNTESTED_CANDIDATES                         # candidate 1, pre-registered
nb3x8_tc_attribution (CLI)                                       # full report
```

Reuses `nb3x8_magnetometry.predicted_transition_K`/`EXPERIMENT`/`MEV_PER_K`,
`nb3x8_gaps.NB3X8_CLUSTERS`/`NB3X8_LT_BULK`, `nb3x8_susceptibility.susceptibility`, and
`odmd_spin.dimer_exchange_analytic` verbatim; no new physics primitives, only the sourced table.

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
- **G5 — honest boundary.** `UNTESTED_CANDIDATES` names exactly candidate 1, with a non-empty
  description, a pre-registered check containing an explicit kill condition, and a stated reason it
  is blocked; candidate 2 no longer appears there — so a future edit cannot silently mark either
  "tested" without the check actually running, nor silently re-mark candidate 2 "untested".
- **G6 — candidate 2 sourced-data regression pin.** `NB3X8_INPLANE_EXCHANGE`'s J1/J2 match the
  primary source exactly (0.12/-0.24 meV Cl, -2.74/-0.37 meV Br).
- **G7 — candidate 2 (DEFINITION OF DONE).** The RPA-required `J'_eff` (19776 meV Cl, 663 meV Br)
  exceeds the most generous sourced bound (2.16 meV Cl, 18.66 meV Br) by 9156× (Cl) / 35.5× (Br), and
  exceeds `J_interlayer` itself (66.2 meV Cl, 119.1 meV Br) by 298.7× (Cl) / 5.6× (Br) — candidate 2
  is **KILLED** for both halides, well past the pre-registered bar ("requires J' ... exceeding
  J_interlayer").

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_tc_attribution_spec.py` encoding G1–G7 (candidate-2 gates initially failing —
   sourced data and functions missing).
2. `nb3x8_tc_attribution.py` (composition of validated primitives + the anchor comparison + sourced
   in-plane exchange table + `UNTESTED_CANDIDATES` reduced to candidate 1).
3. `make gates` (own process; no block2).

## 7. Out of scope

- Actually running candidate 1 — it needs elastic/lattice free-energy data not vendored in this repo
  and not found in any source checked this session.
- A quantitative mean-field model combining candidate 1 with the (killed) candidates 2/3 — moot now
  that 2 and 3 are killed; candidate 1 remains the sole untested mechanism.
- Revisiting the "kagome" terminology in the bead/backlog history: the sourced primary literature
  calls the in-plane net triangular, not kagome; this spec notes the correction (§2) but does not
  rename the historical bead/candidate labels.

## 8. Caveats and risks

- **R1 — the margin (10×) is a generous kill bar, not a tight one (candidate 3).** The needed
  rescaling (36–58%) is 17–28× the 2.04% anchor even at this family's *largest* same-motif spread
  (Nb₃F₈, a different halide from the ones being tested) — the conclusion is robust to the exact
  margin chosen; a 2× or 20× margin gives the same verdict.
- **R2 — "killed" is narrow (candidate 3).** This only rules out a same-structural-motif
  renormalization of t_s⊥ alone. It does not address t_s⊥ changing as part of a genuine structural
  phase shift (candidate 1).
- **R3 — "killed" is also narrow for candidate 2.** The RPA closure tests whether the in-plane
  coupling, acting *alone* via the standard weak-coupling mean-field mechanism, can explain the full
  miss. It does not rule out a small in-plane contribution combined with a lattice relaxation
  (candidate 1) — only that candidate 2 cannot be the *whole* explanation at the sourced coupling
  scale. The margin is large (35.5×–9156× past the generous bound, and the needed coupling exceeds
  even `J_interlayer` itself), so this is not a close call.
- **R4 — single-paper sourcing for candidate 2.** Only one primary-source DFT+U paper
  (`arXiv:2508.17352`) was found and used for the in-plane J1/J2; no second independent calculation
  was cross-checked (unlike the interlayer J, which multiple specs anchor to `arXiv:2501.10320`). The
  kill margin (>35×) is large enough that a factor-of-few revision to J1/J2 would not change the
  verdict, but this is a single-source number.
- **R5 — the earlier "no working internet access" framing was wrong for candidate 2, right for
  candidate 1.** §2 records the correction: `WebSearch` is broken in this sandbox, but bash `curl` to
  `export.arxiv.org` works. Candidate 1 is blocked because the data isn't reported in any source
  found, not because of a tooling limitation — worth re-distinguishing in future sessions.

## 9. Deliverables

- `nb3x8_tc_attribution.py` — `required_t_rescale`, `same_motif_t_spread`,
  `max_family_same_motif_t_spread`, `ht_lt_t_ratio`, `candidate3_verdict`,
  `NB3X8_INPLANE_EXCHANGE`, `required_inplane_coupling`, `sourced_inplane_bound`,
  `candidate2_verdict`, `UNTESTED_CANDIDATES`, CLI.
- `tests/test_nb3x8_tc_attribution_spec.py` — gates G1–G7.
- `specs/SPEC_nb3x8_tc_attribution.md` — this spec.
- `specs/BACKLOG.md` — Nb₃X₈ materials line entry recording the candidate-2 kill verdict and the
  curl-vs-WebSearch process correction.
- Follow-up bead for candidate 1 only, filed as out-of-scope work (data sourcing), not done here:
  chem-cfl (structural-transition elastic/lattice data). chem-b0y (this bead) resolves candidate 2
  and is closed, not left open as a follow-up.
