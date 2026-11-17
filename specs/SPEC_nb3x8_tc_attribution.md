# SPEC: Nb3X8 Tc-overprediction attribution — pre-registering the three open candidates

**Status:** IMPLEMENTED — gates G1–G9 green. All three candidates are tested: candidate 3 killed
(chem-phh), candidate 1 **survives** (chem-cfl, §3b), candidate 2 killed (chem-b0y, §3c) — each in the
narrow form stated.

---

## 1. Goal

`SPEC_nb3x8_magnetometry.md` §G3 found that the isolated Nb₃X₈ interlayer dimer overpredicts the
measured singlet-formation Tc by 5.3× (Nb₃Cl₈) / 2.3× (Nb₃Br₈). Its §7 follow-up named "coordination
/ mean-field reduction" as the likely fix; chem-g78 (`SPEC_nb3x8_gaps.md` G8) tested that and found
it **false** — the spin-channel J_eff does not fall with cluster size for either halide. That left
the attribution open, naming three untested candidates. This spec pre-registers a killable check for
each. All three are now executed:

- **Candidate 3 (chem-phh): a direct lattice renormalization of t_s⊥** — whether shrinking the
  interlayer hopping alone, within the *same* dimerized structural motif, can plausibly account for
  the miss. Falsifiable: if the needed hopping change is no bigger than what this family's own data
  shows between closely related dimerized stacking variants, the candidate survives; if it requires a
  renormalization far outside that range, it is killed.
- **Candidate 1 (chem-cfl): the cooperative/first-order structural transition** — whether a
  physically plausible lattice energy, anchored to a second primary reference (`arXiv:2412.03360`),
  can move the HT/LT free-energy crossing onto the measured Tc (§3b).
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
  candidate 2) — not a new fit, not a plausibility judgment call. And a quantitative plausibility test for candidate 1 (§3b, chem-cfl), using an elastic/lattice
  energy scale from `arXiv:2412.03360`.
- **What we cannot claim:** candidate 1 **surviving** §3b's check is not confirmation that a cooperative structural
  transition *is* the mechanism — see R4 for exactly what the check does and does not show. Killing candidate 3 in its narrow (same-motif) form does **not** kill
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
  actually absent — only the AI-level search tool was broken. chem-cfl reached the same correction independently, with `WebFetch` on direct URLs, and sourced
  candidate 1's data that way (§3b).

## 3. Approach

**Candidate 3 reference = the family's own vendored ab-initio table** (`arXiv:2501.10320` Table
I/IV, already in `nb3x8_gaps.NB3X8_CLUSTERS`/`NB3X8_LT_BULK`), used two ways in the same computation:

1. Solve for the multiplier `eta` on t_s⊥ (U0, Us held fixed) such that the already-validated
   `nb3x8_magnetometry.predicted_transition_K` lands on the measured Tc (`EXPERIMENT[name]["Tc_K"]`).
2. Compare `|1 - eta|` to an *independent* slice of the same table: the relative |t_s⊥| difference
   between the bulk and bilayer stacking variants of the same halide's LT (dimerized) structure — an
   empirical measurement of how much this hopping actually varies between closely related variants of
   the same phase, already present in the table for reasons unrelated to this check.

## 3b. Approach — candidate 1 (chem-cfl, executed this session)

**New data sourced:** Mortazavi, Zhuang & Rabczuk, "A first-principles study on the physical
properties of two-dimensional Nb3Cl8, Nb3Br8 and Nb3I8," `arXiv:2412.03360v1`. Fetched via `WebFetch`
on the arXiv PDF URL (binary saved locally, full text then read with the `Read` tool — multi-page
PDF text extraction, not an abstract/summary) and **verified against the primary text directly**,
per the chem-jiy lesson (`SPEC_nb3x8_magnetometry.md` G5/R1: never carry forward an unverified
search-engine number). Independent of `arXiv:2501.10320` — DFT (VASP, PBE-GGA+D3) structural/elastic
properties only, no electronic or magnetic parameters, and it resolves only the single known (LT,
dimerized) ground-state structure, no HT phase.

Verbatim numbers used (`nb3x8_tc_attribution.NB3X8_EXFOLIATION`): monolayer exfoliation energy
γ = 0.24 J/m² (Nb₃Cl₈) / 0.27 J/m² (Nb₃Br₈); monolayer hexagonal lattice constant a = 6.783 Å
(Nb₃Cl₈) / 7.121 Å (Nb₃Br₈).

**The pre-registered check, executed exactly:**

1. `F_phase(T) = E_mag,phase(T) - T·S_mag,phase(T) + E_lattice,phase` for phase ∈ {HT, LT}, with
   E_mag/S_mag the exact N=2 Boltzmann trace already validated by `nb3x8_susceptibility.n2_spectrum`
   / `nb3x8_thermo.entropy` (new accessor `nb3x8_tc_attribution.mean_energy` exposes `<E>(T)`, the
   same trace `entropy` already computes internally), evaluated on `NB3X8_HT_BULK` vs
   `NB3X8_LT_BULK`.
2. Because `E_lattice` enters `F_HT(T) - F_LT(T)` linearly and only once, the value that makes the
   crossing land exactly on the measured Tc is solved directly (no root-finding):
   `E_lattice = Tc·ΔS_mag(Tc) − ΔE_mag(Tc)` (`required_e_lattice`).
3. The physically-plausible range for `E_lattice` is anchored to the exfoliation energy above,
   converted to meV per Nb₃X₈ formula unit assuming the standard 1-formula-unit hexagonal monolayer
   cell (`exfoliation_energy_per_fu`; area = (√3/2)·a²) — used as a loose **upper** bound: a partial
   HT↔LT structural relaxation within the bound interlayer well cannot plausibly cost more than
   separating the layers to infinity.
4. **KILL** if the solved `E_lattice` is nonpositive (wrong sign) or exceeds the exfoliation-energy
   bound; **SURVIVES** otherwise (`candidate1_verdict`).

**Result:** `E_lattice` needed = 9.94 meV (Cl, 1.7% of the 596.9 meV bound) / 41.61 meV (Br, 5.6% of
the 740.1 meV bound) — both positive, both a small fraction of the bound. **Candidate 1 SURVIVES**:
a modest, physically plausible lattice energy reproduces the measured Tc without extreme tuning. See
R4 for what this does and does not establish.

## 3c. Approach — candidate 2 (chem-b0y)

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

**Result:** see G9 — candidate 2 is **KILLED** for both halides.

## 4. Public interface

```
nb3x8_tc_attribution.required_t_rescale(name) -> float          # eta s.t. predicted Tc = measured Tc
nb3x8_tc_attribution.same_motif_t_spread(name) -> float          # bulk-vs-bilayer |t_s_perp| spread
nb3x8_tc_attribution.max_family_same_motif_t_spread() -> float   # family-wide anchor
nb3x8_tc_attribution.ht_lt_t_ratio(name) -> float                # context: HT/LT t ratio (Cl/Br only)
nb3x8_tc_attribution.candidate3_verdict(name, margin=10.0) -> dict
nb3x8_tc_attribution.mean_energy(U0, t, Us, T) -> float          # <E>(T), N=2 Boltzmann trace
nb3x8_tc_attribution.exfoliation_energy_per_fu(name) -> float    # meV, plausible E_lattice bound
nb3x8_tc_attribution.required_e_lattice(name) -> float           # meV, E_lattice at the measured Tc
nb3x8_tc_attribution.candidate1_verdict(name) -> dict
nb3x8_tc_attribution.NB3X8_INPLANE_EXCHANGE                      # sourced J1/J2, arXiv:2508.17352 Table 1
nb3x8_tc_attribution.required_inplane_coupling(name) -> float    # RPA-required J'_eff (meV)
nb3x8_tc_attribution.sourced_inplane_bound(name, z1=6, z2=6) -> float  # most generous sourced bound
nb3x8_tc_attribution.candidate2_verdict(name) -> dict
nb3x8_tc_attribution.UNTESTED_CANDIDATES                         # empty: all three candidates tested
nb3x8_tc_attribution (CLI)                                       # full report
```

Reuses `nb3x8_magnetometry.predicted_transition_K`/`EXPERIMENT`/`MEV_PER_K`, `nb3x8_gaps.
NB3X8_CLUSTERS`/`NB3X8_LT_BULK`/`NB3X8_HT_BULK`, `nb3x8_susceptibility.n2_spectrum`/`susceptibility`,
`nb3x8_thermo.entropy` and `odmd_spin.dimer_exchange_analytic` verbatim; the only new
physics-adjacent code is `mean_energy` (the trace `entropy` already runs internally, exposed as a
public accessor) and the two sourced tables (arXiv:2412.03360, arXiv:2508.17352).

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
- **G4 — honest boundary.** `UNTESTED_CANDIDATES` is empty: every candidate has a verdict
  function, so none may also be listed as untested (chem-cfl and chem-b0y each removed the other's
  candidate in parallel; the merge pins the combined state).
- **G5 — candidate 1 (DEFINITION OF DONE, chem-cfl).** For both Cl and Br, `required_e_lattice` is
  positive (correct sign) and less than 10% of `exfoliation_energy_per_fu` (the arXiv:2412.03360
  bound) — candidate 1 **SURVIVES**, and by a comfortable margin (robust even to a 2× tightening of
  the bound, per R4).
- **G6 — self-consistency.** Plugging the solved `E_lattice` back into `F_HT(Tc) - F_LT(Tc)`
  reproduces zero (the measured-Tc crossing condition) to <1e-8 meV.
- **G7 — unique crossing.** The `E_lattice` needed to put the crossing at temperature T is
  non-decreasing in T across `[0.01, 3000]` meV (spanning both measured Tc values by orders of
  magnitude) — so a fixed `E_lattice` cannot cross `F_HT - F_LT` more than once, and "the" crossing
  temperature in step 2 of §3b is unambiguous.
- **G8 — candidate 2 sourced-data regression pin.** `NB3X8_INPLANE_EXCHANGE`'s J1/J2 match the
  primary source exactly (0.12/-0.24 meV Cl, -2.74/-0.37 meV Br).
- **G9 — candidate 2 (DEFINITION OF DONE).** The RPA-required `J'_eff` (19776 meV Cl, 663 meV Br)
  exceeds the most generous sourced bound (2.16 meV Cl, 18.66 meV Br) by 9156× (Cl) / 35.5× (Br), and
  exceeds `J_interlayer` itself (66.2 meV Cl, 119.1 meV Br) by 298.7× (Cl) / 5.6× (Br) — candidate 2
  is **KILLED** for both halides, well past the pre-registered bar ("requires J' ... exceeding
  J_interlayer").

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_tc_attribution_spec.py` encoding G1–G5 (initially failing — module missing).
2. `nb3x8_tc_attribution.py` (composition of validated primitives + the anchor comparison +
   `UNTESTED_CANDIDATES`).
3. `make gates` (own process; no block2).
4. **chem-cfl addendum:** source `arXiv:2412.03360` (§3b); add `mean_energy`,
   `exfoliation_energy_per_fu`, `required_e_lattice`, `candidate1_verdict`; move
   `cooperative_structural_transition` out of `UNTESTED_CANDIDATES`; extend the gate file to G5–G7.
5. **chem-b0y addendum** (run in parallel with chem-cfl, both from chem-phh): source
   `arXiv:2508.17352` (§3c); add `NB3X8_INPLANE_EXCHANGE`, `required_inplane_coupling`,
   `sourced_inplane_bound`, `candidate2_verdict`; gates G8–G9. Merged on main: `UNTESTED_CANDIDATES`
   emptied, G4 pins that.

## 7. Out of scope

- Revisiting the "kagome" terminology in the bead/backlog history: the sourced primary literature
  calls the in-plane net triangular, not kagome; this spec notes the correction (§3c) but does not
  rename the historical bead/candidate labels.
- Re-sourcing `arXiv:2501.10320` for the needed data (candidate 1 used a *different* primary
  reference instead, since the vendored one has no elastic/lattice data at all).
- A quantitative mean-field model combining all three candidates — this spec tests them
  independently, as the magnetometry spec's §7 named them.
- Resolving the true HT-phase lattice structure (phonons, actual interlayer relaxation distance) —
  `arXiv:2412.03360` reports only the known LT structure; the exfoliation energy is used as an upper
  *bound*, not a measurement of the HT/LT structural energy difference (R4).

## 8. Caveats and risks

- **R1 — the margin (10×) is a generous kill bar, not a tight one.** The needed rescaling (36–58%)
  is 17–28× the 2.04% anchor even at this family's *largest* same-motif spread (Nb₃F₈, a different
  halide from the ones being tested) — the conclusion is robust to the exact margin chosen; a 2× or
  20× margin gives the same verdict.
- **R2 — "killed" is narrow.** This only rules out a same-structural-motif renormalization of t_s⊥
  alone. It does not address t_s⊥ changing as part of a genuine structural phase shift (candidate 1)
  or in combination with a competing in-plane coupling (candidate 2) — candidate 1 is now tested
  (§3b, survives); candidate 2 is now tested too (§3c, killed).
- **R3 — WebSearch vs WebFetch.** chem-phh/chem-a04's "no working internet access" finding was based
  only on `WebSearch`, which still returns a canned refusal. `WebFetch` on direct URLs (arXiv
  abstract pages, the arXiv API, and PDF URLs) works in this sandbox and is how `arXiv:2412.03360`
  was sourced; chem-b0y found the same with bash `curl` to `export.arxiv.org` and sourced
  `arXiv:2508.17352` that way.
- **R4 — candidate 1 "survives" is a plausibility bound, not confirmation.** Three honest
  limitations, stated explicitly rather than rounded into the "survives" verdict: (a)
  `exfoliation_energy_per_fu` is the full layer-separation energy, a *generous upper bound* on a
  partial HT↔LT structural relaxation, not a measured HT/LT structural energy difference —
  `arXiv:2412.03360` does not resolve an HT phase at all, so no paper directly measures the quantity
  this check needs. (b) The bound assumes 1 formula unit per hexagonal monolayer cell; even a 2×
  tightening (2 formula units/cell) leaves the needed `E_lattice` (9.9/41.6 meV) far below the
  halved bound (298/370 meV), so the verdict is robust to this specific ambiguity, but it is still an
  assumption, not an independently reverified atom count. (c) `NB3X8_HT_BULK`/`NB3X8_LT_BULK` are
  single-dimer electronic parameters with no lattice feedback — the check asks "is a lattice energy
  of this size plausible," not "does this specific lattice model predict Tc=90 K/382 K ab initio."
  Surviving this check rules out only that candidate 1 is *implausible on energetic grounds*; it is
  not a demonstration that the cooperative structural transition is the correct mechanism.
- **R5 — "killed" is also narrow for candidate 2.** The RPA closure tests whether the in-plane
  coupling, acting *alone* via the standard weak-coupling mean-field mechanism, can explain the full
  miss. It does not rule out a small in-plane contribution combined with a lattice relaxation
  (candidate 1) — only that candidate 2 cannot be the *whole* explanation at the sourced coupling
  scale. The margin is large (35.5×–9156× past the generous bound, and the needed coupling exceeds
  even `J_interlayer` itself), so this is not a close call.
- **R6 — single-paper sourcing for candidate 2.** Only one primary-source DFT+U paper
  (`arXiv:2508.17352`) was found and used for the in-plane J1/J2; no second independent calculation
  was cross-checked (unlike the interlayer J, which multiple specs anchor to `arXiv:2501.10320`). The
  kill margin (>35×) is large enough that a factor-of-few revision to J1/J2 would not change the
  verdict, but this is a single-source number.

## 9. Deliverables

- `nb3x8_tc_attribution.py` — `required_t_rescale`, `same_motif_t_spread`,
  `max_family_same_motif_t_spread`, `ht_lt_t_ratio`, `candidate3_verdict`, `mean_energy`,
  `exfoliation_energy_per_fu`, `required_e_lattice`, `candidate1_verdict`, `NB3X8_INPLANE_EXCHANGE`,
  `required_inplane_coupling`, `sourced_inplane_bound`, `candidate2_verdict`, `UNTESTED_CANDIDATES`,
  CLI.
- `tests/test_nb3x8_tc_attribution_spec.py` — gates G1–G9.
- `specs/SPEC_nb3x8_tc_attribution.md` — this spec.
- `specs/BACKLOG.md` — Nb₃X₈ materials line entries recording all three verdicts and the
  network-access correction.
- No follow-up bead: chem-cfl (candidate 1) and chem-b0y (candidate 2) each resolve their candidate.
