# SPEC: Nb3X8 Tc-overprediction attribution — pre-registering the three open candidates

**Status:** IMPLEMENTED — gates G1–G7 green. Candidate 3 is tested (killed, in the narrow form
stated). Candidate 1 is now also tested (chem-cfl, **survives** — see §3b); candidate 2 remains
pre-registered but BLOCKED on data this sandbox has not sourced.

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
  a plausibility judgment call. A quantitative plausibility test for candidate 1 (§3b, chem-cfl),
  using an elastic/lattice energy scale sourced from a second primary reference this sandbox could
  fetch. And an honest, explicit pre-registration of candidate 2's kill condition, for whenever the
  needed data becomes available.
- **What we cannot claim:** candidate 2 (in-plane kagome exchange) is **not tested here** — no
  in-plane exchange constant exists anywhere in `arXiv:2501.10320`'s vendored tables, in
  `arXiv:2412.03360` (candidate 1's source, which reports structural/elastic data only), or
  elsewhere in this repo. Killing candidate 3 in its narrow (same-motif) form does **not** kill the
  broader idea that lattice relaxation affects t_s⊥ — the HT/LT phase ratio shows |t_s⊥| *can* move
  by a large factor under an actual phase change, which is candidate 1's territory, not candidate
  3's. Candidate 1 **surviving** §3b's check is not confirmation that a cooperative structural
  transition *is* the mechanism — see R4 below for exactly what the check does and does not show.
- **Corrected finding, this session (chem-cfl):** chem-phh/chem-a04 concluded "no working internet
  access" from `WebSearch` alone returning a canned refusal every time. `WebFetch` on direct URLs
  (arXiv abstract pages, the `export.arxiv.org/api/query` API, and PDF URLs — whose binary content it
  saves locally for the `Read` tool to extract) **does** work in this sandbox. This is how `arXiv:
  2412.03360` was sourced below; it should be retried before any future session assumes candidate 2
  is permanently blocked.

## 3. Approach

**Reference = the family's own vendored ab-initio table** (`arXiv:2501.10320` Table I/IV, already in
`nb3x8_gaps.NB3X8_CLUSTERS`/`NB3X8_LT_BULK`), used two ways in the same computation:

1. Solve for the multiplier `eta` on t_s⊥ (U0, Us held fixed) such that the already-validated
   `nb3x8_magnetometry.predicted_transition_K` lands on the measured Tc (`EXPERIMENT[name]["Tc_K"]`).
2. Compare `|1 - eta|` to an *independent* slice of the same table: the relative |t_s⊥| difference
   between the bulk and bilayer stacking variants of the same halide's LT (dimerized) structure — an
   empirical measurement of how much this hopping actually varies between closely related variants of
   the same phase, already present in the table for reasons unrelated to this check.

Candidate 2 is pre-registered (description, the check that would kill/support it, and why it is
blocked) in `nb3x8_tc_attribution.UNTESTED_CANDIDATES` — not executed.

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
nb3x8_tc_attribution.UNTESTED_CANDIDATES                         # candidate 2, pre-registered
nb3x8_tc_attribution (CLI)                                       # full report
```

Reuses `nb3x8_magnetometry.predicted_transition_K`/`EXPERIMENT`/`MEV_PER_K`, `nb3x8_gaps.
NB3X8_CLUSTERS`/`NB3X8_LT_BULK`/`NB3X8_HT_BULK`, `nb3x8_susceptibility.n2_spectrum`, and
`nb3x8_thermo.entropy` verbatim; the only new physics-adjacent code is `mean_energy` (the same trace
`entropy` already runs internally, exposed as a public accessor) and the arXiv:2412.03360 constants.

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
- **G4 — honest boundary.** `UNTESTED_CANDIDATES` names exactly candidate 2, with a non-empty
  description, a pre-registered check containing an explicit kill condition, and a stated reason it
  is blocked — so a future edit cannot silently mark it "tested" without the check actually running,
  and cannot leave candidate 1 in this dict now that it has been tested.
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

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_tc_attribution_spec.py` encoding G1–G5 (initially failing — module missing).
2. `nb3x8_tc_attribution.py` (composition of validated primitives + the anchor comparison +
   `UNTESTED_CANDIDATES`).
3. `make gates` (own process; no block2).
4. **chem-cfl addendum:** source `arXiv:2412.03360` (§3b); add `mean_energy`,
   `exfoliation_energy_per_fu`, `required_e_lattice`, `candidate1_verdict`; move
   `cooperative_structural_transition` out of `UNTESTED_CANDIDATES`; extend the gate file to G5–G7.

## 7. Out of scope

- Actually running candidate 2 — it needs an in-plane kagome exchange constant not vendored in this
  repo, in `arXiv:2501.10320`, or in `arXiv:2412.03360`, and not yet sourced from elsewhere.
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
  (§3b, survives); candidate 2 still needs the blocked data.
- **R3 — WebSearch vs WebFetch.** chem-phh/chem-a04's "no working internet access" finding was based
  only on `WebSearch`, which still returns a canned refusal. `WebFetch` on direct URLs (arXiv
  abstract pages, the arXiv API, and PDF URLs) works in this sandbox and is how `arXiv:2412.03360`
  was sourced — a future session should retry it for candidate 2 before assuming that data is
  permanently unreachable.
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

## 9. Deliverables

- `nb3x8_tc_attribution.py` — `required_t_rescale`, `same_motif_t_spread`,
  `max_family_same_motif_t_spread`, `ht_lt_t_ratio`, `candidate3_verdict`, `mean_energy`,
  `exfoliation_energy_per_fu`, `required_e_lattice`, `candidate1_verdict`, `UNTESTED_CANDIDATES`, CLI.
- `tests/test_nb3x8_tc_attribution_spec.py` — gates G1–G7.
- `specs/SPEC_nb3x8_tc_attribution.md` — this spec.
- `specs/BACKLOG.md` — Nb₃X₈ materials line entry recording both verdicts and the WebFetch finding.
- Follow-up bead for candidate 2, filed as out-of-scope work (data sourcing), not done here:
  chem-b0y (in-plane kagome exchange data) — now with a note to retry `WebFetch`, not just
  `WebSearch`.
