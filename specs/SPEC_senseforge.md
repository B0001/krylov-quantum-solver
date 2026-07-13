# SPEC: SenseForge — Nb₃X₈ strain/field sensor screening

**Slug:** `senseforge` · **Status: BUILT (harness), PREMISE BOUNDED (screening).** Gates G1–G9 in
`tests/test_senseforge_spec.py` (30 green). Full PRD `specs/full/spec-nb3x8-sensor-designer.md`,
tasks `specs/tasks/04-senseforge.md`.

> **Supersedes the earlier "BOUNDARY RECORDED — not built" version of this spec.** That disposition
> was written when only the *certified* FoM had been tried. The harness was subsequently built on
> the point-estimate path (which that same note called "viable"). It is real, gated, and useful —
> but the screening premise turns out to be bounded for a *second, independent* reason, recorded in
> §3. The spec had drifted badly out of sync with the code (it said "not built" while a 1,100-line
> package with 25 passing gates sat in the tree); this rewrite closes that gap.

## 1. Goal

Rank Nb₃X₈ strain/field sensor operating points by a figure of merit
`FoM = |sensitivity| / gap-bracket-width`, emitting a ranked `candidates.md` and per-candidate
design cards, each carrying the ADR-0003 cluster-model caveat.

## 2. What was built (and is gated)

A complete, resumable sweep harness on the validated Nb₃X₈ **downfolded interlayer dimer**:

- `senseforge.config` — one validated YAML per run; bad fields fail *by name*; the resolved config
  is SHA-256 hashed into every artifact header (provenance).
- `senseforge.hamiltonian` — strain (`t(ε) = t₀(1+ε)`) and Zeeman (`H_Z = g·μ_B·B·Sz`) perturbations.
- `senseforge.sweep` — grid sweep, **crash-resumable** (kill mid-sweep, resume, get a byte-identical
  CSV; cached points are not recomputed).
- `senseforge.sensitivity` — central-difference `S = d(gap)/dx` with the bracket propagated through.
- `senseforge.candidates` — FoM ranking, `candidates.md`, design cards.
- `senseforge.validation` — closed-form vs exact-diagonalization cross-check on all four halides.

**Three deviations from the PRD, each recorded rather than hidden** (details in
`senseforge/hamiltonian.py`):

1. **No CIF/ab-initio strain geometry.** The repo has no strained-geometry path; strain ε is
   *defined* as the fractional hopping perturbation, reusing `nb3x8_strain.py`'s own "|t| is the
   sole strain proxy" convention. Not an engineering strain with an elastic-tensor mapping.
2. **`certified_gaps.gap_bracket` is not used — it is *actively wrong* here.** The singlet→triplet
   gap J is a **dark excitation** from the closed-shell HF reference (Sz is a good quantum number),
   so `gap_bracket` silently returns ~1117 meV (the *bright* ionic gap) instead of J ≈ 66 meV. Every
   SenseForge gap is exact diagonalization / closed form, wrapped in certchem's contract as an
   explicit **zero-width exact bracket** — honestly labeled `exact`, never as a rigorous bound on an
   approximation.
3. **A uniform field is invisible inside the HF sector** (Sz_tot ≡ 0 there), which is *why* routing
   the field axis through `gap_bracket` could never have worked.

## 3. THE FINDING — the FoM ranking does not discriminate on this model

The screening premise is **vacuous on the only Nb₃X₈ model this repo can reach**, on *both* axes,
for two different reasons. Because the brackets are exact (zero width), `FoM ≡ |S|` — and:

| axis | behaviour of \|S\| | what "rank 1" actually is |
|---|---|---|
| **field** | **exactly constant** — the Zeeman response is exactly linear (G2), so `d(gap)/dB = −g·μ_B` at *every* B | **a tie.** All 19 operating points have identical FoM (0.1158 meV/T). The rank order is **sort noise**. |
| **strain** | **strictly monotone** — J(t) is smooth over the accessible range (\|S\| spans only 123.3 → 126.7, a 2.8% spread) | **the window edge.** Widen the sweep from ±2% to ±5% and "rank 1" obediently moves from +1.75% to +4.75%. |

**Neither axis yields an interior optimum.** A screener that always returns either a tie or the edge
of whatever window you handed it has not screened anything.

**This was shipping as a real recommendation.** The pre-fix `results/senseforge/Nb3Cl8_field/design_card_1.md`
published **"+2 T"** as the #1 operating point — a point provably no better than any other in the
sweep. Fixed: `candidates.ranking_verdict()` classifies every ranking as `degenerate` / `monotone` /
`interior`, and **every artifact now leads with the verdict** ("NO DISCRIMINATION", "NO INTERIOR
OPTIMUM"). Gated by G9 so it cannot silently regress.

**Second defect fixed — false provenance.** `SweepConfig.krylov_dim` was dead (nothing read it;
every `Certificate` is built with `krylov_dim=0` because no Krylov subspace is ever constructed) and
its docstring claimed a use in `validation.py` that did not exist — yet it was hashed into
`content_hash()` and **stamped on every published design card as `krylov_dim=12`**, advertising a
Krylov dimension for an exact-diagonalization result. Removed; a config that still sets it now fails
with the named field.

## 4. Honest disposition

- **The harness is worth keeping.** Config validation, provenance hashing, crash-resume, ADR-0003
  header automation, and the closed-form-vs-exact cross-check are all real and gated. Point them at
  a model with a *non-trivial* sensitivity landscape and they work.
- **The screening premise needs a different system.** It requires a model where (a) the gap bracket
  is tight-but-nontrivial (so the `/width` term does work — needs a system too large for FCI yet
  well-resolved by Krylov), **and** (b) |S| has an interior optimum (so the ranking discriminates).
  The Nb₃X₈ dimer provides neither: it is FCI-trivial, and its response is linear (field) or
  monotone (strain).
- The *sensitivity itself is real and reproduces `nb3x8_strain.py`'s gated Grüneisen result*
  (dJ/d|t| > 0 for every halide, G2) — SenseForge measures a true quantity. It just cannot *rank* on it here.

## 5. Acceptance gates (`tests/test_senseforge_spec.py`, 30 green)

- **G1** config schema: defaults resolve, bad fields fail by name, hash is deterministic + sensitive.
- **G2** physics: ε=0 reproduces the gated `dimer_exchange_analytic` J exactly; strain direction
  matches the Grüneisen sign; Zeeman term is Hermitian; field response is exactly linear below the
  (~572 T) level crossing; a uniform field has **zero** effect inside the HF sector.
- **G3** sweep: one row per grid point; **kill mid-sweep + resume ⇒ identical CSV**; cached points
  are not recomputed.
- **G4** every artifact carries the exact ADR-0003 note.
- **G5** finite differences: the propagated slope bracket contains the analytic slope on a synthetic
  quadratic; non-uniform grids are rejected.
- **G6** FoM: falls back to |S| at zero width; ranking is descending; every row carries its bracket.
- **G7** the committed config runs and the real sweep artifacts exist.
- **G8** closed form == exact diagonalization to < 1e-8 on all four halides.
- **G9 (THE FINDING)** field ranking is `degenerate` (all FoM tied); strain ranking is `monotone`
  (rank 1 is the window edge — widening the window moves it); **no axis yields an interior
  optimum**; every artifact discloses this; `krylov_dim` is gone from config and artifacts.

## 6. Out of scope

- Real strained-geometry CIF→CASCI clusters (the prerequisite for a non-vacuous certified FoM).
- SOC (Cl chosen first precisely to avoid it); uniaxial-vs-biaxial strain scaling.
- Multi-halide / multi-axis joint screening.
