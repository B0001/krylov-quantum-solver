# SPEC: SenseForge — certified Nb₃X₈ strain/field sensor screening

**Slug:** `senseforge` · **Status: BOUNDARY RECORDED — premise unsatisfiable on accessible
models; not built.** Full PRD `specs/full/spec-nb3x8-sensor-designer.md`, tasks
`specs/tasks/04-senseforge.md`.

## Intended goal

Rank Nb₃X₈ strain/field sensor operating points by a **certified** figure of merit
`FoM = |sensitivity| / certified-bracket-width`, using `certified_gap` on strained/field-perturbed
cluster Hamiltonians. Only operating points whose certified gap bracket is both sensitive and
tight would be promoted.

## Why it is not built (the finding)

The certified-FoM premise has **no useful regime** on the Nb₃X₈ models the repo can actually reach
(measured 2026-07-10 with `certified_gaps.gap_bracket` on `nb3x8_gaps` clusters):

| System | Certified spin-gap bracket | Verdict |
|---|---|---|
| Nb₃Cl₈ dimer (validated model, 4 qubits) | width ≈ 0–1e-3 meV | **trivial** — the dimer is exactly solvable, so the bracket carries no information; `FoM = |S|/width → ∞`. |
| 4-site half-filled cluster (8 qubits) | width ≈ 2400 meV, gap ≈ 620 meV | **vacuously loose** — the real-time Krylov space from HF cannot resolve the spin excitation, so the Temple/Weinstein lower certificate is ≫ the gap; `FoM ≈ 0`. |

The *sensitivity itself is real* (dΔ/d|t| ≈ 2.5 meV per % on Cl), so a **point-estimate** sensor
screener is viable — but the **certified** discriminator that is SenseForge's whole reason to
exist does not have a tight-but-nontrivial operating regime on these systems. Certification would
only earn its keep on CIF→CASCI clusters large enough that FCI is out of reach yet well-enough
resolved that the bracket stays tight — a regime this model family does not provide at accessible
Krylov depth, and which requires net-new strained-geometry generation (the repo has only a
`|t|`-proxy strain, not real lattice strain — see recon).

## Disposition

Paused in favor of the portfolio mini-specs (#5–19) that build cleanly on the shipped
`certchem`/`chemcheck` foundations. Revisit only with (a) a strained-geometry CIF→CASCI path and
(b) evidence that the certified gap bracket is tight-but-nontrivial on the resulting clusters.
This mirrors the repo's SDD practice: an unsatisfiable gate is the finding, not a failure
(cf. `SPEC_hchain_tdl.md`).
