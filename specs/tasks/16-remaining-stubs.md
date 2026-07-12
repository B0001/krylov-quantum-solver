# Task Breakdown 16 — Remaining Specs: Seed Tasks Only
The four not yet broken down, each reduced to its single blocking first task —
because each is gated on something above, a full breakdown now would be speculation.

## #2 Virtual Spectrometer — gated on: ODMD spectroscopy modules + certified_gaps hardened
**Seed task:** the blind-exclusion pilot. 5 candidate structures, 1 synthetic
"measured" spectrum, tool must exclude the 4 wrong ones with certificates.
✓ Pass/fail is unambiguous; a pass justifies the full breakdown. (L)

## #3 Strain Camera — gated on: SenseForge task 8 (size-converged calibration curve)
**Seed task:** interval inversion module on the Nb₃Cl₈ calibration curve + one
synthetic strain-field image with per-pixel certified bounds.
✓ Inversion round-trips: synth strain → gap → recovered strain ± bound contains truth. (M)

## #9 Reaction-Network Pruner — gated on: #17 task 1 (charged/open-shell) + certified_thermochem at network scale
**Seed task:** 10-species toy network; prune with path-propagated intervals; report
killed edges (with certificates) AND edges surviving only because brackets are too
loose — the honesty readout that decides whether real networks are in reach.
✓ Report written; width-blowup vs path depth quantified. (M)

## #12 Auto-CAS Selector — gated on: CertChem-M1 + budget for a research arc
**Seed task:** the stabilization-criterion study — on the golden suite, grow CAS
stepwise and test whether bracket convergence reliably tracks true CAS adequacy
(FCI references exist precisely here, so ground truth is available).
✓ Correlation quantified; go/no-go for the selector written at the bottom. (L)

---

# Cross-portfolio dependency spine (read this before starting anything)
```
CertChem-M1 (breakdown 01)
 ├─→ #5 CertLabel (06)      ├─→ #15 Planner (10)     ├─→ #18 Registry (13)
 ├─→ #19 Screening (14)     ├─→ #17 Redox (15)       ├─→ #12 Auto-CAS (seed)
 └─→ SenseForge (04) ──→ #3 Strain Camera (seed)
ChemCheck (02) ──→ #8 Co-Design (11) ; ──→ #16 Noise Thermometer
ShallowForge (03) ──→ feeds 02 + 11
Mikra track (independent): #13 Leyning (08), #14 Corpus (09)
Anytime, no deps: #20 Invariant plugin (05), #10 Playground (07 — needs only one precomputed curve)
```
Two-person-weeks of highest leverage, if choosing today: breakdown 01 end-to-end,
then 05 (plugin) and 07 (playground) as shipping-morale wins while 02 calibration runs.
