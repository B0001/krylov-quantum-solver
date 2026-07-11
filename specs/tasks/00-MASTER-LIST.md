# Master List — All 20 Specs

## Full specs (files delivered earlier)
| # | Name | One-line claim | Status |
|---|------|----------------|--------|
| 1 | **SenseForge** | Nb₃X₈ strain/field sensor operating points via certified gap screening | Full spec |
| 2 | **CertChem** | Chemistry API returning provable two-sided bounds | Full spec |
| 3 | **ChemCheck** | Hardware honesty benchmark: "can this QPU do chemistry yet?" | Full spec |
| 4 | **ShallowForge** | Circuit compiler cutting CX@1.6mHa by 5–10× | Full spec |

## Portfolio I mini-specs (5–12)
| # | Name | One-line claim |
|---|------|----------------|
| 5 | CertLabel | ML-potential training/audit datasets with certified label bounds |
| 6 | Virtual Spectrometer | Structure ID by certified-interval spectral exclusion |
| 7 | Strain Camera | Invert gap-vs-strain calibration into strain imaging |
| 8 | Co-Design Advisor | Target molecule → required device spec (ChemCheck inverted) |
| 9 | Reaction-Network Pruner | Certified pruning of kinetics networks via interval propagation |
| 10 | Multireference Playground | Interactive "watch HF fail" teaching app (precomputed, Dash) |
| 11 | ODMD-Anywhere | ODMD as a standalone classical signal-processing library |
| 12 | Auto-CAS Selector | Auditable automatic active-space selection via bracket stabilization |

## Portfolio II mini-specs (13–20)
| # | Name | One-line claim |
|---|------|----------------|
| 13 | Leyning Coach | Word-level cantillation feedback on Mikra Sync alignment |
| 14 | Open Aligned Tanakh Corpus | "LibriSpeech of Biblical Hebrew" — aligned audio dataset |
| 15 | Shots-to-Certainty Planner | Measurement budget → guaranteed bracket width forecaster |
| 16 | Noise Thermometer | Golden-suite circuits as workload-relevant device-noise probes |
| 17 | Redox-Window Bounder | Certified gas-phase IP/EA pre-screening tables |
| 18 | Certificate Registry | Public verify-by-recompute ledger for computed numbers |
| 19 | Bracket-Aware Screening Loop | Provably-safe candidate elimination via interval dominance |
| 20 | Invariant-Test Framework | Physics-invariants-as-CI-gates pytest plugin |

## Task-breakdown build order (this package)
Foundation first (everything imports it), then leverage, then quick wins, then the rest:
CertChem-M1 → ChemCheck → ShallowForge → SenseForge → #20 → #1 → #10 → #13 → #14 → #15 → …
Each breakdown = numbered tasks with deliverable, acceptance criterion, and rough effort.
Effort key: S = ≤half a day, M = 1–3 days, L = ~1 week.

## Breakdown index (delivered in this package)
01 CertChem-M1 · 02 ChemCheck · 03 ShallowForge · 04 SenseForge · 05 #20 Invariant plugin ·
06 #5 CertLabel · 07 #10 Playground · 08 #13 Leyning Coach · 09 #14 Tanakh Corpus ·
10 #15 Shots Planner · 11 #8 Co-Design · 12 #11 ODMD-Anywhere · 13 #18 Registry ·
14 #19 Screening Loop · 15 #17 Redox Bounder · 16 seed tasks for #2/#3/#9/#12 + dependency spine
