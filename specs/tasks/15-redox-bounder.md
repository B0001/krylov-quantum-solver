# Task Breakdown 15 — #17 Redox-Window Bounder
Goal: certified gas-phase IP/EA table for 20 small fragments, cross-checked vs experiment. Depends: CertChem-M1 (charged-species support is the new work).

1. **Charged/open-shell support audit** — verify the pipeline handles cation/anion + doublet states within caps (PySCF side is fine; check ODMD + floor + Temple path for open shells).
   ✓ H₂⁺ and Li end-to-end with certified brackets; any open-shell gap in the bounds machinery documented before proceeding. (M)
2. **Fragment list** — 20 species with experimental gas-phase IP/EA reference data (NIST WebBook coverage confirmed per species) AND cap compliance.
   ✓ Table: species, CAS choice, reference value + experimental uncertainty. (S)
3. **`certified_ip_ea` wrapper** — energy difference between charge states via `certified_reaction()`; interval propagation as-is.
   ✓ Unit test on H₂ → H₂⁺ against exact values. (S)
4. **Production run + cross-check** — all 20; three-way comparison: bracket vs experiment vs (where available) high-level literature theory. KEY nuance: bracket is active-space-relative — agreement analysis must separate "bracket wrong" (never acceptable) from "active space too small to match experiment" (expected, informative).
   ✓ Zero containment violations vs FCI-in-active-space checks; experiment-gap analysis written per species. (L)
5. **The solvation disclaimer, structurally** — every output row carries "gas-phase intrinsic value; solvation shifts of O(eV) apply" — automated header, ADR-0003 pattern.
   ✓ Grep test. (S)
6. **Release table + memo** — where intrinsic pre-screening would/wouldn't have caught known electrolyte failures (literature mini-review, 3 cases).
   ✓ Memo makes the use-case concrete or honestly reports it doesn't. (M)
