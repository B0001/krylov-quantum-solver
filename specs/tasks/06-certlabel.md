# Task Breakdown 6 — #5 CertLabel: Certified ML-Potential Audit Set
Goal: 1,000 certified labels + the audit-notebook demo. Depends: CertChem-M1.

1. **Fragment/geometry design** — choose H/C/N/O fragments within caps (H₂, H₃⁺, H₄, LiH, CH-type minimal fragments as feasible, N₂ curve, stretched geometries emphasized — multireference points are where ML potentials lie most). Grid: bond scans + angular perturbations.
   ✓ Design doc: N fragments × M geometries ≈ 1,000; every entry cap-checked. (M)
2. **Batch runner** — thin loop over `certified_energy()` writing (geometry, L, U, best, certificate-ref) JSONL; resumable via cache keys.
   ✓ 50-label pilot completes; cost-per-label measured → full-run budget known. (M)
3. **Full generation run** — the 1,000; QC pass: zero floor violations (guaranteed), bracket-width distribution plotted, outliers investigated.
   ✓ Dataset file + stats page; width distribution published, not hidden. (L)
4. **Dataset card + license** — format spec, generation provenance (solver version, config hashes), scope statement (gas-phase, active-space-relative — the ADR-0004 sentence adapted), CC-BY or similar.
   ✓ A stranger can regenerate any label from the card alone. (S)
5. **The audit notebook (the marketing)** — pick one published ML potential covering the fragment space; evaluate on the audit set; count predictions falling outside certified brackets vs the potential's claimed error bars.
   ✓ Honest either way: violations found = headline demo; none found = "potential X passes a certified audit" is also a publishable, friendly result. (M)
6. **Release** — Zenodo/HF dataset + notebook + short writeup.
   ✓ DOI minted; audit reproducible from a fresh clone. (S)
