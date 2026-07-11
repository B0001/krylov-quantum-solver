# Task Breakdown 14 — #19 Bracket-Aware Screening Loop
Goal: interval-dominance screening with provably zero false eliminations, demonstrated on a 50-candidate toy space. Depends: CertChem-M1.

1. **Formalize the pruning rule** — candidate eliminated iff its certified interval strictly excludes the target region. One-page proof note: under bracket containment (ADR-0001 guarantee), elimination is sound.
   ✓ Note reviewed; edge cases (interval touching target boundary) resolved conservatively. (S)
2. **Toy candidate space** — 50 cap-compliant variations (e.g. H₄ geometries / small fragments parameterized by geometry) with exhaustive FCI ground truth computed once.
   ✓ Ground-truth table; target property region chosen so ~10/50 are true hits. (M)
3. **Loop v1: prune-only** — evaluate all certified brackets at pilot precision; eliminate by dominance; refine survivors (tighter target_width via #15's planner if built, else fixed schedule).
   ✓ Zero false eliminations vs ground truth (the invariant claim); total evaluations counted. (M)
4. **Loop v2: acquisition** — order refinement by overlap-with-target heuristic (interval analog of expected improvement).
   ✓ Fewer than half the evaluation cost of exhaustive-at-full-precision, same final hit set. (M)
5. **Ablation + writeup** — v1 vs v2 vs point-estimate BO baseline (which WILL make false eliminations — show it doing so): the demonstration that motivates the method.
   ✓ Plot: false-elimination count (baseline > 0, ours = 0) + cost comparison. Publishable core. (M)
6. **Generalize the oracle interface** — abstract `BoundedOracle` protocol so any interval source (conformal ML models) plugs in — the transfer path past the 16-orbital cap.
   ✓ Second oracle (synthetic noisy-ML mock with conformal intervals) runs through unchanged. (S)
