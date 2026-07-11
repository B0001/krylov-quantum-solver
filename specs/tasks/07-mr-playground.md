# Task Breakdown 7 — #10 Multireference Playground (teaching app)
Goal: slider-driven N₂ dissociation demo, precomputed, deployed on Cloud Run. Reuses your Dash + 0.0.0.0:$PORT experience directly.

1. **Precompute the payload** — N₂ CAS(6,6) dissociation curve, ~40 bond lengths: HF energy, certified bracket, FCI reference per point → one static `n2.json`.
   ✓ JSON < 100 KB; HF−FCI gap reaches ~0.5 Ha at stretch — the money plot exists. (M)
2. **Single-page app** — Dash (or plain HTML+Chart.js — lighter): bond-length slider, three traces (HF, bracket band, FCI), live annotation "HF error here: X mHa". No backend chemistry — serve static JSON.
   ✓ Loads < 1 s; works on a phone (your students will use phones). (M)
3. **The pedagogy layer** — 3 short captions triggered by slider zones: equilibrium ("HF fine here — why?"), intermediate, dissociation ("this is strong correlation"). One "what's a bracket?" expandable.
   ✓ A chemistry undergrad who's had one quantum course can narrate the plot back to you. Test on one real human. (S)
4. **Deploy** — containerize, bind 0.0.0.0:$PORT, Cloud Run scale-to-zero.
   ✓ Public URL; cold start acceptable; costs ≈ $0 idle. (S)
5. **Second molecule toggle (stretch)** — H₄ square geometry for a different flavor of degeneracy.
   ✓ Same JSON pattern; toggle, not a rewrite. (M)
6. **Distribution move (the actual killer risk)** — pair with one concrete outreach: a chem-ed mailing list post, a workshop slot, or a short J. Chem. Educ.-style note.
   ✓ One committed distribution action scheduled before you call it shipped. (S)
