# ADR-0006: ChemCheck circuits are fixed, published, and audit-verified

**Status:** Accepted
**Scope:** ChemCheck

## Context
Benchmarks get gamed. The two cheapest attacks: (a) "compile" the circuit into a
different, easier problem; (b) launder classical computation through unbounded
post-processing so the QPU is decorative.

## Decision
1. Tier Hamiltonians, reference circuits, seeds, and scoring code are published and
   pinned by version hash (benchmark versions: ChemCheck-YYYY.N).
2. Vendor compilation is allowed, but the compiled circuit must be submitted and is
   re-verified against the tier Hamiltonian (unitary/spectral check on simulable
   sizes; structural audit above that). ShallowForge's provenance manifests
   (ADR-0007 scope) satisfy this requirement natively.
3. Error mitigation is permitted but declared; raw and mitigated scores are both
   reported. Classical post-processing budgets are capped and disclosed.
4. Tiers T0–T2 are labeled classically simulable: passing them certifies stack
   correctness, not quantum advantage.

## Consequences
- (+) Scores are comparable within a version and resistant to the cheapest attacks.
- (−) Verification machinery is real work and must ship before any leaderboard.
- (−) Vendors may dispute fairness; mitigation is full methodology publication and
  a PR process, at the cost of governance overhead.
