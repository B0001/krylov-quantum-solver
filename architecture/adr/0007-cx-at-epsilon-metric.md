# ADR-0007: CX@ε is ShallowForge's only headline metric

**Status:** Accepted
**Scope:** ShallowForge; feeds ChemCheck Mode A resource model

## Context
Compiler literature routinely reports gate-count reductions without fixing the
downstream accuracy cost, making results incomparable and often misleading. A
cheaper circuit that ruins the final energy estimate is not an optimization.

## Decision
The headline metric is CX@ε: two-qubit gates per evolution step such that the final
end-to-end ODMD ground-state error on the golden suite stays ≤ ε = 1.6 mHa.
No gate count is ever reported without its ε. Every accuracy-losing transform
(THC rank truncation, Trotter order, randomization) must emit a predicted ε
contribution into the circuit's provenance manifest; total predicted ε must bound
observed ε, or the error model is declared wrong and fixed before proceeding.

## Consequences
- (+) Results are honest, comparable, and directly consumable by ChemCheck's
  headroom model.
- (+) The manifest doubles as ChemCheck's compilation-audit artifact (ADR-0006).
- (−) Every experiment requires a full end-to-end ODMD run, not just a gate count —
  slower iteration, accepted deliberately.
