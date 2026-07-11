# ADR-0001: Variational floor as a hard, non-bypassable invariant

**Status:** Accepted (retroactively documents the rebuild)
**Scope:** krylov-quantum-solver core; inherited by CertChem, ChemCheck

## Context
The original codebase returned ground-state energies hundreds of Hartree below the
theoretical minimum — results that looked numerically confident and were physically
impossible. The failure was silent: nothing in the pipeline knew what "impossible" meant.

## Decision
Every energy estimate produced anywhere in the system is checked against a variational
lower bound before it is returned. A violation raises `FloorViolationError` and aborts
the computation. There is no flag, config option, or code path that returns a
floor-violating number to a caller.

## Consequences
- (+) The old failure mode is structurally impossible, not merely unlikely.
- (+) The same check becomes ChemCheck's public "UNPHYSICAL" fraud detector for free.
- (+) CertChem's certificate contract (ADR-0004) rests on this invariant.
- (−) Small runtime overhead per estimate; requires a floor to be computable for every
  supported system class — this constrains which systems can enter the supported set.
- (−) Legitimate-but-buggy runs fail loudly instead of returning "close" numbers;
  this is intended, but it means CI noise until upstream bugs are fixed.
