# ADR-0002: ODMD/Krylov as the default eigenvalue estimator

**Status:** Accepted
**Scope:** krylov-quantum-solver core; downstream consequence for ShallowForge

## Context
Candidate estimators: VQE (variational, barren-plateau and optimizer-fragility risks),
textbook QPE (deep circuits, brittle to noise), and ODMD/quantum-Krylov methods
(spectral estimation from noisy time-series data, documented noise robustness).

## Decision
ODMD over a Krylov time-evolution signal is the default and reference estimator.
All accuracy claims, certificates, and benchmark scores are defined with respect to
this pipeline.

## Consequences
- (+) Robustness to stochastic error is inherited by every downstream product.
- (+) Opens ShallowForge's most original research angle: cheap randomized Trotter
  formulas may be tolerable here where stricter pipelines fail (compiler spec, R2).
- (−) Requires real-time evolution circuits, whose two-qubit gate depth (~6,500 CX/step
  for N₂ CAS(6,6)) is the current hardware wall. Depth reduction is therefore a
  first-class project (ShallowForge), not an optimization afterthought.
- (−) Signal-processing hyperparameters (Krylov dimension, filtering) must be recorded
  in every certificate for reproducibility (see ADR-0004).
