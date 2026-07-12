# ADR-0003: Materials claims are scoped to finite-cluster models

**Status:** Accepted
**Scope:** Nb₃X₈ modules; SenseForge

## Context
The CIF ingestion path builds finite (dimer-scale) clusters, not periodic solids.
Periodic embedding is a large, separate research effort. Claims that quietly conflate
cluster predictions with bulk-material properties are a known failure pattern in the
field and would repeat the overclaiming this repo was rebuilt to escape.

## Decision
1. Every artifact derived from the CIF path carries an automatic header:
   "cluster-model prediction; not validated for the periodic solid."
2. SenseForge results require two gates before external claims: cluster-size
   convergence of the candidate ranking, and cross-method agreement (DMRG/AFQMC)
   within combined error bars.
3. Periodic support, if ever added, arrives as a new system class with its own
   validation suite — not as a widening of the current one.

## Consequences
- (+) Publishable claims are defensible; the sensor "invention" is phrased as a
  design prediction with stated model scope.
- (−) Slower path to headline materials results; some reviewers/readers will want
  bulk numbers we deliberately refuse to provide.
- (−) Nb₃I₈ additionally requires spin–orbit coupling; until SOC exists, iodide
  outputs are labeled indicative-only (SenseForge spec §8).
