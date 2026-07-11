# ADR-0004: "certified" and "fast" are separate, never-blurred contracts

**Status:** Accepted
**Scope:** CertChem API; solver library

## Context
The product's value is the guarantee, not the number. Any ambiguity about whether a
response is guaranteed destroys that value: one user burned by a "certified-ish"
answer ends trust in the whole service.

## Decision
- `mode=certified` responses always contain a two-sided bracket [L, U], a floor-check
  result, and a machine-readable certificate (method, Krylov dim, convergence,
  solver version). The containment promise: FCI for the requested active space and
  basis lies in [L, U].
- `mode=fast` responses carry `"certified": false`, no bracket, and share no wording
  with certified responses.
- The guarantee is always stated relative to the requested (active space, basis) —
  never the exact non-relativistic limit. This sentence appears verbatim in API docs.
- Bounds are never tightened heuristically. Wide-but-true beats tight-but-hopeful.

## Consequences
- (+) The certificate is auditable and reproducible; disputes are resolvable from
  the certificate alone.
- (−) Brackets on small-gap systems can be embarrassingly wide (Temple bounds
  degrade); we expose `bracket_width` and let users filter rather than hide it.
- (−) Two response schemas to maintain and test.
