# ADR-0005: Library-first; the HTTP service is a thin shell

**Status:** Accepted
**Scope:** CertChem; SenseForge

## Context
Two consumers need the same capability: external users (HTTP) and the SenseForge
pipeline (in-process, thousands of calls per sweep). Building logic into the web
layer would force SenseForge through HTTP or duplicate the logic.

## Decision
All science lives in a Python library entry point:
`certified_energy(molecule, basis, cas, mode) -> Bracket` (plus `certified_gap`,
`certified_reaction`). FastAPI translates HTTP ⇄ library calls and owns only:
validation of caps, auth, caching, queuing, serialization. The library has zero
knowledge of HTTP.

## Consequences
- (+) SenseForge imports the library directly; sweep throughput is not bottlenecked
  by a network hop.
- (+) The golden regression suite tests the library once; the API layer needs only
  contract tests.
- (−) Caps must be enforced twice conceptually (library raises `CapExceededError`;
  API pre-validates to fail fast) — kept consistent by sharing one limits module.
