# ADR-0008: Deterministic content-hash caching; async for anything slow

**Status:** Accepted
**Scope:** CertChem service layer

## Context
Quantum-chemistry results are deterministic functions of
(molecule geometry, charge, multiplicity, basis, active space, mode, solver version).
Compute is the dominant cost; identical requests are common (regression suites,
sweeps, retries).

## Decision
- Cache key = SHA-256 over the canonicalized request tuple **including solver
  version**; cache hit returns the stored result byte-identically.
- Requests estimated > ~30 s run through a Redis-backed queue (worker pool);
  the API returns 202 + job_id, results retrieved via GET /v1/jobs/{id}.
- One job = one worker process, with BLAS/OMP thread counts pinned explicitly.

## Consequences
- (+) Sweeps and CI re-runs are nearly free; retries are safe (idempotent).
- (+) Including solver version in the key means upgrades never serve stale physics.
- (−) Cache is useless across solver releases by design — accepted.
- (−) Two response flows (sync/async) in the OpenAPI contract.
