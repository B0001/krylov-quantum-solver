# ADR-0009: Stateless containers on Cloud Run; bind 0.0.0.0:$PORT

**Status:** Accepted
**Scope:** CertChem deployment; pattern reused for any hosted component

## Context
Prior deployment experience (the Dash bible app) surfaced the classic trap: dev
servers default to binding 127.0.0.1, which is unreachable behind a container
platform's ingress. Candidate platforms evaluated: Google Cloud Run, Render,
Railway. Cloud Run chosen for scale-to-zero pricing, native job/queue story
(Cloud Run Jobs), and managed Redis (Memorystore) adjacency; the container is
platform-portable regardless.

## Decision
- All services run as stateless containers listening on `0.0.0.0:$PORT` (port read
  from the environment; never hardcoded, never loopback).
- State lives only in managed services: Redis (queue + cache index), object storage
  (result blobs, benchmark bundles). Containers are disposable.
- Infrastructure is declared in Terraform (see deployment/); no console-clicked
  resources.
- Ingress: HTTPS only, API-key auth at the application layer, per-key rate limits
  from day one. Workers have no public ingress; egress restricted to Redis/storage.

## Consequences
- (+) Scale-to-zero keeps idle cost near nil; horizontal scaling is configuration.
- (+) The container runs unchanged on Render/Railway if platform migration is needed.
- (−) Cold starts affect first-request latency; acceptable for a compute-bound API
  (success criterion is 60 s end-to-end, not 60 ms).
- (−) Long jobs must respect platform timeouts → queue tier (ADR-0008) is mandatory,
  not optional.
