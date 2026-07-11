# C4 Level 2 — Containers: CertChem deployment

```mermaid
flowchart TB
    client["API client<br/>(HTTPS + API key)"]

    subgraph gcp["Cloud Run project (Terraform-managed)"]
        api["API Service<br/>[Cloud Run service: FastAPI]<br/>validation, auth, rate limits,<br/>cache lookup, job dispatch<br/>binds 0.0.0.0:$PORT"]
        redis["Redis<br/>[Memorystore]<br/>job queue + cache index<br/>+ rate-limit counters"]
        workers["Worker Pool<br/>[Cloud Run jobs]<br/>one job = one process,<br/>BLAS threads pinned<br/>runs Solver Core"]
        gcs["Object storage<br/>[GCS bucket]<br/>result blobs, certificates,<br/>golden-suite references"]
    end

    client -->|"POST /v1/energy|reaction|gap"| api
    client -->|"GET /v1/jobs/{id}"| api
    api <-->|"enqueue / status"| redis
    api <-->|"cache key lookup"| redis
    workers -->|"BRPOP jobs"| redis
    workers -->|"write results"| gcs
    api -->|"read results"| gcs
```

Container responsibilities
| Container | Owns | Explicitly does not own |
|---|---|---|
| API Service | HTTP contract, caps, auth, cache/queue orchestration | any chemistry |
| Worker Pool | executing `certified_energy()` et al. | HTTP, auth |
| Redis | queue, cache index, rate counters | result payloads (blobs go to GCS) |
| GCS | durable results + certificates | mutable state |

Sync path: cache hit or job estimated <30 s → response inline.
Async path: 202 + job_id → worker → GET /v1/jobs/{id} (ADR-0008).
