# Deployment Topology

Only CertChem is a hosted service. ChemCheck, ShallowForge, and SenseForge are
batch/CLI workloads that run locally or in CI and publish artifacts — they need
no ingress at all, which is itself a security decision.

```mermaid
flowchart LR
    subgraph internet["Internet"]
        user["API clients"]
    end

    subgraph gcp["GCP project: certchem-prod (Terraform)"]
        subgraph ingress["Ingress"]
            lb["Cloud Run HTTPS endpoint<br/>TLS termination, managed cert"]
        end
        api["certchem-api<br/>Cloud Run service<br/>min=0, max=10 instances<br/>concurrency=20<br/>listens 0.0.0.0:$PORT"]
        redis["Memorystore Redis<br/>private VPC IP only<br/>queue + cache index + rate counters"]
        worker["certchem-worker<br/>Cloud Run job, 2 vCPU / 4 GiB<br/>max parallelism 5<br/>no public ingress"]
        gcs["GCS bucket certchem-results<br/>uniform access, versioned"]
        sa1["SA: api-sa<br/>redis rw, gcs read"]
        sa2["SA: worker-sa<br/>redis rw, gcs write"]
    end

    ci["GitHub Actions CI<br/>golden suite gate +<br/>ChemCheck/ShallowForge batch runs"]

    user -->|"HTTPS 443, X-Api-Key"| lb --> api
    api <-->|"6379, VPC connector"| redis
    worker <-->|"6379, VPC connector"| redis
    worker -->|"write blobs"| gcs
    api -->|"read blobs"| gcs
    ci -->|"deploy on green<br/>(containment failure blocks)"| api
    ci -->|"deploy image"| worker
```

## Traffic and trust rules
| Path | Allowed | Denied |
|---|---|---|
| Internet → api | HTTPS only, valid API key, per-key rate limit | HTTP, anonymous (except GET /v1/limits) |
| Internet → worker / redis / gcs | — | everything (no public surface) |
| worker egress | Redis + GCS via VPC | general internet (chemistry needs no egress) |
| api → worker | none directly | direct RPC (all coordination via Redis) |

## Scaling and failure behavior
- api scales to zero when idle (ADR-0009); cold start accepted (success criterion
  is 60 s end-to-end, not milliseconds).
- Long jobs go through the queue; Cloud Run request timeouts therefore never
  truncate a computation — a timeout kills at most a poll.
- Redis loss = queue/cache loss only: results in GCS survive, cache repopulates
  by recomputation (determinism, ADR-0008). Nothing irreplaceable lives in Redis.
- Deploy gate: CI runs the golden regression suite; one bracket-containment
  failure blocks the deploy (CertChem spec §7.1).
