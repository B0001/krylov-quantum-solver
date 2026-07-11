# Architecture Documentation — krylov-quantum-solver ecosystem

Covers the solver repo plus the four spec'd products: CertChem (certified chemistry
API), ChemCheck (hardware honesty benchmark), ShallowForge (shallow-evolution
compiler), SenseForge (Nb₃X₈ sensor design pipeline).

If it is not written down, it did not happen. This tree is the write-down.

## Layout
```
architecture/
├── adr/            Architectural Decision Records (context → decision → consequences)
│   ├── 0001  Variational floor as a hard invariant        [repo]
│   ├── 0002  ODMD as default estimator                    [repo]
│   ├── 0003  Finite-cluster scope for materials           [repo, SenseForge]
│   ├── 0004  certified vs fast: never blurred             [CertChem]
│   ├── 0005  Library-first, service-second                [CertChem, SenseForge]
│   ├── 0006  Fixed circuits + compilation audit           [ChemCheck]
│   ├── 0007  CX@ε as the only headline metric             [ShallowForge]
│   ├── 0008  Content-hash cache + async jobs              [CertChem]
│   └── 0009  Cloud Run stateless deploy, 0.0.0.0:$PORT    [deployment]
├── c4/             Structural models (Mermaid; renders on GitHub)
│   ├── c1-system-context.md          who uses what, external systems
│   ├── c2-containers-certchem.md     API / Redis / workers / storage
│   ├── c3-components-solver-core.md  the library's internal decomposition
│   └── c3-components-chemcheck.md    the scorer's decomposition
├── interfaces/     Contractual, immutable module boundaries
│   ├── certchem-openapi.yaml            HTTP contract (OpenAPI 3.1)
│   ├── solver-library-contract.md       in-process Python contract + invariants
│   ├── chemcheck-submission.schema.json vendor → scorer boundary
│   ├── chemcheck-scorecard.schema.json  scorer → world boundary
│   └── compiler-manifest.schema.json    compiler → everything provenance
├── data/
│   └── erd.md      Entity-relationship diagram + schema rules with teeth
└── deployment/
    ├── topology.md network diagram, trust rules, failure behavior
    ├── Dockerfile  one image, api|worker entrypoints
    └── main.tf     Terraform: Cloud Run, Redis, GCS, least-privilege SAs
```

## The five load-bearing decisions (read these first)
1. **ADR-0001** — no code path returns a floor-violating energy. Everything honest
   about this system descends from that.
2. **ADR-0004/0005** — the certificate is the product; the library is the boundary;
   HTTP is a shell.
3. **ADR-0006/0007** — benchmark scores and compiler claims are audit-linked through
   one artifact: the provenance manifest.
4. **ADR-0003** — materials claims are cluster-scoped until validated otherwise.
5. **ADR-0008/0009** — determinism → content-hash caching → stateless disposable
   containers.

## Maintenance rules
- A structural change without an ADR is not merged. New ADRs supersede rather than
  edit (mark old ones `Superseded by ADR-NNNN`).
- interfaces/ files are versioned contracts: breaking changes bump versions, never
  mutate in place.
- Diagrams are Mermaid-in-Markdown so they diff, review, and render like code.
