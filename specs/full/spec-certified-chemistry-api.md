# Spec: Certified Chemistry API

**Working name:** CertChem (rename freely)
**Status:** Draft v0.1
**Depends on:** krylov-quantum-solver (`temple_bounds.py`, `certified_gaps.py`, `certified_thermochem.py`, variational floor, core eigensolver)

---

## 1. Goal

A web service that returns quantum-chemistry quantities with **rigorous two-sided bounds** instead of a bare number: "the ground-state energy is between L and U, provably." The unique selling point is the certificate, not the solver — most tools return a point estimate with no guarantee; this one returns a bracket you can act on.

## 2. Target users

- ML-potential developers who need trusted labels to validate cheaper models against
- Battery/catalyst screeners who want a "no surprises" pre-filter on small critical fragments
- Educators/researchers who want FCI-quality answers on small systems without setting anything up

## 3. Scope limits (v1, hard caps)

- Active space ≤ 16 spin orbitals (i.e., within the sizes the repo's test suite validates)
- Molecules specified by XYZ (SMILES accepted, converted via RDKit, flagged as geometry-not-optimized unless user opts into a cheap pre-optimization)
- Non-periodic, gas-phase, no solvation
- Per-job wall-clock timeout (default 10 min; queued tier up to 2 h)

Caps are features: everything inside the fence is *validated*; everything outside returns HTTP 422 with a clear reason, never a degraded silent answer.

## 4. API design

### `POST /v1/energy`
```json
{
  "molecule": {"xyz": "...", "charge": 0, "multiplicity": 1},
  "basis": "sto-3g | 6-31g",
  "active_space": {"electrons": 6, "orbitals": 6} ,
  "mode": "fast | certified"
}
```
Response:
```json
{
  "best_estimate_hartree": -109.1234,
  "lower_bound_hartree": -109.1241,
  "upper_bound_hartree": -109.1229,
  "certificate": {
    "method": "temple_bound + variational_floor",
    "floor_check": "pass",
    "krylov_dim": 24,
    "convergence": "converged"
  },
  "job_id": "…", "wall_time_s": 41.2
}
```

### `POST /v1/reaction`
List of species with stoichiometric coefficients → reaction energy with **propagated brackets** (`certified_thermochem.py`). Bracket arithmetic: ΔE ∈ [Σν·L or U chosen worst-case]. Interval widths add — document this loudly so users size their active spaces accordingly.

### `POST /v1/gap`
Excitation / spin gap with certified bracket (`certified_gaps.py`).

### `GET /v1/jobs/{id}`
Async status + result retrieval for queued jobs.

### `GET /v1/limits`
Machine-readable statement of current caps and validated system classes.

## 5. Architecture

- **FastAPI** app + worker pool (Celery or RQ + Redis) for anything > ~30 s
- Solver runs in worker processes; one job = one process (NumPy/BLAS thread caps set explicitly)
- **Result cache** keyed on hash(molecule + basis + active space + mode) — chemistry is deterministic, cache aggressively
- Containerized; deploy on Cloud Run / Render / Railway. Same lesson as the bible app: bind `0.0.0.0:$PORT`, not `127.0.0.1`, and read the port from the environment
- Rate limiting per API key from day one (compute is the cost center)

## 6. The certificate contract (the actual product)

Every `certified`-mode response promises:

1. **Containment:** the true FCI energy for the requested active space lies in [L, U].
2. **Floor check:** the estimate was verified above the variational floor; any violation aborts the job with an error rather than returning a number (the exact failure mode the old codebase had — hundreds of Hartree below the minimum — must be structurally impossible to return).
3. **Traceability:** certificate lists method, Krylov dimension, convergence status, and solver version; any response is reproducible from its certificate.

`fast` mode returns a point estimate with no bracket and is labeled `"certified": false`. Never blur the two.

## 7. Validation plan

1. **Golden regression suite as CI gate:** H₂, H₄ chain, LiH, N₂ CAS — every deploy must reproduce FCI references within documented tolerance and every bracket must contain the reference. One containment failure blocks deploy.
2. **Fuzz set:** ~200 random small molecules/geometries within caps; check floor violations = 0, bracket sanity (L ≤ estimate ≤ U), no NaNs, timeout behavior clean.
3. **Interval-arithmetic tests** for `/v1/reaction`: hand-checked worst-case propagation on 3 known reactions.
4. **Load test:** 50 concurrent `fast` jobs on the smallest instance tier; document the real throughput before writing any pricing/limits page.

## 8. Risks

- **Loose brackets on hard systems.** Temple-type bounds can be wide when the gap is small. Mitigation: report the bracket honestly and expose `bracket_width` so users can filter; never tighten by heuristics.
- **Auto active-space selection is a research problem.** v1 requires the user to specify it; an `auto` mode is explicitly out of scope until there's a validated selector.
- **Compute cost.** Cache + caps + rate limits from day one; measure before offering anything free-tier.
- **Overclaim risk.** All docs state: certified with respect to *the requested active space and basis*, not the exact non-relativistic limit. This sentence appears in the API response docs verbatim.

## 9. Milestones

1. **M1 — Solver-as-library (1 week):** clean Python entry point `certified_energy(mol, basis, cas) -> Bracket` with tests; no web layer yet.
2. **M2 — FastAPI wrapper + cache (1 week):** `/v1/energy` sync-only, golden suite in CI.
3. **M3 — Async queue + `/v1/reaction`, `/v1/gap` (1–2 weeks).**
4. **M4 — Deploy + fuzz + load test (1 week):** public endpoint behind API keys.
5. **M5 — Docs site + example notebooks:** "validate your ML potential against certified labels" as the flagship tutorial.

## 10. Success criteria

- Zero bracket-containment failures across golden + fuzz suites.
- Cold-start-to-answer under 60 s for an H₂/LiH-class `certified` job on the deployed tier.
- One external user (or one of your own projects, e.g. feeding the Nb₃X₈ sensor pipeline) consuming the API in anger.
