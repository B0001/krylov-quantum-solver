# SPEC: ChemCheck — hardware chemistry honesty benchmark (M1 + Mode A)

**Slug:** `chemcheck` · **Depends on:** `certchem` (floor guard), `hybrid_quantum_solver`
(`build_molecular_hamiltonian`, `HardwareKrylovSolver.resource_report`), full PRD
`specs/full/spec-hardware-honesty-benchmark.md`, tasks `specs/tasks/02-chemcheck.md`.

## Goal

A neutral, reproducible scorecard answering "can this QPU do real chemistry yet, and if not how
far off?" This spec covers **M1** (frozen tier registry) and the **Mode A** paper scorer
(device spec sheet → PASS/FAIL + fidelity-headroom per tier), plus the Mode B **floor detector**
(the anti-fraud core). The M2 calibration harness against noisy simulation and the M4 launch
report are out of scope here (see Caveats).

## Interface

```python
from chemcheck import (
    TIERS, BENCHMARK_VERSION, build_tier_hamiltonian, canonical_hamiltonian_sha256,
    routing_overhead, expected_total_error, required_two_qubit_error, headroom_factor,
    validate_submission, score_mode_a, render_markdown, mode_b_energy_verdict,
)
```

## Acceptance gates (`tests/test_chemcheck_spec.py`)

1. **Registry loads + hashes stable.** `TIERS` imports with zero solver deps; each T0–T3 frozen
   `hamiltonian_sha256` / `fci_reference_hartree` / `two_qubit_gates_per_trotter_step` matches a
   live recompute. T4 is present and `aspirational == True` with no frozen reference.
2. **Submission validation.** Valid spec sheets pass; malformed ones raise with a JSON-pointer
   path to the offending field. Mode-A submissions need only `device_spec`.
3. **Routing overhead ordering.** `overhead(all_to_all) == 1.0`; `heavy_hex > grid > all_to_all`.
4. **Error budget is a pure function** matching hand-computed depolarizing cases.
5. **Headroom monotonicity.** Better (smaller) two-qubit error → smaller headroom;
   `headroom == 1.0` exactly at the required threshold; `headroom <= 1` ⇔ PASS.
6. **Floor detector.** Known-bad energies (hundreds of Ha below FCI) → 100% `UNPHYSICAL`; exact
   golden FCI energies → zero false positives (`PASS`).
7. **Scorecard emitter.** `score_mode_a` output validates against
   `architecture/interfaces/chemcheck-scorecard.schema.json`; `render_markdown` carries the
   `classically_simulable` disclaimer on T0–T2.

## Out of scope / honest caveats

- **Depolarizing-only, single-Trotter-step budget.** The v1 error model counts one reference
  Trotter step of two-qubit gates × a published routing multiplier, and uses the "≤ 1 expected
  two-qubit error per circuit" coherence heuristic as the PASS threshold. It is deliberately
  crude and **uncalibrated** — M2 (noisy-sim sweep) is what would validate the pass/fail
  transition. Headroom factors are order-of-magnitude, not certified.
- Routing multipliers are literature lookup values, not compiled-circuit measurements.
- Mode B is only the floor detector here; full ODMD energy extraction from submitted `runs`,
  the reproducibility check, and the compilation audit are M3.
