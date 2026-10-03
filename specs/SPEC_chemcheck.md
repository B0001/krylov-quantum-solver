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

1. **Registry loads + frozen values reproduce.** `TIERS` imports with zero solver deps; each T0–T3
   frozen `fci_reference_hartree` and Pauli-term count matches a live recompute **on every
   platform**; `two_qubit_gates_per_trotter_step` matches on every platform for T0–T2 and on Linux
   for T3; the frozen `hamiltonian_sha256` matches **on Linux only** (platform note below). T4 is
   present and `aspirational == True` with no frozen reference.
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

## Platform note — the frozen hash is Linux-reference-only (chem-7rb, 2026-10-02)

`canonical_hamiltonian_sha256` was documented as "stable across machines". It is not. On macOS 27
(Apple M3, scipy 1.15.3 on Accelerate, qiskit-nature 0.8.0, pyscf 2.13.1)
`test_frozen_tier_values_match_live_recompute` failed for T1 and T3 on `hamiltonian_sha256` (the
first assertion, so the later ones never ran). Re-measured assertion by assertion: FCI energy (1e-6
Ha) and Pauli-term count match for all four tiers, the T0 and T2 hashes match, and the CX count
matches for T0–T2 but **not T3 (live 6530, frozen 6528; identical for PYTHONHASHSEED 0, 1, 2)** (the
hand-off that scoped this work said CX passes everywhere; the original single test asserted the hash
first, so a failing run never reached the CX assertion):

| tier | live (macOS) | frozen (Linux) | what reproduces it |
|---|---|---|---|
| T1 H₄ chain | `fff8a314e967…` | `2fe671eae4ae…` | an MO **sign flip**: 4 of the 16 sign patterns reproduce the frozen hash exactly (e.g. flipping the 0-based active MO 3) |
| T3 N₂ CAS(6,6) | `ccde9d59df68…` | `382de579ca32…` | **none** of the 64 sign patterns; the active space holds two degenerate π pairs (SCF ε = −0.57139 ×2 at 0-based MO 4–5, 0.28019 ×2 at MO 7–8) |

The hash is a bitwise hash of SCF-derived integrals in whatever MO gauge the platform's eigensolver
returned: eigenvector signs (T1) and, by inference, rotations inside degenerate orbital pairs (T3 —
inferred from the degeneracy, the absence of any sign-gauge match, and the CX difference, but not
verified; it would need the Linux orbitals). The CX evidence: `build_trotter_step` orders terms by
|coefficient| (`trotter_krylov.canonical_term_order`), which a sign gauge cannot change — T1's CX
matches (2170) — while T3's differs, so T3's |coefficients| differ; near-tied |coefficients| breaking
differently at the 1e-16 level would give the same symptom and cannot be separated here. Flipping the sign of spatial MO k is conjugation of the Jordan–Wigner Hamiltonian by
Z_{k↑}Z_{k↓}, which flips the sign of each Pauli term with an odd number of X/Y on qubits k and k+n.
(Source: throwaway probe on this Mac that enumerated all sign patterns of the live Hamiltonian
against the frozen hashes; not committed.)

**Decision.** The hash equality is split into its own test, `test_frozen_tier_hash_matches_live_recompute`,
`skipif(sys.platform != "linux")` with these numbers in the reason; the T3 CX equality is split into
`test_frozen_tier_cx_count_matches_live_recompute[T3]`, skipped off Linux with its own numbers; FCI
energy, Pauli-term count and the T0–T2 CX counts stay on every platform. The
`canonical_hamiltonian_sha256` docstring is corrected. No frozen value changed.
**Not verified:** Linux (the skipped test is the unchanged original assertion).

**Follow-up (out of scope).** A gauge-invariant hash would need re-freezing all four hashes and keeping
`shallowforge/ir.py`'s `hamiltonian_hash` in sync (`tests/test_shallowforge_spec.py::test_hash_matches_chemcheck_for_same_operator`
asserts the two are the same function). `specs/tasks/02-chemcheck.md` line 5 still says "hashes
stable across machines".
