# SPEC: ShallowForge — shallow-evolution circuit compiler (M1: metric harness)

**Slug:** `shallowforge` · **Depends on:** `hybrid_quantum_solver.trotter_krylov.build_trotter_
step`, qiskit transpile/Operator, `jsonschema`. Full PRD `specs/full/spec-shallow-evolution-
compiler.md`, tasks `specs/tasks/03-shallowforge.md`.

## Goal

Cut two-qubit gates per evolution step for chemistry Hamiltonians by 5–10× at fixed downstream
energy error. **This spec covers M1 only** — the measurement infrastructure that makes any later
claim trustworthy: a term-stream IR + content hash, a schema-valid provenance manifest, a
correctness verifier (compiled step vs exact `exp(-iHt)`), and the frozen first-order-Trotter
baselines. The R1–R5 technique rungs and the full ODMD CX@ε binary-search harness are M2+.

## The metric (defined once)

**CX@ε:** two-qubit gates per step such that the final ODMD energy error stays ≤ ε = 1.6 mHa. A
gate count is never reported without its ε (ADR-0007). M1 freezes the *baseline* CX/step; the
end-to-end ODMD ε loop that turns a raw count into a certified CX@ε is M2.

## Interface

```python
from shallowforge import (
    TermStreamIR, hamiltonian_hash,                      # IR + content hash (task 1)
    TransformEntry, build_manifest, validate_manifest,   # provenance manifest (task 2)
    compile_ir_to_circuit, count_cx, step_fidelity, baseline_cx_per_step,  # verifier (task 3-4)
)
```

## Acceptance gates (`tests/test_shallowforge_spec.py`)

1. **IR hash + round-trip (task 1).** `hamiltonian_hash` is invariant under reordering of the
   input terms; `TermStreamIR` round-trips through `to_dict`/`from_dict` and through
   `SparsePauliOp`; the hash matches `chemcheck.canonical_hamiltonian_sha256` for the same
   operator (cross-package agreement).
2. **Manifest emitter (task 2).** `build_manifest` output validates against
   `compiler-manifest.schema.json`; a lossless `TransformEntry` with non-zero ε is rejected;
   `totals.cx_at_epsilon_claim.epsilon_mha == 1.6`.
3. **Correctness verifier (task 3).** For T0–T2, a first-order step at the reference `dt` has
   `step_fidelity` well below 1, and increasing Suzuki order/reps drives fidelity ≥ 1 − 1e-8 —
   proving the verifier actually resolves Trotter error (no exact-evolution-in-disguise).
4. **Baseline frozen (task 4 headline).** `baseline_cx_per_step` reproduces the frozen
   `shallowforge/baselines.json` CX counts for all four golden systems, including the N₂
   CAS(6,6) figure quoted next to the ~6,500 wall.

## Out of scope / honest caveats

- **No compression yet.** M1 measures; it does not reduce gates. The 5–10× target is M2–M4.
- The full **CX@ε** harness (compiled circuit → full ODMD → binary-search step count to ε ≤ 1.6
  mHa) is deferred to M2; M1 freezes raw first-order baselines and the correctness gate only.
- `step_fidelity` is statevector-only (T0–T2); T3 (N₂, 12 qubits) is verified by CX count, not
  by full-unitary fidelity.
