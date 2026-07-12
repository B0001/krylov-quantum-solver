"""ShallowForge — shallow-evolution circuit compiler for chemistry Hamiltonians (M1).

M1 is the measurement infrastructure: a canonical term-stream IR with a content hash, a
provenance manifest emitter (the ChemCheck compilation-audit artifact), and a verification
harness that proves a compiled Trotter step reproduces exact ``exp(-iHt)`` before any gate-count
claim is made. The metric is **CX@ε** — two-qubit gates per step at a fixed downstream energy
error — and a gate count is never reported without its ε (ADR-0007).

See specs/SPEC_shallowforge.md, tasks specs/tasks/03-shallowforge.md, full PRD
specs/full/spec-shallow-evolution-compiler.md. The R1–R5 technique rungs are M2+.
"""

from .ir import TermStreamIR, hamiltonian_hash
from .manifest import (
    MANIFEST_VERSION,
    TransformEntry,
    build_manifest,
    validate_manifest,
)
from .verify import (
    baseline_cx_per_step,
    compile_ir_to_circuit,
    count_cx,
    step_fidelity,
)

__all__ = [
    "TermStreamIR",
    "hamiltonian_hash",
    "MANIFEST_VERSION",
    "TransformEntry",
    "build_manifest",
    "validate_manifest",
    "compile_ir_to_circuit",
    "count_cx",
    "step_fidelity",
    "baseline_cx_per_step",
]
