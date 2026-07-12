"""ChemCheck — hardware chemistry honesty benchmark (M1 + Mode A + floor detector).

A neutral, versioned scorecard: device spec sheet → PASS/FAIL + fidelity headroom per tier,
plus the anti-fraud floor detector. See specs/SPEC_chemcheck.md and the full PRD
specs/full/spec-hardware-honesty-benchmark.md.

The tier registry (:data:`TIERS`) imports with zero solver dependencies; anything that touches
a Hamiltonian (:func:`build_tier_hamiltonian`, :func:`canonical_hamiltonian_sha256`) imports
PySCF/qiskit lazily.
"""

from .budget import (
    ROUTING_OVERHEAD,
    expected_total_error,
    headroom_factor,
    required_two_qubit_error,
    routing_overhead,
)
from .scorecard import mode_b_energy_verdict, render_markdown, score_mode_a, scoring_code_hash
from .submission import validate_submission
from .tiers import (
    BENCHMARK_VERSION,
    TIERS,
    Tier,
    build_tier_hamiltonian,
    canonical_hamiltonian_sha256,
    recompute_tier_reference,
)

__all__ = [
    "BENCHMARK_VERSION",
    "TIERS",
    "Tier",
    "build_tier_hamiltonian",
    "canonical_hamiltonian_sha256",
    "recompute_tier_reference",
    "ROUTING_OVERHEAD",
    "routing_overhead",
    "expected_total_error",
    "required_two_qubit_error",
    "headroom_factor",
    "validate_submission",
    "score_mode_a",
    "render_markdown",
    "mode_b_energy_verdict",
    "scoring_code_hash",
]
