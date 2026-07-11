"""CertChem — certified chemistry solver-as-library (CertChem-M1).

The public boundary consumed by CertChem workers and SenseForge (architecture ADR-0005).
This package is the foundation everything downstream imports; see
``specs/tasks/01-certchem-m1.md`` and ``architecture/interfaces/solver-library-contract.md``.

Only the contract layer (types + limits) is frozen here so far — it imports with **zero
solver dependencies** by design, so lightweight callers can validate inputs and construct
results without pulling in PySCF/qiskit.
"""

from .contract import (
    Bracket,
    CapExceededError,
    Certificate,
    CertifiedResult,
    ConvergenceError,
    FloorViolationError,
    Mode,
)
from .core import certified_energy, certified_gap, floor_guard, solver_version
from .limits import ALLOWED_BASES, MAX_SPIN_ORBITALS, check_caps

__all__ = [
    "Mode",
    "certified_energy",
    "certified_gap",
    "floor_guard",
    "solver_version",
    "Bracket",
    "Certificate",
    "CertifiedResult",
    "FloorViolationError",
    "CapExceededError",
    "ConvergenceError",
    "MAX_SPIN_ORBITALS",
    "ALLOWED_BASES",
    "check_caps",
]
