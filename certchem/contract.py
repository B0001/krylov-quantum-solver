"""Frozen public contract types for the CertChem solver core.

A verbatim freeze of ``architecture/interfaces/solver-library-contract.md``. This module
must import with **zero solver dependencies** (no numpy/pyscf/qiskit) so any caller can
build and inspect results, and so the type layer can be reasoned about in isolation.

Invariants the types enforce (architecture ADR-0001/0004):
  * A ``CertifiedResult`` only exists if the floor check passed — the physics layer must
    never construct one otherwise (enforced by the pipeline, not this module).
  * ``bracket.lower_hartree <= best_estimate_hartree <= upper_hartree`` — enforced here at
    construction time, so an out-of-order bracket is impossible to hold.
  * ``Mode.FAST`` never yields a ``Bracket`` — the return type union makes blurring a
    point estimate and a certified result a type error.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any


class Mode(Enum):
    """Estimation mode. ``FAST`` returns a bare float; ``CERTIFIED`` a floor-checked bracket."""

    FAST = "fast"
    CERTIFIED = "certified"


@dataclass(frozen=True)
class Bracket:
    """A two-sided energy bracket in Hartree, enclosing a best point estimate."""

    lower_hartree: float
    upper_hartree: float
    best_estimate_hartree: float

    def __post_init__(self) -> None:
        # Invariant 2: lower <= best <= upper, always. An out-of-order bracket is a bug
        # upstream; refuse to construct one rather than hand a caller a false guarantee.
        if not (self.lower_hartree <= self.best_estimate_hartree <= self.upper_hartree):
            raise ValueError(
                "Bracket violates lower <= best_estimate <= upper: "
                f"lower={self.lower_hartree}, best={self.best_estimate_hartree}, "
                f"upper={self.upper_hartree}"
            )

    @property
    def width(self) -> float:
        """Bracket width in Hartree (upper - lower); a non-negative uncertainty measure."""
        return self.upper_hartree - self.lower_hartree


@dataclass(frozen=True)
class Certificate:
    """Provenance for a certified result. Holding one implies ``floor_check == "pass"``."""

    method: str
    floor_check: str
    krylov_dim: int
    convergence: str
    solver_version: str
    manifest: dict[str, Any] | None = None


@dataclass(frozen=True)
class CertifiedResult:
    """A floor-checked bracket plus its certificate. Its existence is the guarantee."""

    bracket: Bracket
    certificate: Certificate


class FloorViolationError(Exception):
    """Estimate fell below the variational floor. Never returns a number.

    Carries diagnostics so callers can see how far below the floor the estimate landed.
    """

    def __init__(self, message: str, *, diagnostics: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.diagnostics: dict[str, Any] = diagnostics or {}


class CapExceededError(Exception):
    """System is outside the validated envelope. Names the violated cap."""

    def __init__(self, message: str, *, cap: str) -> None:
        super().__init__(message)
        self.cap = cap


class ConvergenceError(Exception):
    """ODMD signal insufficient at the configured Krylov dim. No estimate produced.

    Carries whatever partial data was gathered so the caller can retry with more resources.
    """

    def __init__(self, message: str, *, partial: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.partial: dict[str, Any] = partial or {}
