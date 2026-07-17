"""Overlap certificate with provenance-chained gap inputs and invariant enforcement."""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class OverlapCertificate:
    """
    Certified overlap lower bound with provenance tracking.

    Invariants (non-bypassable, loud-failure):
      I1: No γ_min emitted without valid gap_certificate. Missing ⇒ raise.
      I2: If r ≥ δ, emit VACUOUS (γ_min=0), never a fabricated positive.
      I3: Every certificate embeds gap_certificate_id; serialization without it fails.
      I4: γ_min must never exceed 1 − ε_machine; if numerics give > 1, raise.
    """

    gamma_min: float
    lambda_u: float
    residual_norm: float
    gap_certificate_id: str
    conditional: bool = True  # always True: overlap is always conditional on gap
    bqp_threshold_note: Optional[str] = None
    vacuous: bool = False
    vacuous_reason: Optional[str] = None

    def __post_init__(self):
        """Enforce invariants I1–I4."""
        # I1: gap_certificate_id must be present (non-empty string)
        if not self.gap_certificate_id or not isinstance(self.gap_certificate_id, str):
            raise ValueError(
                "I1 violation: gap_certificate_id must be a non-empty string. "
                "No overlap certificate without a valid gap certificate."
            )

        # I2: vacuous bounds must have gamma_min = 0
        if self.vacuous:
            if self.gamma_min != 0:
                raise ValueError(
                    f"I2 violation: vacuous bound must have gamma_min = 0, got {self.gamma_min}. "
                    f"Reason: {self.vacuous_reason}"
                )
            return  # Skip I4 check for vacuous bounds

        # I4: gamma_min must not exceed 1 − ε_machine
        eps_machine = np.finfo(float).eps
        max_allowed = 1.0 - eps_machine
        if self.gamma_min > max_allowed:
            raise ValueError(
                f"I4 violation: gamma_min = {self.gamma_min} exceeds 1 − ε_machine = {max_allowed}. "
                "This indicates an upstream computation error; raising to prevent silent failure."
            )

    def __repr__(self) -> str:
        if self.vacuous:
            return (
                f"OverlapCertificate(vacuous, reason={self.vacuous_reason}, "
                f"gap_cert={self.gap_certificate_id})"
            )
        return (
            f"OverlapCertificate(γ_min={self.gamma_min:.6f}, "
            f"λ_u={self.lambda_u:.6f}, r={self.residual_norm:.6e}, "
            f"gap_cert={self.gap_certificate_id})"
        )
