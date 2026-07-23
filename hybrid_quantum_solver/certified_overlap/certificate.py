"""Overlap certificate with provenance-chained gap inputs and invariant enforcement.

SPEC-21 (specs/SPEC_certified_overlap_bounds.md), invariants I1-I4. House style:
invariants are non-bypassable and failures are loud -- raise, never warn, never default.
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class GapCertificate:
    """A certified lower bound on E_1, consumed (never derived) by this module.

    The floor comes from the solver's existing gap machinery (temple_bounds /
    certified_gaps): oracle mode (caller-supplied exact E_1) or the premise-gated
    self mode eps_1 = theta_1 - sigma_1 (valid for M >= 6 -- see
    specs/SPEC_certified_gaps.md). This module only consumes the result; the
    derivation stays exclusively in the gap machinery (SPEC-21 section 7).

    ``e1_floor`` must be in the SAME energy frame as the Hamiltonian handed to
    ``certify_overlap`` (electronic vs total -- mixing frames is a silent-wrongness
    bug, which is why the frame is the caller's stated responsibility here).
    """

    e1_floor: float          # certified lower bound on E_1 (same frame as H)
    certificate_id: str      # provenance id, e.g. "oracle:fci" or "krylov_self_eps:M=8"
    source: str              # "oracle" | "krylov_self_eps"

    def __post_init__(self):
        if not self.certificate_id or not isinstance(self.certificate_id, str):
            raise ValueError(
                "GapCertificate requires a non-empty certificate_id (provenance is mandatory)."
            )
        if not np.isfinite(self.e1_floor):
            raise ValueError(
                f"GapCertificate.e1_floor must be finite, got {self.e1_floor}. "
                "A vacuous gap input must be rejected here, not propagated."
            )
        if self.source not in ("oracle", "krylov_self_eps"):
            raise ValueError(
                f"GapCertificate.source must be 'oracle' or 'krylov_self_eps', got {self.source!r}."
            )


@dataclass
class OverlapCertificate:
    """
    Certified overlap lower bound with provenance tracking.

    Invariants (non-bypassable, loud-failure):
      I1: No γ_min emitted without valid gap_certificate. Missing ⇒ raise.
      I2: If r ≥ δ, emit VACUOUS (γ_min=0), never a fabricated positive.
      I3: Every certificate embeds gap_certificate_id; serialization without it fails.
      I4: γ_min must never exceed 1 − ε_machine; if numerics give > 1, raise.

    ``e0_lower_temple`` is the Temple lower bound on E_0 computed from the SAME
    (u, gap-input) pair -- overlap floor and energy floor share provenance and are
    reported together (SPEC-21 section 2). None only when the Temple premise
    λ_u < e1_floor fails.

    ``cluster_size`` d: γ_min bounds |⟨u|ψ₀⟩| for a simple ground state (d=1, the
    SPEC-21 default) or ‖P_S u‖ for the lowest-d eigenspace (d>1, block Davis–Kahan,
    SPEC-21b). Invariant Ib: d must be an integer ≥ 1.
    """

    gamma_min: float
    lambda_u: float
    residual_norm: float
    gap_certificate_id: str
    conditional: bool = True  # always True: overlap is always conditional on gap
    bqp_threshold_note: Optional[str] = None
    vacuous: bool = False
    vacuous_reason: Optional[str] = None
    e0_lower_temple: Optional[float] = None
    cluster_size: int = 1

    def __post_init__(self):
        """Enforce invariants I1-I4 and Ib."""
        # I1: gap_certificate_id must be present (non-empty string)
        if not self.gap_certificate_id or not isinstance(self.gap_certificate_id, str):
            raise ValueError(
                "I1 violation: gap_certificate_id must be a non-empty string. "
                "No overlap certificate without a valid gap certificate."
            )

        # Ib: cluster_size must be an integer >= 1 (bool is not an acceptable int here)
        if isinstance(self.cluster_size, bool) or not isinstance(self.cluster_size, int):
            raise ValueError(
                f"Ib violation: cluster_size must be an int, got {type(self.cluster_size).__name__}."
            )
        if self.cluster_size < 1:
            raise ValueError(
                f"Ib violation: cluster_size must be >= 1, got {self.cluster_size}. "
                "A ground eigenspace of dimension < 1 is meaningless."
            )

        # I2: vacuous bounds must have gamma_min = 0 and a stated reason
        if self.vacuous:
            if self.gamma_min != 0:
                raise ValueError(
                    f"I2 violation: vacuous bound must have gamma_min = 0, got {self.gamma_min}. "
                    f"Reason: {self.vacuous_reason}"
                )
            if not self.vacuous_reason:
                raise ValueError(
                    "I2 violation: a vacuous certificate must state its reason."
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
        if self.gamma_min < 0:
            raise ValueError(
                f"gamma_min = {self.gamma_min} is negative; a non-vacuous certificate "
                "must carry a bound in [0, 1)."
            )

    def to_dict(self) -> dict:
        """Serialize. I3: refuses to serialize without the gap provenance chain."""
        if not self.gap_certificate_id or not isinstance(self.gap_certificate_id, str):
            raise ValueError(
                "I3 violation: cannot serialize an OverlapCertificate without its "
                "gap_certificate_id provenance chain."
            )
        return {
            "gamma_min": self.gamma_min,
            "lambda_u": self.lambda_u,
            "residual_norm": self.residual_norm,
            "gap_certificate_id": self.gap_certificate_id,
            "conditional": self.conditional,
            "bqp_threshold_note": self.bqp_threshold_note,
            "vacuous": self.vacuous,
            "vacuous_reason": self.vacuous_reason,
            "e0_lower_temple": self.e0_lower_temple,
            "cluster_size": self.cluster_size,
        }

    def __repr__(self) -> str:
        target = "|⟨u|ψ₀⟩|" if self.cluster_size == 1 else f"‖P_S u‖ (d={self.cluster_size})"
        if self.vacuous:
            return (
                f"OverlapCertificate(vacuous, target={target}, "
                f"reason={self.vacuous_reason}, gap_cert={self.gap_certificate_id})"
            )
        return (
            f"OverlapCertificate(γ_min={self.gamma_min:.6f} ≤ {target}, "
            f"λ_u={self.lambda_u:.6f}, r={self.residual_norm:.6e}, "
            f"E0≥{self.e0_lower_temple}, gap_cert={self.gap_certificate_id})"
        )
