"""Davis–Kahan sin-θ bounds for guiding-state overlap certification (SPEC-21).

Mathematical core. For Hermitian H with simple ground state ψ₀, normalized trial u,
λ_u = ⟨u|H|u⟩, r = ‖(H − λ_u)u‖, and any certified δ with
0 < δ ≤ dist(λ_u, spec(H) \\ {E₀}):

    expanding u in the eigenbasis, r² = Σᵢ |cᵢ|²(λᵢ − λ_u)² ≥ δ² Σ_{i≥1} |cᵢ|²
    = δ² sin²θ(u, ψ₀), hence sin θ ≤ r/δ and

        γ_min = sqrt(1 − r²/δ²) ≤ cos θ = |⟨u|ψ₀⟩|.

When λ_u < E₁, dist(λ_u, spec(H) \\ {E₀}) = E₁ − λ_u, so any certified E₁ floor β
gives the valid (conservative) separation δ = β − λ_u.
"""

from typing import Optional, Union

import numpy as np
from scipy.sparse import spmatrix

from .certificate import GapCertificate, OverlapCertificate
from .rayleigh import rayleigh_quotient
from .residual import residual_norm
from .temple import temple_lower_bound

_NORMALIZATION_TOL = 1e-10


def gamma_min(residual: float, delta: float) -> float:
    """
    Certified overlap lower bound γ_min = sqrt(1 − r²/δ²) from Davis–Kahan sin-θ.

    Pure math on the (r, δ) pair; premise checks are loud. Callers wanting the
    I1/I2-enforced pipeline use ``certify_overlap`` instead.

    Args:
        residual: r = ‖(H − λ_u)u‖, must be ≥ 0
        delta: certified separation δ ≤ dist(λ_u, spec(H) \\ {E₀}), must be > 0

    Returns:
        γ_min ∈ [0, 1]: certified lower bound on |⟨u|ψ₀⟩|

    Raises:
        ValueError: if δ ≤ 0, r < 0, or r ≥ δ (the vacuous regime — callers must
            emit an explicit VACUOUS certificate instead, invariant I2).
    """
    if delta <= 0:
        raise ValueError(f"delta must be positive, got {delta}.")
    if residual < 0:
        raise ValueError(f"residual norm cannot be negative, got {residual}.")
    if residual >= delta:
        raise ValueError(
            f"vacuous regime: r = {residual} >= delta = {delta}. "
            "Emit a VACUOUS certificate (I2), never call gamma_min here."
        )
    ratio = residual / delta
    return float(np.sqrt(1.0 - ratio * ratio))


def certify_overlap(
    H: Union[np.ndarray, spmatrix],
    u: np.ndarray,
    gap_certificate: Optional[GapCertificate],
    n_qubits: Optional[int] = None,
) -> OverlapCertificate:
    """
    Certify a lower bound on |⟨u|ψ₀⟩| for a guiding state u, conditional on a gap input.

    The single entry point of SPEC-21. Enforces the invariants:
      I1 — no gap certificate ⇒ raise (never warn, never default a gap);
      I2 — r ≥ δ (or λ_u ≥ β) ⇒ explicit VACUOUS certificate with γ_min = 0;
      I4 — numerically impossible γ_min ⇒ raise (checked in OverlapCertificate).

    Also reports the Temple lower bound on E₀ from the same (u, gap-input) pair —
    shared provenance, reported together (SPEC-21 section 2).

    Args:
        H: Hermitian matrix (dense or sparse), same energy frame as the gap input
        u: normalized trial state (raises if ‖u‖ deviates from 1)
        gap_certificate: certified E₁ floor from the existing gap machinery
        n_qubits: optional system size for the guided-LH 1/poly framing note

    Returns:
        OverlapCertificate (possibly VACUOUS — check ``.vacuous`` before quoting γ_min).
    """
    # I1 — non-bypassable, loud
    if gap_certificate is None:
        raise ValueError(
            "I1 violation: certify_overlap called without a gap certificate. "
            "A certified overlap bound is conditional on a certified gap input; "
            "there is no default and no warning path."
        )
    if not isinstance(gap_certificate, GapCertificate):
        raise TypeError(
            f"gap_certificate must be a GapCertificate, got {type(gap_certificate).__name__}. "
            "Derive it from the existing gap machinery (temple_bounds / certified_gaps)."
        )

    norm = float(np.linalg.norm(u))
    if abs(norm - 1.0) > _NORMALIZATION_TOL:
        raise ValueError(
            f"trial state must be normalized: got ||u|| = {norm}. "
            "Normalize explicitly at the call site; silent renormalization hides bugs."
        )

    lam = rayleigh_quotient(H, u)
    r = residual_norm(H, u, lam)
    beta = gap_certificate.e1_floor
    delta = beta - lam

    if delta <= 0:
        # Temple premise (lambda_u < beta) also fails here: no energy floor either.
        return OverlapCertificate(
            gamma_min=0.0,
            lambda_u=lam,
            residual_norm=r,
            gap_certificate_id=gap_certificate.certificate_id,
            vacuous=True,
            vacuous_reason=(
                f"no positive separation: lambda_u = {lam} >= certified E1 floor {beta}"
            ),
            e0_lower_temple=None,
        )

    # Temple holds whenever lambda_u < beta, independent of the r-vs-delta comparison.
    e0_temple = temple_lower_bound(lam, r, beta)

    if r >= delta:
        return OverlapCertificate(
            gamma_min=0.0,
            lambda_u=lam,
            residual_norm=r,
            gap_certificate_id=gap_certificate.certificate_id,
            vacuous=True,
            vacuous_reason=f"residual r = {r} >= separation delta = {delta}",
            e0_lower_temple=e0_temple,
        )

    gm = gamma_min(r, delta)

    if n_qubits is not None:
        threshold = 1.0 / n_qubits
        rel = ">=" if gm >= threshold else "<"
        note = (
            f"gamma_min {rel} 1/n for n = {n_qubits} qubits "
            "(guided-LH 1/poly framing, arXiv:2111.09079)"
        )
    else:
        note = "n unstated: no guided-LH 1/poly comparison made"

    return OverlapCertificate(
        gamma_min=gm,
        lambda_u=lam,
        residual_norm=r,
        gap_certificate_id=gap_certificate.certificate_id,
        bqp_threshold_note=note,
        e0_lower_temple=e0_temple,
    )
