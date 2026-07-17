"""Temple/Lehmann-type lower energy bound (shared provenance with overlap certificate)."""


def temple_lower_bound(
    rayleigh_quotient: float,
    residual_norm: float,
    upper_bound_e1: float,
) -> float:
    """
    Compute Temple/Lehmann lower bound on E₀.

        E₀ ≥ λ_u − r² / (β − λ_u)    for certified β ≤ E₁ with λ_u < β

    Args:
        rayleigh_quotient: λ_u = ⟨u|H|u⟩
        residual_norm: r = ‖(H − λ_u)u‖
        upper_bound_e1: β = certified upper bound on E₁ (first excited)

    Returns:
        E₀_lower: lower bound on ground-state energy

    Raises:
        NotImplementedError: This stub wraps the existing solver's certified
            lower-bound path; implementation details TBD pending integration.
    """
    raise NotImplementedError(
        "temple_lower_bound: wraps existing certified lower-bound path; implementation pending"
    )
