"""Temple lower energy bound sharing provenance with the overlap certificate.

The inequality is the same one the repo's existing certified-bracket path applies to
Krylov Ritz states (temple_bounds.krylov_bracket, specs/SPEC_temple_bracket.md); here it
is exposed as a pure function of the (λ_u, r, β) triple so the overlap certificate can
report the energy floor derived from the SAME gap input (SPEC-21 section 2). The
derivation of β itself stays exclusively in the gap machinery.
"""


def temple_lower_bound(
    lambda_u: float,
    residual_norm: float,
    e1_floor: float,
) -> float:
    """
    Temple (1928) lower bound on E₀ from a trial state's mean and residual.

        E₀ ≥ λ_u − r² / (β − λ_u)    for any certified β ≤ E₁ with λ_u < β

    using that the trial-state variance ⟨u|(H − λ_u)²|u⟩ equals r² exactly.

    Args:
        lambda_u: λ_u = ⟨u|H|u⟩ (Rayleigh quotient of the normalized trial state)
        residual_norm: r = ‖(H − λ_u)u‖
        e1_floor: β = certified lower bound on E₁ (same energy frame as λ_u)

    Returns:
        Certified lower bound on the ground-state energy E₀.

    Raises:
        ValueError: if the Temple premise λ_u < β fails. Loud, never a default:
            with λ_u ≥ β the inequality is not valid and no number may be emitted.
    """
    if not lambda_u < e1_floor:
        raise ValueError(
            f"Temple premise violated: lambda_u = {lambda_u} must lie strictly below "
            f"the certified E1 floor beta = {e1_floor}. No bound is emitted."
        )
    return lambda_u - residual_norm**2 / (e1_floor - lambda_u)
