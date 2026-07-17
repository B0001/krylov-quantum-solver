"""Davis–Kahan sin-θ bounds for guiding-state overlap certification."""


def gamma_min(residual_norm: float, spectral_gap: float) -> float:
    """
    Compute certified overlap lower bound via Davis–Kahan sin-θ theorem.

    Given residual r = ‖(H − λ_u)u‖ and certified spectral gap δ from the
    solver's existing gap machinery, compute:

        sin θ(u, ψ₀) ≤ r / δ
        ⇒ γ_min = sqrt(1 − r²/δ²) ≤ |⟨u|ψ₀⟩|

    Args:
        residual_norm: r = ‖(H − λ_u)u‖
        spectral_gap: δ = lower bound on dist(λ_u, spec(H) \\ {E₀})

    Returns:
        γ_min: lower bound on |⟨u|ψ₀⟩|

    Raises:
        NotImplementedError: This is a stub pending review of the full
            Davis–Kahan derivation and integration with the existing gap machinery.
    """
    raise NotImplementedError("gamma_min: Davis–Kahan bound implementation pending")


def vacuous_bound() -> tuple:
    """
    Return an explicit VACUOUS result when r ≥ δ.

    Returns:
        (gamma_min=0, reason="residual_norm >= spectral_gap")

    Raises:
        NotImplementedError: This is a stub pending full implementation.
    """
    raise NotImplementedError("vacuous_bound: implementation pending")
