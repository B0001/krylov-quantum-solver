"""
Certified guiding-state overlap bounds (SPEC-21).

Extends the solver's two-sided certification machinery to certify overlap
γ_min ≤ |⟨u|ψ₀⟩| given a guiding state u and a certified gap input, using
Davis–Kahan sin-θ theory.

Reference: arXiv:2111.09079, arXiv:2207.10097, arXiv:2207.10250 (guided-LH);
           O'Gorman–Irani–Whitfield–Fefferman arXiv:2103.08215 (ES hardness).
"""

from .rayleigh import rayleigh_quotient
from .residual import residual_norm
from .davis_kahan import gamma_min, vacuous_bound
from .certificate import OverlapCertificate

__all__ = [
    "rayleigh_quotient",
    "residual_norm",
    "gamma_min",
    "vacuous_bound",
    "OverlapCertificate",
]
