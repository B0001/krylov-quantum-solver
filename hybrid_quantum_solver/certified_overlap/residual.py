"""Residual norm computation for trial states."""

import numpy as np
from typing import Union
from scipy.sparse import spmatrix


def residual_norm(
    H: Union[np.ndarray, spmatrix], u: np.ndarray, lambda_u: float
) -> float:
    """
    Compute the residual norm r = ‖(H − λ_u)u‖.

    Args:
        H: Hermitian matrix (dense ndarray or sparse matrix)
        u: normalized trial state (complex array)
        lambda_u: Rayleigh quotient ⟨u|H|u⟩

    Returns:
        r: float, the Euclidean norm of the residual vector (H − λ_u)u
    """
    # (H - lambda_u * I) @ u
    residual = H @ u - lambda_u * u
    # ‖residual‖_2
    r = np.linalg.norm(residual)
    return r
