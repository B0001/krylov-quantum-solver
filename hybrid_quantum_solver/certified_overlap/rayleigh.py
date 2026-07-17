"""Rayleigh quotient computation for trial states."""

import numpy as np
from typing import Union
from scipy.sparse import spmatrix


def rayleigh_quotient(
    H: Union[np.ndarray, spmatrix], u: np.ndarray
) -> float:
    """
    Compute the Rayleigh quotient λ_u = ⟨u|H|u⟩ / ⟨u|u⟩.

    Assumes u is normalized: ⟨u|u⟩ = 1.

    Args:
        H: Hermitian matrix (dense ndarray or sparse matrix)
        u: normalized trial state (complex array)

    Returns:
        λ_u: float, the Rayleigh quotient
    """
    # H @ u
    Hu = H @ u
    # ⟨u|H|u⟩ = u† @ (H @ u)
    lambda_u = np.real(np.conj(u) @ Hu)
    return lambda_u
