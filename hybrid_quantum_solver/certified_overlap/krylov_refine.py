"""Lanczos-chained refinement of the SPEC-21 overlap bound.

The direct Davis-Kahan certificate bounds |<u|psi_0>| using u's OWN residual, so it goes VACUOUS
whenever r_u >= delta_u -- which for a Hartree-Fock guiding state happens exactly on the
multireference systems the certificate is wanted for.

Chaining fixes that. With v the Krylov ground Ritz vector, the angle triangle inequality on the
unit sphere gives

    theta(u, psi_0)  <=  theta(u, v) + theta(v, psi_0)

and Davis-Kahan applied to v (not u) bounds the second term by arcsin(r_v / delta_v). Hence

    |<u|psi_0>|  =  cos theta(u, psi_0)  >=  cos( theta(u,v) + arcsin(r_v/delta_v) )

valid whenever the angle sum stays below pi/2. Because v is a Ritz vector of a Krylov space that
already contains u, its residual is orders of magnitude smaller than u's, so the chained bound is
dramatically tighter -- and it costs NO extra measurements: <u|v> is a linear combination of the
Krylov overlap-matrix entries <u|e^{-ik dt H}|u> the solver has already formed.

Measured on a symmetry-clean set (specs/SPEC_symmetry_reachability.md defines what that means):
valid everywhere, strictly tighter than the direct bound wherever the direct bound is non-vacuous,
and it rescues every vacuous case -- linear H6 R=1.0 at M=12 goes direct VACUOUS -> chained 0.945
against an exact 0.950.

TWO PROPERTIES CALLERS MUST NOT ASSUME:
  * It SATURATES. When v is converged to machine precision (r_v/delta_v ~ 1e-16) the bound equals
    the exact overlap, and float rounding puts it 1-3 ulp either side. Compare with slack, never
    with a bare ``<=``.
  * It is NOT monotone in the Krylov dimension M. More depth does not always mean a better bound.
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
from scipy.sparse import spmatrix

# Saturation slack: at a machine-converged Ritz vector the bound IS the exact overlap, so rounding
# can put it a few ulp above. Anything beyond this is a real violation, not arithmetic.
SATURATION_SLACK = 1e-14


def refine_via_lanczos(
    H: Union[np.ndarray, spmatrix],
    u: np.ndarray,
    v: np.ndarray,
    e1_floor: float,
) -> Optional[float]:
    """Chained lower bound on |<u|psi_0>| via the Krylov ground Ritz vector ``v``.

    Args:
        H: Hermitian matrix, same energy frame as ``e1_floor``.
        u: normalized guiding state (typically Hartree-Fock).
        v: normalized Krylov ground Ritz vector -- the refinement anchor.
        e1_floor: certified lower bound on E_1, the first level above the ground state. Oracle or
            the repo's premise-gated self-mode floor; the caller owns that provenance.

    Returns:
        gamma_chain in [0, 1], or ``None`` when the bound is VACUOUS -- either the Temple premise
        fails (lambda_v >= e1_floor), or the Ritz residual does not resolve the gap
        (r_v >= delta_v), or the chained angle sum reaches pi/2. ``None`` means "no statement",
        never a fabricated zero-information positive (invariant I2 in spirit).

    Raises:
        ValueError: if ``u`` or ``v`` is not normalized -- a silently unnormalized input would
            corrupt both the Rayleigh quotient and the overlap, so it fails loudly.
    """
    u = np.asarray(u, dtype=complex).ravel()
    v = np.asarray(v, dtype=complex).ravel()
    for name, vec in (("u", u), ("v", v)):
        norm = float(np.linalg.norm(vec))
        if not np.isclose(norm, 1.0, atol=1e-8):
            raise ValueError(f"{name} must be normalized, got ||{name}|| = {norm}.")

    Hv = H @ v
    lambda_v = float(np.real(np.vdot(v, Hv)))
    delta_v = float(e1_floor) - lambda_v
    if delta_v <= 0:
        return None                                  # Temple premise fails: no separation to use
    r_v = float(np.linalg.norm(Hv - lambda_v * v))
    if r_v >= delta_v:
        return None                                  # residual does not resolve the gap

    theta_uv = float(np.arccos(min(1.0, abs(complex(np.vdot(u, v))))))
    total = theta_uv + float(np.arcsin(r_v / delta_v))
    if total >= 0.5 * np.pi:
        return None                                  # chained angle exceeds a right angle
    return float(np.cos(total))
