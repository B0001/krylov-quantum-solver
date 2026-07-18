"""Lanczos-chained refinement of overlap bounds (stub)."""


def refine_via_lanczos():
    """
    Refine overlap bound by running Lanczos from the guiding state u.

    Ritz vectors produced by Lanczos have monotonically non-increasing residuals.
    Each iterate yields a possibly-tighter γ_min for the Ritz vector, chained back
    to u via computable inner products. This makes the certificate improvable at
    the cost of matvecs, matching the solver's existing convergence loop.

    Raises:
        NotImplementedError: Stub pending Lanczos integration.
    """
    raise NotImplementedError("refine_via_lanczos: implementation pending")
