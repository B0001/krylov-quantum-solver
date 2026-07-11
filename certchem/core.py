"""CertChem core — the certified_energy entry point (CertChem-M1 tasks 3-4, 7).

No new physics: this is plumbing that routes a molecule through the repo's validated
primitives (PySCF Hamiltonian → real-time quantum Krylov → Temple/Weinstein certified
bracket, ``temple_bounds.krylov_bracket``) and assembles the frozen contract types.

The variational floor is a hard chokepoint (architecture ADR-0001): the ONLY path from an
energy estimate to a ``CertifiedResult`` passes through :func:`floor_guard`, which raises
:class:`FloorViolationError` — it never returns a number below the certified lower bound.
"""

from __future__ import annotations

import math
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from .contract import (
    Bracket,
    Certificate,
    CertifiedResult,
    ConvergenceError,
    FloorViolationError,
    Mode,
)
from .limits import check_caps

#: Default Krylov dimension. >=6 is the gated boundary for the self-mode Temple premise
#: (see temple_bounds.py module docstring); 12 gives finite brackets on the golden suite.
DEFAULT_KRYLOV_DIM = 12

#: Variance below this counts as fully converged in the certificate.
_CONVERGED_VARIANCE = 1e-6

#: Energies are quantized to this many decimal places (Ha) before being placed in a result.
#: BLAS/LAPACK reduction order is not bit-reproducible across runs (~1e-16 jitter); quantizing
#: at 1e-12 Ha — six orders below chemical accuracy (1.6e-3 Ha) and below the solver tolerance —
#: makes the serialized result byte-identical, which is what the content-hash cache keys on
#: (ADR-0008). Honest limit: a value landing within ~1e-16 of a 1e-12 grid line could still
#: straddle it; this stabilizes caching, it is not a claim of bit-level BLAS reproducibility.
_ENERGY_QUANTUM_DIGITS = 12
_ENERGY_GRID = 10**_ENERGY_QUANTUM_DIGITS


def _q_nearest(x: float) -> float:
    """Quantize a point estimate onto the determinism grid (round to nearest)."""
    return round(x, _ENERGY_QUANTUM_DIGITS) if math.isfinite(x) else x


def _q_lower(x: float) -> float:
    """Quantize a certified LOWER bound: round toward -inf so it stays a valid lower bound."""
    return math.floor(x * _ENERGY_GRID) / _ENERGY_GRID if math.isfinite(x) else x


def _q_upper(x: float) -> float:
    """Quantize a certified UPPER bound: round toward +inf so it stays a valid upper bound."""
    return math.ceil(x * _ENERGY_GRID) / _ENERGY_GRID if math.isfinite(x) else x


def solver_version() -> str:
    """Version string pinned into every certificate (determinism / cache key, ADR-0008)."""
    for pkg in ("certchem", "hybrid_quantum_solver"):
        try:
            return version(pkg)
        except PackageNotFoundError:
            continue
    return "0+local"


def floor_guard(best_estimate_hartree: float, floor_hartree: float, *, tol: float = 1e-9) -> None:
    """The variational-floor chokepoint. Raises if the estimate sank below the certified floor.

    The Ritz estimate is a variational upper bound on E0, and E0 is itself >= the Temple lower
    bound, so a physically valid estimate is always >= ``floor_hartree``. Anything below it is a
    broken estimator (the "-hundreds of Hartree" failure this repo was rebuilt to catch) and
    must never be dressed up as a result.
    """
    if best_estimate_hartree < floor_hartree - tol:
        raise FloorViolationError(
            f"estimate {best_estimate_hartree} Ha is below the certified floor "
            f"{floor_hartree} Ha",
            diagnostics={
                "best_estimate_hartree": best_estimate_hartree,
                "floor_hartree": floor_hartree,
            },
        )


def _build_mh(molecule: Any, basis: str, cas: tuple[int, int]):
    """molecule (PySCF atom string) + basis + CAS → validated MolecularHamiltonian."""
    from hybrid_quantum_solver import build_molecular_hamiltonian

    n_electrons, n_orbitals = cas
    return build_molecular_hamiltonian(
        atom=molecule,
        basis=basis,
        active_electrons=n_electrons,
        active_orbitals=n_orbitals,
    )


def _estimate(mh, krylov_dim: int, eps: float | None) -> tuple[float, float, float, float]:
    """Run the certified estimator. Returns (best_estimate, lower, upper, variance), all Ha.

    Isolated so the floor-guard chokepoint has a single seam to test against.
    """
    from temple_bounds import krylov_bracket

    br = krylov_bracket(mh, krylov_dim, eps=eps)
    # Best estimate is the Ritz value == the variational upper bound.
    return br.upper, br.lower, br.upper, br.variance


def _estimate_gap(mh, krylov_dim: int, e1: float | None) -> tuple[float, float, float]:
    """Certified reachable-sector gap. Returns (best_gap, gap_lower, gap_upper), all Ha."""
    from certified_gaps import gap_bracket

    gb = gap_bracket(mh, krylov_dim, e1=e1)
    best = gb.theta1 - gb.theta0  # Ritz gap, provably inside [gap_lower, gap_upper]
    return best, gb.gap_lower, gb.gap_upper


def certified_energy(
    molecule: Any,
    basis: str,
    cas: tuple[int, int],
    mode: Mode = Mode.CERTIFIED,
    *,
    krylov_dim: int = DEFAULT_KRYLOV_DIM,
    eps: float | None = None,
) -> CertifiedResult | float:
    """Certified ground-state energy for ``molecule`` in ``basis`` over active space ``cas``.

    ``Mode.CERTIFIED`` returns a floor-checked :class:`CertifiedResult`; ``Mode.FAST`` returns a
    bare ``float`` point estimate with no guarantee (the return type makes the two impossible to
    confuse, ADR-0004). Raises :class:`CapExceededError` outside the validated envelope,
    :class:`FloorViolationError` on a sub-floor estimate, :class:`ConvergenceError` if no finite
    lower bound could be certified.
    """
    check_caps(molecule, basis, cas)
    mh = _build_mh(molecule, basis, cas)

    if mode is Mode.FAST:
        from hybrid_quantum_solver import QuantumKrylovSolver

        # eigenstates() already returns TOTAL energies (offset included).
        energies, _ = QuantumKrylovSolver(mh).eigenstates(krylov_dim, n_states=1)
        return float(energies[0])

    best, lower, upper, variance = _estimate(mh, krylov_dim, eps)

    if lower == float("-inf") or upper == float("inf"):
        raise ConvergenceError(
            "Temple premise gave a vacuous (infinite) bound; raise krylov_dim or supply eps",
            partial={"krylov_dim": krylov_dim, "upper": upper, "lower": lower},
        )

    # Chokepoint: nothing reaches Bracket construction without clearing the floor.
    floor_guard(best, lower)

    convergence = "converged" if variance < _CONVERGED_VARIANCE else "converged_marginal"
    return CertifiedResult(
        bracket=Bracket(
            lower_hartree=_q_lower(lower),
            upper_hartree=_q_upper(upper),
            best_estimate_hartree=_q_nearest(best),
        ),
        certificate=Certificate(
            method="temple_bound + variational_floor",
            floor_check="pass",
            krylov_dim=krylov_dim,
            convergence=convergence,
            solver_version=solver_version(),
            manifest=None,
        ),
    )


#: Gap types this module can certify. Only the fundamental reachable-sector gap in M1.
_SUPPORTED_GAP_TYPES = frozenset({"fundamental"})


def certified_gap(
    molecule: Any,
    basis: str,
    cas: tuple[int, int],
    gap_type: str = "fundamental",
    mode: Mode = Mode.CERTIFIED,
    *,
    krylov_dim: int = DEFAULT_KRYLOV_DIM,
    e1: float | None = None,
) -> CertifiedResult | float:
    """Certified bracket on a spectral gap (excitation energy) of ``molecule``.

    ``gap_type="fundamental"`` is Δ = E₁ − E₀ of the HF-reachable sector (the only gap M1
    certifies). ``Mode.CERTIFIED`` returns a floor-checked bracket built from interlacing +
    Temple + Weinstein (``certified_gaps.gap_bracket``); ``Mode.FAST`` returns the bare Ritz gap.
    The lower certificate's premise is only valid for ``krylov_dim >= 6`` (see certified_gaps.py).
    """
    if gap_type not in _SUPPORTED_GAP_TYPES:
        raise ValueError(
            f"gap_type {gap_type!r} not supported; M1 certifies {sorted(_SUPPORTED_GAP_TYPES)}"
        )
    check_caps(molecule, basis, cas)
    mh = _build_mh(molecule, basis, cas)

    best, lower, upper = _estimate_gap(mh, krylov_dim, e1)

    if mode is Mode.FAST:
        return float(best)

    if lower == float("-inf") or upper == float("inf"):
        raise ConvergenceError(
            "gap certificate is vacuous (infinite bound); raise krylov_dim (>=6) or supply e1",
            partial={"krylov_dim": krylov_dim, "gap_lower": lower, "gap_upper": upper},
        )

    # Chokepoint: a certified gap can never sink below its own certified lower bound.
    floor_guard(best, lower)

    return CertifiedResult(
        bracket=Bracket(
            lower_hartree=_q_lower(lower),
            upper_hartree=_q_upper(upper),
            best_estimate_hartree=_q_nearest(best),
        ),
        certificate=Certificate(
            method="certified_gap: interlacing + temple + weinstein",
            floor_check="pass",
            krylov_dim=krylov_dim,
            convergence="converged",
            solver_version=solver_version(),
            manifest=None,
        ),
    )
