"""Interval-source oracles (task 6). A ``BoundedOracle`` is the only thing the loop depends on.

The soundness guarantee is *relative to* the oracle's brackets actually containing truth — that is
CertChem's certified contract for :class:`CertchemEnergyOracle`, and a construction invariant for
:class:`SyntheticOracle` (bracket half-width always exceeds the injected estimate bias). Any other
interval source — a conformal-ML model, a cheaper surrogate — plugs in by implementing the
protocol.
"""

from __future__ import annotations

import hashlib
from typing import Any, Protocol, runtime_checkable

from .pruning import Interval


@runtime_checkable
class BoundedOracle(Protocol):
    """Yields a certified bracket (and a point estimate) for a candidate at a precision level."""

    def bracket(self, candidate: Any, precision: int) -> Interval:
        """Certified interval enclosing the candidate's true property. Tightens with precision."""
        ...

    def point_estimate(self, candidate: Any, precision: int) -> float:
        """Best single-number guess — NOT guaranteed to be the true value (the baseline uses it)."""
        ...

    def cost(self, precision: int) -> float:
        """Relative evaluation cost at this precision (higher precision costs more)."""
        ...


def _sign(candidate: Any, seed: int) -> float:
    """Deterministic +/-1 bias direction for a candidate."""
    h = hashlib.sha256(f"{seed}:{candidate}".encode()).digest()[0]
    return 1.0 if h & 1 else -1.0


class SyntheticOracle:
    """Deterministic mock oracle over a known ground truth (for gates + the baseline contrast).

    At precision ``p`` the certified half-width is ``w0 * ratio**p`` and the point estimate is
    biased by ``b0 * ratio**p`` (sign per candidate). Since ``w0 > b0``, the bracket ALWAYS
    contains the truth — soundness holds — while the point estimate can sit outside the target at
    low precision, which is exactly what trips the point-estimate baseline into false eliminations.
    """

    def __init__(
        self,
        truths: dict[Any, float],
        *,
        w0: float = 1.0,
        b0: float = 0.6,
        ratio: float = 0.4,
        seed: int = 0,
    ) -> None:
        if not w0 > b0:
            raise ValueError("w0 must exceed b0 so brackets always contain truth")
        self.truths = dict(truths)
        self.w0, self.b0, self.ratio, self.seed = w0, b0, ratio, seed

    def _bias(self, candidate: Any, precision: int) -> float:
        return self.b0 * self.ratio**precision * _sign(candidate, self.seed)

    def point_estimate(self, candidate: Any, precision: int) -> float:
        return self.truths[candidate] + self._bias(candidate, precision)

    def bracket(self, candidate: Any, precision: int) -> Interval:
        pe = self.point_estimate(candidate, precision)
        w = self.w0 * self.ratio**precision
        return Interval(pe - w, pe + w)

    def cost(self, precision: int) -> float:
        return 2.0**precision  # tighter precision is exponentially more expensive


class CertchemEnergyOracle:
    """Real oracle: CertChem certified energy of ``(molecule, basis, cas)`` candidates.

    ``precision`` maps to Krylov dimension. Smoke-tested only — the headline gates use
    :class:`SyntheticOracle` so the zero-false-elimination invariant is checked over a large space
    deterministically and fast.
    """

    def __init__(self, *, base_krylov_dim: int = 8, step: int = 4) -> None:
        self.base_krylov_dim = base_krylov_dim
        self.step = step

    def _krylov_dim(self, precision: int) -> int:
        return self.base_krylov_dim + self.step * precision

    def bracket(self, candidate: Any, precision: int) -> Interval:
        from certchem import certified_energy

        molecule, basis, cas = candidate
        result = certified_energy(molecule, basis, cas, krylov_dim=self._krylov_dim(precision))
        return Interval(result.bracket.lower_hartree, result.bracket.upper_hartree)

    def point_estimate(self, candidate: Any, precision: int) -> float:
        from certchem import Mode, certified_energy

        molecule, basis, cas = candidate
        return float(
            certified_energy(molecule, basis, cas, Mode.FAST, krylov_dim=self._krylov_dim(precision))
        )

    def cost(self, precision: int) -> float:
        return float(self._krylov_dim(precision))
