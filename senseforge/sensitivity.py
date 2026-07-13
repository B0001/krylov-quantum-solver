#!/usr/bin/env python3
"""
senseforge.sensitivity -- certified central-difference sensitivities + plateau detection.

S = d(gap)/d(axis) via central finite differences over a uniform grid, with the bracket bounds
propagated through the difference by interval arithmetic (worst-case corners of the two
brackets, not just the point-estimate slope): for a central difference (y_plus - y_minus)/(2h),
interval subtraction gives [lower_plus - upper_minus, upper_plus - lower_minus] / (2h). This is
gated on synthetic data (PRD sec 5 / task 5): a known quadratic gap curve's recovered slope must
have its analytic derivative inside the propagated interval.

Second derivative (central, 3-point) flags plateaus: PRD sec 5, "prefer plateaus of high |S| over
knife-edge points" -- a small |d2(gap)/dx2| near a high-|S| point is the stable operating point.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

from certchem.contract import CertifiedResult


@dataclass(frozen=True)
class Sensitivity:
    """Central-difference sensitivity at grid point ``x`` (an interior point, needs neighbours)."""

    x: float
    slope: float          # point-estimate d(gap)/dx from best_estimate values
    slope_lower: float    # propagated lower bound
    slope_upper: float    # propagated upper bound
    second_derivative: float  # 3-point central second difference around x itself

    @property
    def width(self) -> float:
        return self.slope_upper - self.slope_lower


def _central_slope_bracket(lo_minus: float, up_minus: float, lo_plus: float, up_plus: float,
                           h: float) -> tuple:
    """Interval-arithmetic central difference: [lo_plus - up_minus, up_plus - lo_minus] / (2h)."""
    lo = (lo_plus - up_minus) / (2.0 * h)
    hi = (up_plus - lo_minus) / (2.0 * h)
    return lo, hi


def certified_central_differences(xs: Sequence[float],
                                  results: Sequence[CertifiedResult]) -> List[Sensitivity]:
    """Sensitivities at every INTERIOR point of a uniform grid ``xs`` (length >= 3).

    ``results[i]`` is the certified gap bracket at ``xs[i]``. Requires a uniform step (as
    ``SweepConfig.grid()`` produces); raises ``ValueError`` otherwise.
    """
    n = len(xs)
    if n < 3:
        raise ValueError(f"need >= 3 grid points for a central difference, got {n}")
    if len(results) != n:
        raise ValueError(f"xs ({n}) and results ({len(results)}) length mismatch")
    steps = [xs[i + 1] - xs[i] for i in range(n - 1)]
    h = steps[0]
    if any(abs(s - h) > 1e-9 * max(abs(h), 1.0) for s in steps):
        raise ValueError("certified_central_differences requires a uniform grid step")

    out = []
    for i in range(1, n - 1):
        b_minus, b_plus = results[i - 1].bracket, results[i + 1].bracket
        slope = (b_plus.best_estimate_hartree - b_minus.best_estimate_hartree) / (2.0 * h)
        lo, hi = _central_slope_bracket(b_minus.lower_hartree, b_minus.upper_hartree,
                                        b_plus.lower_hartree, b_plus.upper_hartree, h)
        # 3-point central second difference around x_i itself (x_{i-1}, x_i, x_{i+1} are all
        # in bounds for every i in this loop's range, so this is never partial).
        y_minus = results[i - 1].bracket.best_estimate_hartree
        y0 = results[i].bracket.best_estimate_hartree
        y_plus = results[i + 1].bracket.best_estimate_hartree
        second = (y_plus - 2.0 * y0 + y_minus) / (h * h)
        out.append(Sensitivity(x=xs[i], slope=slope, slope_lower=lo, slope_upper=hi,
                               second_derivative=second))
    return out
