"""ScreenLoop — bracket-aware candidate screening with zero false eliminations (portfolio #19).

Eliminate candidates by certified interval dominance: a candidate whose certified bracket
strictly excludes the target region cannot be a hit (its true value is inside the bracket), so it
is safe to drop — and a genuine hit is provably never eliminated. A point-estimate screener has no
such guarantee. See specs/SPEC_screenloop.md and tasks specs/tasks/14-screening-loop.md.

The core (:func:`classify`, :class:`Interval`) is pure interval logic with zero solver deps; the
:class:`BoundedOracle` protocol adapts any interval source (CertChem certified brackets, a
conformal-ML model, a synthetic mock).
"""

from .oracle import BoundedOracle, CertchemEnergyOracle, SyntheticOracle
from .pruning import Interval, Verdict, classify, overlaps
from .screen import ScreenResult, point_estimate_screen, screen

__all__ = [
    "Interval",
    "Verdict",
    "classify",
    "overlaps",
    "BoundedOracle",
    "SyntheticOracle",
    "CertchemEnergyOracle",
    "screen",
    "point_estimate_screen",
    "ScreenResult",
]
