"""The interval-dominance rule (task 1). Pure, zero solver deps.

Soundness proof (one paragraph): a candidate's certified bracket ``[lo, hi]`` encloses its true
property value ``v`` (CertChem containment, ADR-0001): ``lo <= v <= hi``. The target region is
``[t_lo, t_hi]``; a *hit* is a candidate with ``v in [t_lo, t_hi]``. We ELIMINATE iff the bracket
is strictly disjoint from the target (``hi < t_lo`` or ``lo > t_hi``). If a candidate is
eliminated, its whole bracket — and therefore ``v`` — lies strictly outside the target, so it is
not a hit. Contrapositive: every hit survives. Hence **zero false eliminations**. Boundary contact
(``hi == t_lo`` or ``lo == t_hi``) is treated as overlap (kept), because ``v`` could equal the
boundary — conservative, preserving soundness.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum


@dataclass(frozen=True)
class Interval:
    """A closed interval ``[lo, hi]`` with ``lo <= hi``."""

    lo: float
    hi: float

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError(f"interval lo ({self.lo}) > hi ({self.hi})")

    @property
    def width(self) -> float:
        return self.hi - self.lo

    def contains(self, x: float) -> bool:
        return self.lo <= x <= self.hi


class Verdict(Enum):
    """Outcome of classifying one candidate bracket against the target region."""

    ELIMINATED = "eliminated"  # bracket strictly excludes the target -> provably not a hit
    CONFIRMED = "confirmed"    # bracket fully inside the target -> provably a hit
    UNDECIDED = "undecided"    # bracket straddles a target boundary -> refine to decide


def overlaps(bracket: Interval, target: Interval) -> bool:
    """True if the closed intervals share any point (boundary contact counts as overlap)."""
    return bracket.lo <= target.hi and target.lo <= bracket.hi


def classify(bracket: Interval, target: Interval) -> Verdict:
    """Classify a candidate's certified ``bracket`` against the ``target`` region.

    ELIMINATED iff strictly disjoint; CONFIRMED iff the bracket is contained in the target;
    otherwise UNDECIDED. Never eliminates a bracket that overlaps the target — the soundness core.
    """
    if not overlaps(bracket, target):
        return Verdict.ELIMINATED
    if target.lo <= bracket.lo and bracket.hi <= target.hi:
        return Verdict.CONFIRMED
    return Verdict.UNDECIDED


def target_from_threshold(*, at_most: float | None = None, at_least: float | None = None) -> Interval:
    """Build a target region from a one-sided threshold (e.g. "property <= at_most")."""
    lo = -math.inf if at_least is None else at_least
    hi = math.inf if at_most is None else at_most
    return Interval(lo, hi)
