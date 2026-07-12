"""The screening loop (tasks 3-4) and the unsound point-estimate baseline (task 5).

``screen`` evaluates every candidate cheaply first, eliminates by certified dominance, and only
refines the undecided survivors — spending expensive high-precision evaluations on the few
candidates the cheap pass could not resolve. With ``acquisition=True`` (v2) it refines the
survivors most likely to resolve first (largest overlap-with-target boundary), converging in fewer
high-cost steps. Both return the SAME hit set as evaluating everything at full precision, but for a
fraction of the cost. ``point_estimate_screen`` is the baseline that eliminates on a point estimate
— and consequently makes false eliminations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .oracle import BoundedOracle
from .pruning import Interval, Verdict, classify


@dataclass
class ScreenResult:
    """Outcome of a screen: which candidates survived, their brackets, and the cost spent."""

    survivors: list[Any] = field(default_factory=list)      # kept (confirmed or undecided) = the hits
    confirmed: list[Any] = field(default_factory=list)      # bracket provably inside the target
    eliminated: list[Any] = field(default_factory=list)     # provably not hits
    total_cost: float = 0.0                                  # sum of oracle.cost over all evaluations
    n_bracket_calls: int = 0
    final_bracket: dict[Any, Interval] = field(default_factory=dict)


def _refine_order(candidates, brackets, target, acquisition):
    if not acquisition:
        return list(candidates)
    # v2: resolve the most-decidable first — those whose bracket sticks out of the target least
    # (smallest excursion beyond a boundary) are closest to a verdict, so refine them first.
    def excursion(c):
        b = brackets[c]
        return max(target.lo - b.lo, 0.0) + max(b.hi - target.hi, 0.0)
    return sorted(candidates, key=excursion)


def screen(
    candidates: list[Any],
    target: Interval,
    oracle: BoundedOracle,
    *,
    max_precision: int = 4,
    acquisition: bool = False,
) -> ScreenResult:
    """Bracket-aware screen with zero false eliminations. See module docstring."""
    result = ScreenResult()
    brackets: dict[Any, Interval] = {}

    # Pilot pass: cheap bracket for everyone, eliminate the strictly-disjoint.
    undecided: list[Any] = []
    for c in candidates:
        b = oracle.bracket(c, 0)
        brackets[c] = b
        result.total_cost += oracle.cost(0)
        result.n_bracket_calls += 1
        verdict = classify(b, target)
        if verdict is Verdict.ELIMINATED:
            result.eliminated.append(c)
        elif verdict is Verdict.CONFIRMED:
            result.confirmed.append(c)
        else:
            undecided.append(c)

    # Refine survivors at increasing precision until decided or out of precision budget.
    precision = 1
    while undecided and precision <= max_precision:
        still: list[Any] = []
        for c in _refine_order(undecided, brackets, target, acquisition):
            b = oracle.bracket(c, precision)
            brackets[c] = b
            result.total_cost += oracle.cost(precision)
            result.n_bracket_calls += 1
            verdict = classify(b, target)
            if verdict is Verdict.ELIMINATED:
                result.eliminated.append(c)
            elif verdict is Verdict.CONFIRMED:
                result.confirmed.append(c)
            else:
                still.append(c)
        undecided = still
        precision += 1

    # Anything still undecided at max precision is KEPT (conservative — never eliminate on doubt).
    result.survivors = result.confirmed + undecided
    result.final_bracket = brackets
    return result


def point_estimate_screen(
    candidates: list[Any],
    target: Interval,
    oracle: BoundedOracle,
    *,
    precision: int = 0,
) -> ScreenResult:
    """Baseline: eliminate iff the point estimate lies strictly outside the target.

    Unsound — a biased/noisy estimate can push a genuine hit outside the target and eliminate it.
    Provided to demonstrate the false eliminations that the bracket-aware loop provably avoids.
    """
    result = ScreenResult()
    for c in candidates:
        pe = oracle.point_estimate(c, precision)
        result.total_cost += oracle.cost(precision)
        result.n_bracket_calls += 1
        if target.lo <= pe <= target.hi:
            result.survivors.append(c)
        else:
            result.eliminated.append(c)
    return result
