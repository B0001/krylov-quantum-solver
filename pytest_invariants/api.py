"""The ``@invariant.*`` decorator API (pytest-invariants v1: decorate test functions).

A test decorated with an invariant *returns an observation* (the quantity under test — a
scalar, a tuple, a bracket object, whatever). Each stacked ``@invariant.<name>(...)``
decorator extracts what it needs from that observation via an optional ``key`` and asserts
its relation, raising :class:`InvariantViolation` (an ``AssertionError`` subclass, so pytest
reports it as a normal failure) with a message that stands on its own — expected relation,
observed values, and the *provenance* of the bound.

Why a single composed wrapper: decorators apply innermost-first, so the innermost one builds
the wrapper and every outer one just appends its spec onto the same wrapper's
``__invariants__`` list. The underlying test body therefore runs exactly once, its return is
the shared observation, and the wrapper returns ``None`` to pytest (so no
``PytestReturnNotNoneWarning``). The plugin reads ``__invariants__`` for the coverage report.
"""

from __future__ import annotations

import functools
import math
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any


class InvariantViolation(AssertionError):
    """Raised when an invariant fails. Subclasses ``AssertionError`` so pytest treats it as a
    plain test failure (and so ``pytest.raises(AssertionError)`` catches it in demo tests)."""


def _resolve(x: Any) -> Any:
    """A bound/reference may be a literal or a zero-arg callable computing it (e.g. an FCI
    lookup). Callables are resolved lazily here, at test time, not at decoration time."""
    return x() if callable(x) else x


def _fmt(x: Any, units: str) -> str:
    """Render a numeric value with optional units, keeping non-numerics (inf, tuples) legible."""
    suffix = f" {units}" if units else ""
    if isinstance(x, float):
        return f"{x:.12g}{suffix}"
    return f"{x}{suffix}"


@dataclass
class InvariantSpec:
    """One invariant attached to a test: its kind, the check to run, and reporting metadata."""

    name: str
    provenance: str
    check: Callable[[Any], None]  # (observation) -> None, raises InvariantViolation on failure
    extras: dict[str, Any] = field(default_factory=dict)


def _attach(func: Callable[..., Any], spec: InvariantSpec) -> Callable[..., Any]:
    """Attach ``spec`` to ``func``, composing onto an existing wrapper if one is already there."""
    if getattr(func, "_pytest_invariants_wrapper", False):
        func.__invariants__.append(spec)  # type: ignore[attr-defined]
        return func

    specs: list[InvariantSpec] = [spec]

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> None:
        observation = func(*args, **kwargs)
        # Check innermost-declared invariant first (source order, bottom-up).
        for s in wrapper.__invariants__:  # type: ignore[attr-defined]
            s.check(observation)
        # Return None so pytest does not warn about a non-None test return.
        return None

    wrapper._pytest_invariants_wrapper = True  # type: ignore[attr-defined]
    wrapper.__invariants__ = specs  # type: ignore[attr-defined]
    return wrapper


class _Invariant:
    """Namespace object exposing the decorators as ``invariant.lower_bound`` / ``.contains``."""

    def lower_bound(
        self,
        fn_or_value: Callable[[], float] | float,
        *,
        key: Callable[[Any], float] | None = None,
        tol: float = 0.0,
        provenance: str = "caller-supplied lower bound",
        units: str = "Ha",
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        """Assert the observed scalar is ``>= fn_or_value`` (within ``tol``).

        This is the variational-floor check (``certchem.floor_guard``) as a decorator: the test
        returns its estimate, ``key`` extracts the scalar to compare (identity by default), and
        ``fn_or_value`` is the floor — a literal or a zero-arg callable (e.g. an FCI reference).
        An estimate below the floor is the "-hundreds of Hartree" pathology this repo exists to
        catch. ``provenance`` names where the bound came from, so the failure explains itself.
        """

        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            def check(observation: Any) -> None:
                observed = key(observation) if key else observation
                bound = _resolve(fn_or_value)
                if observed < bound - tol:
                    tol_note = f" (tol {tol:g})" if tol else ""
                    raise InvariantViolation(
                        f"Invariant `lower_bound` violated in {func.__name__}:\n"
                        f"  expected:  observed >= bound{tol_note}\n"
                        f"  observed:  {_fmt(observed, units)}\n"
                        f"  bound:     {_fmt(bound, units)}  "
                        f"(provenance: {provenance})\n"
                        f"  margin:    observed - bound = "
                        f"{_fmt(observed - bound, units)}  (must be >= 0)"
                    )

            spec = InvariantSpec(
                name="lower_bound",
                provenance=provenance,
                check=check,
                extras={"tol": tol},
            )
            return _attach(func, spec)

        return decorator

    def contains(
        self,
        reference: Callable[[], float] | float,
        *,
        key: Callable[[Any], tuple[float, float]] | None = None,
        provenance: str = "caller-supplied reference",
        units: str = "Ha",
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        """Assert a ``[lower, upper]`` interval brackets ``reference``.

        The certified-bracket-contains-FCI check as a decorator: the test returns its bracket,
        ``key`` maps it to ``(lower, upper)`` (identity expects a 2-tuple), and ``reference`` is
        the value that must fall inside — a literal or zero-arg callable (e.g. an FCI lookup).
        ``provenance`` names the reference so a failure reads without opening this file.
        """

        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            def check(observation: Any) -> None:
                lower, upper = key(observation) if key else observation
                ref = _resolve(reference)
                if not (lower <= ref <= upper):
                    if math.isfinite(ref) and math.isfinite(lower) and ref < lower:
                        rel = f"reference is BELOW the interval by {_fmt(lower - ref, units)}"
                    elif math.isfinite(ref) and math.isfinite(upper) and ref > upper:
                        rel = f"reference is ABOVE the interval by {_fmt(ref - upper, units)}"
                    else:
                        rel = "reference lies outside the interval"
                    raise InvariantViolation(
                        f"Invariant `contains` violated in {func.__name__}:\n"
                        f"  expected:  lower <= reference <= upper\n"
                        f"  reference: {_fmt(ref, units)}  (provenance: {provenance})\n"
                        f"  interval:  [{_fmt(lower, units)}, {_fmt(upper, units)}]\n"
                        f"  {rel}"
                    )

            spec = InvariantSpec(
                name="contains",
                provenance=provenance,
                check=check,
            )
            return _attach(func, spec)

        return decorator


#: The public decorator namespace: ``from pytest_invariants import invariant``.
invariant = _Invariant()
