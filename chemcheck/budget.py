"""Mode A error-budget model (tasks 3–5). Pure functions, zero solver deps.

**v1 crudeness, stated loudly (SPEC_chemcheck.md caveats):** depolarizing-only, a single
reference Trotter step of two-qubit gates scaled by a published routing multiplier, and a
"≤ 1 expected two-qubit error per circuit" coherence heuristic as the PASS threshold. This is
uncalibrated — M2's noisy-sim sweep is what would validate it. Headroom is order-of-magnitude.
"""

from __future__ import annotations

#: Routing/SWAP overhead multiplier per connectivity class (v1 literature lookup, not measured).
#: Ordering that MUST hold: all_to_all (1.0) < grid < heavy_hex < linear. Sources: typical SWAP
#: blow-up for mapping dense fermionic circuits onto sparse hardware graphs (e.g. heavy-hex needs
#: ~3–4× the two-qubit gates of an all-to-all layout; a 1-D line is worse still).
ROUTING_OVERHEAD: dict[str, float] = {
    "all_to_all": 1.0,
    "grid": 2.5,
    "heavy_hex": 3.5,
    "linear": 5.0,
}

#: PASS coherence heuristic: a circuit is scored PASS when its expected number of two-qubit
#: errors stays at or below this. "One error per circuit" is the standard order-of-magnitude line.
MAX_EXPECTED_TWO_QUBIT_ERRORS = 1.0


def routing_overhead(connectivity: str) -> float:
    """Two-qubit-gate overhead multiplier for a connectivity class."""
    try:
        return ROUTING_OVERHEAD[connectivity]
    except KeyError:
        raise ValueError(
            f"unknown connectivity {connectivity!r}; known: {sorted(ROUTING_OVERHEAD)} "
            "(custom graphs need an explicit overhead)"
        ) from None


def expected_total_error(
    two_qubit_count: int, overhead: float, two_qubit_error: float
) -> float:
    """Depolarizing probability that the routed circuit suffers ≥1 two-qubit error.

    ``1 - (1 - p)^(N · overhead)`` — bounded in [0, 1). Pure; unit-tested against hand cases.
    """
    routed_gates = two_qubit_count * overhead
    return 1.0 - (1.0 - two_qubit_error) ** routed_gates


def required_two_qubit_error(
    two_qubit_count: int, overhead: float, max_expected_errors: float = MAX_EXPECTED_TWO_QUBIT_ERRORS
) -> float:
    """Two-qubit error at which the routed circuit hits exactly the PASS threshold.

    From ``routed_gates · p = max_expected_errors`` (linear small-error regime, the inverse of
    the coherence heuristic): ``p_required = max_expected_errors / (N · overhead)``.
    """
    routed_gates = two_qubit_count * overhead
    if routed_gates <= 0:
        raise ValueError("two_qubit_count · overhead must be positive")
    return max_expected_errors / routed_gates


def headroom_factor(current_two_qubit_error: float, required_two_qubit_error: float) -> float:
    """current / required. ``> 1`` ⇒ FAIL (too error-prone); ``1.0`` ⇒ exactly at threshold.

    Equivalently the expected number of two-qubit errors in the routed circuit. Monotone
    decreasing in device quality: a smaller current error → smaller headroom.
    """
    if required_two_qubit_error <= 0:
        raise ValueError("required_two_qubit_error must be positive")
    return current_two_qubit_error / required_two_qubit_error
