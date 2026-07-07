"""
Acceptance gates G1-G4 for specs/SPEC_certified_thermochem.md.

Claim: a relative energy Delta = E(B) - E(A) (reaction / dissociation / stretch) gets a certified
interval from Krylov data alone (no FCI), by composing the Temple/Ritz brackets at each geometry:
Delta in [tau_B - rho_A, rho_B - tau_A]. The exact in-basis FCI relative energy lies inside at every
depth (zero escapes); the interval closes with depth; and it is dominated by the strongly-correlated
endpoint, whose vacuous Temple lower bound at intermediate depth makes Delta one-sided there (the
temple premise regime).

Exact statevector, sector-restricted. PySCF/qiskit, no block2; `make gates` runs it in its own
process.
"""
import math

import pytest

from certified_thermochem import (
    HARTREE_TO_EV,
    certified_relative_energy_ladder,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian


def _chain(r):
    return f"H 0 0 0; H 0 0 {r}; H 0 0 {2 * r}; H 0 0 {3 * r}"


_DIMS = (6, 8, 10, 12, 16, 20)


@pytest.fixture(scope="module")
def h4_stretch():
    mh_eq = build_molecular_hamiltonian(atom=_chain(0.9))
    mh_st = build_molecular_hamiltonian(atom=_chain(2.3))
    de_fci = mh_st.ground_state_energy() - mh_eq.ground_state_energy()   # Ha
    ladder = certified_relative_energy_ladder(mh_eq, mh_st, _DIMS)
    return de_fci, ladder


def test_G1_zero_escapes(h4_stretch):
    """The exact FCI relative energy lies inside [Delta_lo, Delta_hi] at EVERY depth (treating a
    vacuous side as +/-inf -- still a valid one-sided certificate). One escape kills the claim."""
    de_fci, ladder = h4_stretch
    for rb in ladder:
        assert rb.delta_lower - 1e-9 <= de_fci <= rb.delta_upper + 1e-9, (
            rb.m, rb.delta_lower, de_fci, rb.delta_upper)


def test_G2_two_sided_interval_closes(h4_stretch):
    """DEFINITION OF DONE: at sufficient depth the certified interval is finite, two-sided, tight
    (< 0.05 eV) and contains the exact relative energy (~8.23 eV) -- a genuinely useful certified
    reaction energy without FCI."""
    de_fci, ladder = h4_stretch
    deep = ladder[-1]                                                     # M=20
    assert math.isfinite(deep.width), deep.width
    assert deep.width * HARTREE_TO_EV < 0.05, deep.width * HARTREE_TO_EV
    assert deep.delta_lower - 1e-9 <= de_fci <= deep.delta_upper + 1e-9
    assert abs(de_fci * HARTREE_TO_EV - 8.2255) < 1e-2                    # in-basis FCI sanity


def test_G3_correlated_endpoint_dominates(h4_stretch):
    """THE FINDING: the certified error bar is dominated by the strongly-correlated (stretched)
    endpoint -- the equilibrium bracket closes far faster. At M=6 (both finite) the equilibrium
    width is < 0.1x the stretched width, so the uncertainty is localized at the hard geometry."""
    _, ladder = h4_stretch
    rb = ladder[0]                                                       # M=6
    assert math.isfinite(rb.width_a) and math.isfinite(rb.width_b)
    assert rb.width_a < 0.1 * rb.width_b, (rb.width_a, rb.width_b)


def test_G4_inherits_temple_premise_but_upper_always_holds(h4_stretch):
    """The temple-premise boundary: at intermediate depth the stretched endpoint's Temple lower
    bound is vacuous, so Delta_lower is -inf (a one-sided certificate) -- yet the UPPER certificate
    (Delta_upper, from the Ritz values + the easy endpoint's Temple) is finite and holds at EVERY
    depth. The correlated endpoint sets the premise regime, exactly as certified_gaps charts."""
    de_fci, ladder = h4_stretch
    assert any(rb.delta_lower == -math.inf for rb in ladder), "expected a vacuous-lower depth"
    for rb in ladder:
        assert math.isfinite(rb.delta_upper), rb.m                       # upper always finite
        assert de_fci <= rb.delta_upper + 1e-9, (rb.m, de_fci, rb.delta_upper)
