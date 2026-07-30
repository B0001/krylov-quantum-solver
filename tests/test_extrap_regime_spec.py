"""Gates for specs/SPEC_extrap_regime.md.

The DMRG extrapolation `method` label conflates two opposite states -- *converged* (nothing left to
extrapolate) and *uncontrolled* (non-monotone truncation) -- because ``np.allclose(dws, 0.0)``
reduces to ``max(dws) <= 1e-8`` (numpy's default atol). Four gate files assert
``method == "dweight"`` as a quality guard, so a schedule converged past that floor fails a gate
written to demand quality.

These gates are deliberately **pure**: synthetic ``per_D`` triples, no block2 and no pyscf import.
They therefore run in environments where the four DMRG gates skip -- which is the point, since that
is exactly where a broken regime predicate would otherwise ship green.
"""
import re
from pathlib import Path

import numpy as np
import pytest

from hybrid_quantum_solver.dmrg_reference import (
    DISCARD_WEIGHT_FLOOR,
    REGIMES,
    truncation_regime,
)

REPO = Path(__file__).resolve().parents[1]

# The four gate files that assert on the extrapolation regime label.
GATE_FILES = (
    "tests/test_hchain_tdl_spec.py",
    "tests/test_hchain_largen2_spec.py",
    "tests/test_nbn_dmrg_reference_spec.py",
    "tests/test_singleramp_spec.py",
)


def _per_D(dws, Ds=None, Es=None):
    """Build (D, discarded_weight, E) triples; only the weights drive the regime."""
    n = len(dws)
    Ds = Ds if Ds is not None else [100 * 2**i for i in range(n)]
    Es = Es if Es is not None else [-1.0 - 1e-3 * i for i in range(n)]
    return list(zip(Ds, dws, Es))


def _legacy_method(per_D):
    """The pre-change expression from dmrg_reference.py:193-198, recomputed verbatim.

    This is the oracle for G3: the refactor must not move a single label.
    """
    if not per_D:
        return "invD"
    dws = np.array([p[1] for p in per_D], dtype=float)
    usable = (len(per_D) >= 2 and not np.allclose(dws, 0.0)
              and np.all(np.diff(dws) <= 1e-12))
    return "dweight" if usable else "invD"


# --- G1: the three regimes actually separate (DEFINITION OF DONE) --------------------------------

@pytest.mark.parametrize("name,dws,expected", [
    # The recorded NbN headline profile, SPEC_nbn_dmrg_reference.md:55 -- "weights ~1e-9-1e-13:
    # nothing left to extrapolate". Today this is demoted and indistinguishable from failure.
    ("nbn_headline_converged", [1e-9, 1e-11, 1e-13], "converged"),
    # A healthy ladder: weights resolvable, monotone non-increasing.
    ("healthy_truncation", [1e-3, 3e-4, 1e-4], "truncation"),
    # The failure that killed SPEC_hchain_largen.md: above the floor and non-monotone.
    ("non_monotone_uncontrolled", [1e-4, 5e-4, 2e-4], "uncontrolled"),
])
def test_G1_three_regimes_separate(name, dws, expected):
    assert truncation_regime(_per_D(dws)) == expected, name


def test_G1_degenerate_ladder_is_uncontrolled():
    """Fewer than two points cannot support any extrapolation."""
    assert truncation_regime(_per_D([1e-4])) == "uncontrolled"
    assert truncation_regime([]) == "uncontrolled"


def test_G1_regimes_are_exactly_the_declared_set():
    assert set(REGIMES) == {"converged", "truncation", "uncontrolled"}


# --- G2: the floor boundary is pinned, not inherited from a library default ----------------------

def test_G2_floor_boundary_is_explicit():
    """The physical threshold is asserted. It must not be an implicit numpy atol again."""
    assert DISCARD_WEIGHT_FLOOR == 1e-8
    # Just above the floor -> still extrapolatable.
    assert truncation_regime(_per_D([2e-8, 1.5e-8, 1.1e-8])) == "truncation"
    # Just below -> nothing left to extrapolate.
    assert truncation_regime(_per_D([2e-9, 1.5e-9, 1.1e-9])) == "converged"


def test_G2_floor_is_a_max_not_a_last_element():
    """R3: using dws[-1] would label a normally-converging ladder 'converged', abandon the
    discarded-weight fit and silently change extrapolated energies."""
    assert truncation_regime(_per_D([1e-3, 1e-5, 1e-9])) == "truncation"


# --- G3: bit-for-bit `method` compatibility -- the no-numerics-changed oracle --------------------

@pytest.mark.parametrize("dws", [
    [1e-9, 1e-11, 1e-13],
    [1e-3, 3e-4, 1e-4],
    [1e-4, 5e-4, 2e-4],
    [2e-8, 1.5e-8, 1.1e-8],
    [2e-9, 1.5e-9, 1.1e-9],
    [1e-3, 1e-5, 1e-9],
    [1e-12, 3e-12, 1e-12],
    [1e-4],
    [0.0, 0.0],
])
def test_G3_method_is_bit_for_bit_unchanged(dws):
    """Every regime maps back to the legacy two-valued label with no drift.

    One mismatch kills the claim that this change is label-only.
    """
    per_D = _per_D(dws)
    regime = truncation_regime(per_D)
    new_method = "dweight" if regime == "truncation" else "invD"
    assert new_method == _legacy_method(per_D), (dws, regime)


# --- G4: ordering -- converged beats non-monotone ------------------------------------------------

def test_G4_converged_wins_over_float_noise_non_monotonicity():
    """Below the floor AND non-monotone by float noise: this is a converged run, not a broken one.

    Kills any predicate that tests monotonicity before the floor.
    """
    assert truncation_regime(_per_D([1e-12, 3e-12, 1e-12])) == "converged"


# --- G5: no gate compares against a string outside REGIMES ---------------------------------------

def test_G5_gate_files_only_compare_declared_regimes():
    """Guard against shipping a green, dead assertion in a test that skips locally.

    `regime != "invD"` is a plausible slip -- that value moved to `method` -- and it would pass
    forever. So would a typo. Any literal compared against `.regime` must be a declared regime.
    """
    pattern = re.compile(r"""\.regime\s*[!=]=\s*["']([^"']+)["']""")
    checked = 0
    for rel in GATE_FILES:
        path = REPO / rel
        if not path.exists():
            continue
        for literal in pattern.findall(path.read_text()):
            checked += 1
            assert literal in REGIMES, f"{rel}: .regime compared against {literal!r}"
    assert checked > 0, "no .regime assertions found in the gate files -- migration incomplete"
