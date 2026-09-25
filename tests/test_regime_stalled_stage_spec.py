"""Gates for specs/SPEC_regime_stalled_stage.md (bead chem-4e9).

`truncation_regime` read discarded weights only, so a ramp stage that never converged -- the n=40
D=100 stall of SPEC_hchain_largen2 §11.2, 0.71 Ha above the D=200/400 energies but with small,
monotone weights -- was certified "truncation", and its fit landed 0.82 mHa BELOW its own D=400
energy. These gates make the regime read the energies too.

Pure: synthetic and vendored ``per_D`` triples, no block2 and no pyscf import.
"""
import csv
from pathlib import Path

import numpy as np
import pytest

from hybrid_quantum_solver import dmrg_reference as dr
from hybrid_quantum_solver.dmrg_reference import (
    DISCARD_WEIGHT_FLOOR,
    ENERGY_NOISE,
    STAGE_SLOPE_RATIO,
    UNDERSHOOT_FACTOR,
    extrapolate_ladder,
    truncation_regime,
)

REPO = Path(__file__).resolve().parents[1]
TABLE = REPO / "specs" / "hchain_tdl_localized_table.csv"

# §11.2 recorded: stalled D=100 at -20.924 Ha; D=400 at -21.6340616891 Ha (same digits as the
# 200/400/800 rerun). The weights were not recorded; these are the "small, monotone" profile that
# reproduces the recorded symptom (fit 0.82 mHa below E(D=400), stderr 8.0e-4) -- see G1.
E_D400 = -21.6340616891
STALLED = [(100, 1e-7, -20.924), (200, 2.3e-10, E_D400 + 1e-9), (400, 1e-12, E_D400)]


def _legacy_weight_regime(per_D):
    """The pre-chem-4e9 predicate (SPEC_extrap_regime), recomputed verbatim: weights only."""
    dws = np.array([p[1] for p in per_D], dtype=float)
    if len(dws) < 2:
        return "uncontrolled"
    if float(np.max(dws)) <= DISCARD_WEIGHT_FLOOR:
        return "converged"
    return "truncation" if np.all(np.diff(dws) <= 1e-12) else "uncontrolled"


def _table_ladders():
    for r in csv.DictReader(TABLE.open()):
        Ds = [int(x) for x in r["bond_dims"].split("/")]
        dws = [float(x) for x in r["dw_per_D"].split("/")]
        Es = [float(x) for x in r["e_per_D"].split("/")]
        yield r, list(zip(Ds, dws, Es))


# --- G1: the stalled ladder is uncontrolled (DEFINITION OF DONE) ---------------------------------

def test_G1_stalled_ladder_is_uncontrolled():
    # The fixture is faithful: weights alone call it "truncation", and the fit reproduces the
    # recorded symptom (0.82 mHa below E(D_max), stderr 8.0e-4).
    assert _legacy_weight_regime(STALLED) == "truncation"
    res = extrapolate_ladder(STALLED)
    assert abs((E_D400 - res.energy) - 0.82e-3) < 0.01e-3, E_D400 - res.energy
    assert abs(res.stderr - 8.0e-4) < 0.3e-4, res.stderr
    assert res.regime == "uncontrolled"
    assert truncation_regime(STALLED) == "uncontrolled"


@pytest.mark.parametrize("dws", [
    (1e-7, 2.3e-10, 1e-12),
    (1e-5, 2.3e-8, 1e-10),
    (1e-4, 1e-6, 1e-8),
    (3e-8, 1e-8, 1e-12),
    (1e-7, 1e-9, 1e-11),
])
def test_G1_stall_is_caught_across_weight_scales(dws):
    """Not tuned to one weight profile: any small monotone weights under a 0.71 Ha stage gap."""
    per_D = list(zip((100, 200, 400), dws, (-20.924, E_D400 + 1e-9, E_D400)))
    assert _legacy_weight_regime(per_D) == "truncation"
    assert truncation_regime(per_D) == "uncontrolled"


def test_G1_each_energy_check_catches_the_stall_on_its_own():
    """The stage-gap and undershoot checks are independent: disabling either still rejects."""
    _, dws, Es = dr._ladder_arrays(STALLED)
    assert dr._stalled_stage(dws, Es, ENERGY_NOISE, STAGE_SLOPE_RATIO)
    energy, _ = dr._linear_extrapolate(dws, Es)
    assert Es[-1] - energy > UNDERSHOOT_FACTOR * (Es[-2] - Es[-1]) + ENERGY_NOISE


# --- G2: non-variational ladders ------------------------------------------------------------------

def test_G2_energy_rising_with_D_is_uncontrolled():
    healthy = [(100, 1e-3, -1.000), (200, 3e-4, -1.001), (400, 1e-4, -1.0015)]
    assert truncation_regime(healthy) == "truncation"
    risen = [(100, 1e-3, -1.000), (200, 3e-4, -1.001), (400, 1e-4, -1.001 + 10 * ENERGY_NOISE)]
    assert truncation_regime(risen) == "uncontrolled"


def test_G2_rise_within_noise_is_tolerated():
    """A converged ladder that wobbles upward by float noise is still converged."""
    wobble = [(100, 1e-12, -1.0), (200, 1e-13, -1.0 - 1e-9), (400, 1e-14, -1.0 - 5e-10)]
    assert truncation_regime(wobble) == "converged"


# --- G3: weights below the floor cannot certify a stalled stage ----------------------------------

def test_G3_stall_below_the_floor_is_not_converged():
    per_D = [(100, 5e-9, -20.924), (200, 1e-10, E_D400 + 1e-9), (400, 1e-12, E_D400)]
    assert _legacy_weight_regime(per_D) == "converged"
    assert truncation_regime(per_D) == "uncontrolled"


# --- G4: no false positives on any recorded real ladder -------------------------------------------

def test_G4_every_vendored_ladder_keeps_its_recorded_regime_and_numbers():
    rows = 0
    for r, per_D in _table_ladders():
        rows += 1
        res = extrapolate_ladder(per_D)
        assert res.regime == r["regime"], (r["n"], res.regime)
        assert res.method == r["extrap_method"], r["n"]
        # e_per_D is printed to 1e-10, so the refit agrees to that, not bit-for-bit.
        assert abs(res.energy - float(r["e_dmrg_extrap"])) < 1e-9, r["n"]
    assert rows == 9


def test_G4_tolerances_have_headroom_on_recorded_truncation_ladders():
    """The constants are justified by the recorded ladders. Measured margins: slope ratio <= 1.2
    vs 1e3 (~800x); drop/gap_last <= 0.102 vs 10 (~98x). Asserted with >= 50x."""
    n_trunc = 0
    for r, per_D in _table_ladders():
        if r["regime"] != "truncation":
            continue
        n_trunc += 1
        _, dws, Es = dr._ladder_arrays(per_D)
        slopes = -np.diff(Es) / -np.diff(dws)
        assert slopes.max() / slopes.min() < STAGE_SLOPE_RATIO / 50, (r["n"], slopes)
        energy, _ = dr._linear_extrapolate(dws, Es)
        assert (Es[-1] - energy) / (Es[-2] - Es[-1]) < UNDERSHOOT_FACTOR / 50, r["n"]
    assert n_trunc == 3


# --- G5: numbers unchanged, only the label moves --------------------------------------------------

def test_G5_stalled_ladder_keeps_its_legacy_fit():
    """`method`, energy and stderr still follow the weights alone (SPEC_extrap_regime G3)."""
    res = extrapolate_ladder(STALLED)
    dws = np.array([p[1] for p in STALLED])
    Es = np.array([p[2] for p in STALLED])
    coef, cov = np.polyfit(dws, Es, 1, cov=True)
    assert res.method == "dweight"
    assert res.energy == float(coef[1])
    assert res.stderr == float(np.sqrt(cov[1, 1]))


def test_G5_tolerances_are_pinned():
    assert ENERGY_NOISE == 1e-6
    assert STAGE_SLOPE_RATIO == 1e3
    assert UNDERSHOOT_FACTOR == 10.0


# --- G6: the n=40 stall, reproduced (SPEC §8) -----------------------------------------------------
# `benchmark_hchain_tdl.py --localize` integrals, n=40, protocol="ramp", bond_dims 100/200/400,
# 4 threads, stack_mem 6 GB, unseeded -- the §11.2 run's settings. Recorded 2026-09-25, 947 s.
REPRODUCED_N40 = [
    (100, 7.37332779285622e-07, -20.924206023061657),
    (200, 1.6764842710361932e-09, -21.634061152184756),
    (400, 2.6902054150132614e-11, -21.634061689075224),
]


def test_G6_reproduced_n40_stall_is_rejected_by_the_regime_alone():
    # It is the §11.2 row: same stalled D=100, same D=400 digits, same undershoot and stderr.
    assert abs(REPRODUCED_N40[0][2] - (-20.924)) < 1e-3
    assert abs(REPRODUCED_N40[2][2] - E_D400) < 1e-9
    res = extrapolate_ladder(REPRODUCED_N40)
    assert abs((REPRODUCED_N40[2][2] - res.energy) - 0.82e-3) < 0.01e-3
    assert abs(res.stderr - 8.0e-4) < 0.1e-4
    # Weights alone certified it; the energies reject it. The numbers themselves are unchanged.
    assert _legacy_weight_regime(REPRODUCED_N40) == "truncation"
    assert res.method == "dweight"
    assert res.regime == "uncontrolled"
