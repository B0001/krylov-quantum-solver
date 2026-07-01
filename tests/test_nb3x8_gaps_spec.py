"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_gaps.md (exact Nb3X8 cluster gaps vs Hubbard-I).

Test-first: ``nb3x8_gaps`` does not exist yet, so this file is RED until the spec is implemented. The
Nb3X8 bilayer cluster (a generalized Hubbard dimer, from the arXiv:2501.10320 cRPA parameters) is
exactly diagonalizable; we compute the exact charge gap and the Hubbard-I gap and quantify where the
Hubbard-I approximation the source paper uses breaks down. Reference: exact diagonalization + the
t->0 atomic limit (both -> U0).

PySCF FCI / NumPy only (no block2); `make gates` runs it in its own process.
"""
import numpy as np

from nb3x8_gaps import NB3X8_LT_BULK, exact_charge_gap, hubbard_i_gap


def test_G1_atomic_limit_validation():
    """t -> 0: both the exact and Hubbard-I gaps reduce to the atomic Mott gap U0."""
    for p in NB3X8_LT_BULK.values():
        assert abs(exact_charge_gap(p["U0"], 1e-6, p["Us"]) - p["U0"]) < 1e-3
        assert abs(hubbard_i_gap(p["U0"], 1e-6, p["Us"]) - p["U0"]) < 1e-3


def test_G2_exact_gaps_are_the_new_numbers():
    """The exact charge gaps (meV), positive/insulating for every compound."""
    expected = {"Nb3I8": 842.4, "Nb3Br8": 1086.0, "Nb3Cl8": 1311.8, "Nb3F8": 2580.8}
    for name, gap in expected.items():
        computed = exact_charge_gap(**NB3X8_LT_BULK[name])
        assert abs(computed - gap) < 1.0, (name, computed, gap)
        assert computed > 0.0


def test_G3_hubbard_i_error_grows_toward_weak_coupling():
    """DEFINITION OF DONE: Hubbard-I is exact for strong correlation but underestimates the
    weakly-correlated Nb3I8 gap by >20%; the error grows monotonically as U0/|t| falls."""
    order = ["Nb3F8", "Nb3Cl8", "Nb3Br8", "Nb3I8"]                  # decreasing U0/|t|
    rel_err = []
    for name in order:
        p = NB3X8_LT_BULK[name]
        ge, gh = exact_charge_gap(**p), hubbard_i_gap(**p)
        rel_err.append(abs(gh - ge) / ge)
    assert all(rel_err[i] <= rel_err[i + 1] for i in range(len(rel_err) - 1)), rel_err   # monotone
    assert rel_err[0] < 0.01 and rel_err[1] < 0.01, rel_err        # F, Cl ~ exact
    assert rel_err[-1] > 0.20, rel_err                             # Nb3I8 > 20% (measured 29%)
    pI = NB3X8_LT_BULK["Nb3I8"]
    assert hubbard_i_gap(**pI) < exact_charge_gap(**pI)            # underestimate


def test_G4_error_correlates_with_correlation_strength():
    """The Hubbard-I error is ordered by U0/|t| (weaker correlation -> larger error)."""
    ratios, errs = [], []
    for p in NB3X8_LT_BULK.values():
        ratios.append(p["U0"] / abs(p["t"]))
        ge, gh = exact_charge_gap(**p), hubbard_i_gap(**p)
        errs.append(abs(gh - ge) / ge)
    # rank correlation between U0/|t| and error must be perfectly negative (Spearman = -1)
    order_by_ratio = np.argsort(ratios)
    err_ranks = np.argsort(np.argsort(np.array(errs)[order_by_ratio]))
    assert list(err_ranks) == sorted(err_ranks, reverse=True), (ratios, errs)
