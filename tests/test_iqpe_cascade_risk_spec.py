"""
Acceptance gates G1-G4 for specs/SPEC_iqpe_cascade_risk.md (iterative QPE -- a fat-tailed bit-flip
cascade the median error hides).

Deliberately no new library code: the checks live entirely here, reusing `iterative_qpe.py`'s
existing `iqpe_ground_energy` unmodified -- a genuine external characterization of behavior the
module already exhibits stochastically (per-run) but never summarizes across seeds.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from iterative_qpe import iqpe_ground_energy

N_SEEDS = 40


@pytest.fixture(scope="module")
def h2_system():
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, 2, 2)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), 2)
    return h1, eri, float(e_core), (1, 1), float(cas.e_tot)


def _errors(h2_system, n_bits, shots_per_bit, n_seeds=N_SEEDS):
    h1, eri, e_core, nelec, casci = h2_system
    errs = [
        abs(iqpe_ground_energy(h1, eri, e_core, nelec, 2, n_bits=n_bits,
                               shots_per_bit=shots_per_bit, seed=s)[0] - casci)
        for s in range(n_seeds)
    ]
    return np.array(errs)


def test_G1_median_precision_monotonically_non_increasing_with_bits(h2_system):
    """At the module's own default shots_per_bit=15, the median error over 30 seeds is
    non-increasing across n_bits -- the docstring's directional claim holds, checked not assumed."""
    medians = [np.median(_errors(h2_system, n, 15, n_seeds=30))
              for n in (4, 6, 8, 10, 12, 14, 16, 18, 20)]
    for m_next, m_prev in zip(medians[1:], medians[:-1]):
        assert m_next <= m_prev + 1e-15, medians


def test_G2_fat_tail_at_low_shots_per_bit_vanishes_at_the_default(h2_system):
    """THE FINDING / definition of done: at n_bits=8, shots_per_bit=3, the worst-case (max) error
    over 40 seeds is more than 5x the median -- a real cascade from a flipped low-significance bit
    corrupting the feedback for every subsequent bit. At the SAME n_bits with the module's own
    default shots_per_bit=15, the ratio collapses to exactly 1 (fully deterministic) -- the default
    already avoids the pathology, but a naive shot-cutting choice reintroduces it."""
    errs_low = _errors(h2_system, 8, 3)
    ratio_low = np.max(errs_low) / np.median(errs_low)
    assert ratio_low > 5.0, ratio_low

    errs_default = _errors(h2_system, 8, 15)
    assert np.max(errs_default) == np.median(errs_default), (
        np.max(errs_default), np.median(errs_default)
    )


def test_G3_the_cascade_has_a_measured_shots_per_bit_threshold(h2_system):
    """The tail is not open-ended: at n_bits=12, shots_per_bit=1, the cascade is present
    (max/median > 5); at shots_per_bit=3 and above, it is gone (max == median, deterministic)."""
    errs_spb1 = _errors(h2_system, 12, 1)
    assert np.max(errs_spb1) / np.median(errs_spb1) > 5.0, errs_spb1

    for spb in (3, 5, 7, 9):
        errs = _errors(h2_system, 12, spb)
        assert np.max(errs) == np.median(errs), (spb, np.max(errs), np.median(errs))


def test_G4_smallest_bit_counts_never_reach_chemical_accuracy(h2_system):
    """Boundary, recorded not smoothed over: even at the module's own default shots_per_bit=15,
    n_bits=8 stays above chemical accuracy (10 mHa) on EVERY tested seed -- the __main__ demo's
    own printed n_bits=4,6,8 rows are honestly still in the gross-error regime."""
    errs = _errors(h2_system, 8, 15)
    assert np.all(errs > 1e-2), errs
