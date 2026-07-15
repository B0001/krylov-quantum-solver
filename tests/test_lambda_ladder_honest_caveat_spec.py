"""
Acceptance gates G1-G4 for specs/SPEC_lambda_ladder_honest_caveat.md (lambda_ladder -- the
docstring's honest caveat, and an unmonotonic accuracy trend it doesn't mention).

Deliberately no new library code: every gate calls `lambda_ladder.py`'s and `df_factorization.py`'s
existing functions directly (`lambda_and_terms`, `fit_thc`, `fci_energy_error`,
`double_factorize`, `reconstruct_eri`) -- the same building blocks `lambda_ladder()`'s print loop
already uses, captured as return values instead of parsed from printed output.
"""
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from df_factorization import double_factorize, reconstruct_eri
from lambda_ladder import fci_energy_error, fit_thc, lambda_and_terms


@pytest.fixture(scope="module")
def n2_system():
    mol = gto.M(atom="N 0 0 0; N 0 0 1.10", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    norb, ne = 3, 4
    cas = mcscf.CASCI(mf, norb, ne)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    nelec = (ne // 2, ne // 2)
    return h1, eri, norb, nelec, e_core, float(cas.e_tot)


@pytest.fixture(scope="module")
def df_sweep(n2_system):
    h1, eri, norb, nelec, e_core, casci = n2_system
    _, _, full_rank = double_factorize(eri, norb)
    rows = []
    for R in range(1, full_rank + 1):
        leaves, _, _ = double_factorize(eri, norb, rank=R)
        eriR = reconstruct_eri(leaves, norb)
        lam, terms = lambda_and_terms(h1, eriR, norb)
        err = fci_energy_error(h1, eriR, norb, nelec, e_core, casci)
        rows.append({"rank": R, "lambda": lam, "terms": terms, "err_mHa": err})
    return full_rank, rows


@pytest.fixture(scope="module")
def thc_sweep(n2_system):
    h1, eri, norb, nelec, e_core, casci = n2_system
    rows = []
    for M in range(2, 7):
        eriT = fit_thc(eri, norb, M)
        lam, terms = lambda_and_terms(h1, eriT, norb)
        err = fci_energy_error(h1, eriT, norb, nelec, e_core, casci)
        rows.append({"M": M, "lambda": lam, "terms": terms, "err_mHa": err})
    return rows


def test_G1_full_rank_reconstruction_is_exact_for_both_methods(n2_system, df_sweep, thc_sweep):
    """DF at its own reported full_rank and THC at M=6 (matching full_rank on this system) both
    reproduce the naive lambda and give zero FCI error."""
    h1, eri, norb, _nelec, _e_core, _casci = n2_system
    lam_naive, _terms_naive = lambda_and_terms(h1, eri, norb)

    full_rank, df_rows = df_sweep
    df_full = next(r for r in df_rows if r["rank"] == full_rank)
    assert abs(df_full["lambda"] - lam_naive) < 1e-6, (df_full, lam_naive)
    assert df_full["err_mHa"] < 1e-6, df_full

    thc_full = next(r for r in thc_sweep if r["M"] == full_rank)
    assert abs(thc_full["lambda"] - lam_naive) < 1e-6, (thc_full, lam_naive)
    assert thc_full["err_mHa"] < 1e-6, thc_full


def test_G2_the_honest_caveat_is_true_not_just_prose(df_sweep, thc_sweep):
    """THE FINDING / definition of done: at the cheapest tested rank, THC's lambda exceeds DF's
    by more than 20%, and THC's term count exceeds DF's by at least 3x -- "comparable to or denser
    than DF" at small CAS is a checked inequality, not an assertion."""
    _full_rank, df_rows = df_sweep
    df_cheapest = df_rows[0]  # R=1
    thc_cheapest = thc_sweep[0]  # M=2 (the module's own lowest thc_ranks value)

    assert thc_cheapest["lambda"] > 1.2 * df_cheapest["lambda"], (thc_cheapest, df_cheapest)
    assert thc_cheapest["terms"] >= 3 * df_cheapest["terms"], (thc_cheapest, df_cheapest)


def test_G3_df_rank_truncation_accuracy_is_not_monotonic(df_sweep):
    """Boundary, recorded not smoothed over: somewhere in the DF rank sweep, a HIGHER rank gives a
    WORSE FCI error than a lower rank -- the docstring never mentions this, and G1's endpoint
    behavior alone would suggest a smooth "more factors = better" story the data doesn't support."""
    _full_rank, df_rows = df_sweep
    errs = [r["err_mHa"] for r in df_rows]
    non_monotonic = any(errs[i + 1] > errs[i] + 1e-9 for i in range(len(errs) - 1))
    assert non_monotonic, errs


def test_G4_thc_worst_accuracy_is_at_the_cheapest_rank(thc_sweep):
    """Sanity: THC's error at the lowest tested M is the worst of the swept range -- confirms
    THC's fit isn't accidentally non-monotonic in a way that would confound G2's "cheapest rank"
    comparison (i.e. "cheapest" and "least accurate" align for THC here)."""
    errs = [r["err_mHa"] for r in thc_sweep]
    assert errs[0] == max(errs), errs
    assert thc_sweep[-1]["err_mHa"] < 1e-6, thc_sweep[-1]
