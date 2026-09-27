"""
Acceptance gates G1-G5 for specs/SPEC_lambda_ladder_honest_caveat.md (lambda_ladder -- the
docstring's honest caveat, plus the orbital-gauge sensitivity found while gating it).

Deliberately no new library code: every gate calls `lambda_ladder.py`'s and `df_factorization.py`'s
existing functions directly (`lambda_and_terms`, `fit_thc`, `fci_energy_error`,
`double_factorize`, `reconstruct_eri`) -- the same building blocks `lambda_ladder()`'s print loop
already uses, captured as return values instead of parsed from printed output.

G3 history (chem-1yr): this gate originally asserted DF's rank-truncation error sequence was
NON-monotonic on N2 CAS(3,4). That was true under the *unpinned* orbital gauge (plain `gto.M`
with no symmetry) but was itself a gauge artifact, not a property of double factorization: N2's
HOMO pi pair is exactly degenerate, CASCI(norb=3, ne=4) freezes only one member of that pair into
"core", and which arbitrary rotation of the degenerate pair RHF happened to return determined
whether the resulting (frozen-orbital-dependent) active-space integrals gave a monotonic or
non-monotonic error sequence -- a ~50/50 split confirmed by explicitly sweeping the rotation angle
(G5 below). Pinning the gauge with `symmetry="D2h"` makes the sequence monotonic every time; G3 is
revised below to assert that, and G5 is added to keep the retracted claim's evidence runnable
rather than just prose. See SPEC_lambda_ladder_honest_caveat.md Sec 8 R2.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from df_factorization import double_factorize, reconstruct_eri
from lambda_ladder import fci_energy_error, fit_thc, lambda_and_terms


@pytest.fixture(scope="module")
def n2_system():
    # symmetry='D2h' pins the gauge: N2/sto-3g's HOMO pi pair (indices 4,5) and LUMO
    # pi* pair (7,8) are EXACTLY degenerate (D_inf,h), and CASCI(norb=3, ne=4)'s default
    # active-space selection (freeze the 5 lowest orbitals by energy) splits the HOMO
    # pair across the core/active boundary -- orbital 4 frozen, orbital 5 active. Without
    # symmetry, RHF returns an arbitrary (run-to-run varying) orthonormal basis for that
    # degenerate pair, and because freezing breaks the pair's rotational invariance, which
    # specific rotation lands in "core" vs "active" is a real physical choice, not a null
    # gauge move -- it changes CASCI's energy and every downstream number (see chem-1yr,
    # SPEC_lambda_ladder_honest_caveat.md Sec 8 R2). symmetry='D2h' makes PySCF solve each
    # irrep block independently, giving a reproducible, canonical split.
    mol = gto.M(atom="N 0 0 0; N 0 0 1.10", basis="sto-3g", symmetry="D2h")
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


def test_G3_df_rank_truncation_accuracy_is_monotonic_once_gauge_is_pinned(df_sweep):
    """REVISED (chem-1yr): under the canonical symmetry='D2h' gauge, the DF rank sweep's FCI error
    is monotonically non-increasing in rank -- the original "non-monotonic" claim was an artifact
    of an unpinned, run-to-run-arbitrary orbital gauge in N2's degenerate HOMO pi pair (see G5),
    not a property of double-factorization rank truncation. Retracted, not smoothed over: this
    test now encodes the corrected finding, and G5 encodes the gauge-freedom evidence that forced
    the retraction, so the retraction itself stays falsifiable rather than becoming prose."""
    _full_rank, df_rows = df_sweep
    errs = [r["err_mHa"] for r in df_rows]
    non_monotonic = [(i, errs[i], errs[i + 1]) for i in range(len(errs) - 1)
                     if errs[i + 1] > errs[i] + 1e-9]
    assert not non_monotonic, (non_monotonic, errs)


def test_G4_thc_worst_accuracy_is_at_the_cheapest_rank(thc_sweep):
    """Sanity: THC's error at the lowest tested M is the worst of the swept range -- confirms
    THC's fit isn't accidentally non-monotonic in a way that would confound G2's "cheapest rank"
    comparison (i.e. "cheapest" and "least accurate" align for THC here)."""
    errs = [r["err_mHa"] for r in thc_sweep]
    assert errs[0] == max(errs), errs
    assert thc_sweep[-1]["err_mHa"] < 1e-6, thc_sweep[-1]


def test_G5_the_retracted_G3_claim_was_an_unpinned_orbital_gauge_artifact():
    """THE ROOT CAUSE, kept runnable (chem-1yr): N2/sto-3g's RHF HOMO pi pair (mo indices 4, 5) is
    exactly degenerate, so ANY orthonormal rotation of that pair is an equally valid RHF solution
    with the same total energy. CASCI(norb=3, ne=4) freezes only orbital 4 into "core" and keeps
    orbital 5 active -- freezing breaks the pair's rotational symmetry, so the specific rotation
    RHF happens to return (arbitrary/run-dependent without a symmetry constraint) is a physically
    consequential choice for the resulting active-space integrals, not a null gauge move.

    This test builds the unpinned (no-symmetry) RHF once, then explicitly rotates the degenerate
    pair through a swept angle and reruns CASCI + the DF rank sweep at each angle -- a controlled,
    deterministic demonstration (not incidental flakiness) that BOTH a monotonic and a
    non-monotonic error sequence are reachable by equally-valid, equal-energy gauges of the SAME
    physical system. That mix is why G3's original "non-monotonic" claim was retracted rather than
    just reported as flaky, and why the fixture used elsewhere in this file pins the gauge with
    symmetry='D2h'."""
    mol = gto.M(atom="N 0 0 0; N 0 0 1.10", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    norb, ne = 3, 4
    nelec = (ne // 2, ne // 2)
    base_coeff = mf.mo_coeff.copy()

    # Degenerate pair straddling the core(idx 4)/active(idx 5) boundary at CAS(3,4).
    energy_gap = abs(mf.mo_energy[5] - mf.mo_energy[4])
    assert energy_gap < 1e-6, ("mo indices 4,5 are not degenerate on this system/geometry -- "
                               "the gauge argument below no longer applies", mf.mo_energy)

    outcomes = []
    for theta in np.linspace(0, np.pi, 13):
        c, s = np.cos(theta), np.sin(theta)
        coeff = base_coeff.copy()
        v4, v5 = base_coeff[:, 4].copy(), base_coeff[:, 5].copy()
        coeff[:, 4] = c * v4 - s * v5
        coeff[:, 5] = s * v4 + c * v5

        cas = mcscf.CASCI(mf, norb, ne)
        cas.verbose = 0
        cas.kernel(mo_coeff=coeff)
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        casci = float(cas.e_tot)

        _, _, full_rank = double_factorize(eri, norb)
        errs = []
        for R in range(1, full_rank + 1):
            leaves, _, _ = double_factorize(eri, norb, rank=R)
            eriR = reconstruct_eri(leaves, norb)
            errs.append(fci_energy_error(h1, eriR, norb, nelec, e_core, casci))
        non_monotonic = any(errs[i + 1] > errs[i] + 1e-9 for i in range(len(errs) - 1))
        outcomes.append(non_monotonic)

    # The rotation is a pure gauge freedom (same RHF solution, same energy at theta=0 vs any
    # theta), yet it flips the DF rank sweep between monotonic and non-monotonic -- proof the
    # original G3 result depended on WHICH gauge RHF happened to return, not on double
    # factorization itself. Exact counts vary run to run (the sweep's phase depends on the
    # unpinned baseline), so only the qualitative mix is asserted.
    assert any(outcomes), "expected at least one swept gauge angle to be non-monotonic"
    assert not all(outcomes), "expected at least one swept gauge angle to be monotonic"
