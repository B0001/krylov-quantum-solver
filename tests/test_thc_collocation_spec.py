"""
Acceptance gates G1-G4 for specs/SPEC_thc_collocation.md (chem-5oj): scores `lambda_ladder.fit_thc`'s
nonlinear THC fit -- and a new deterministic non-random baseline, `pair_indicator_collocation` --
with the native `thc_lambda`, at matched rank M=norb(norb+1)/2, against random collocation
(`tensor_hypercontraction`) and `df_lambda` as the reference. Revisits SPEC_thc_lambda's G4
(locked: unoptimized/random collocation does not beat df_lambda) with a check that spec's own text
invited: does a non-random collocation change that?

System: LiH/STO-3G, norb=6 (the CI-gate cap chem-5oj's cost note requires; norb=7/H2O is reported
in the handoff, not gated here). PySCF/NumPy/SciPy only, no block2.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from df_factorization import double_factorize, df_lambda
from thc_factorization import (
    pair_indicator_collocation,
    reconstruct_thc,
    tensor_hypercontraction,
    thc_lambda,
    thc_rank,
)
from lambda_ladder import fit_thc


@pytest.fixture(scope="module")
def lih_system():
    mol = gto.M(atom="Li 0 0 0; H 0 0 1.6", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    norb = mol.nao_nr()
    assert norb == 6, "chem-5oj's CI-gate cap assumes LiH/STO-3G gives norb=6"
    na = nb = mol.nelectron // 2
    cas = mcscf.CASCI(mf, norb, (na, nb))
    cas.verbose = 0
    cas.kernel()
    h1, ecore = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return dict(h1=h1, eri=eri, norb=norb, nelec=(na, nb), ecore=ecore, e_fci=cas.e_tot)


@pytest.fixture(scope="module")
def df_baseline(lih_system):
    c = lih_system
    leaves, _, full_rank = double_factorize(c["eri"], c["norb"])
    lam_df = df_lambda(leaves, c["h1"], c["norb"])
    return lam_df, full_rank


@pytest.fixture(scope="module")
def random_collocation(lih_system):
    c = lih_system
    M = thc_rank(c["norb"])
    chi, zeta = tensor_hypercontraction(c["eri"], c["norb"], n_thc=M, seed=0)
    return chi, zeta


@pytest.fixture(scope="module")
def structured_collocation(lih_system):
    chi, zeta = pair_indicator_collocation(lih_system["eri"], lih_system["norb"])
    return chi, zeta


@pytest.fixture(scope="module")
def nonlinear_collocation(lih_system):
    """fit_thc at its own fast CI defaults (restarts=4, max_nfev=4000), seed pinned to 0 --
    matches the module's existing default so this reproduces what `lambda_ladder()`'s print loop
    already runs, just scored differently."""
    c = lih_system
    M = thc_rank(c["norb"])
    X, Z = fit_thc(c["eri"], c["norb"], M, restarts=4, seed=0, return_factors=True)
    return X, Z


def test_G1_random_and_structured_collocation_reconstruct_exactly_at_matched_rank(
    lih_system, random_collocation, structured_collocation
):
    """Both linear (non-nonlinear-fit) collocations hit machine precision at M=thc_rank(norb) --
    matched rank and matched (exact) reconstruction is achievable in closed form for these two."""
    c = lih_system
    M = thc_rank(c["norb"])
    chi_r, zeta_r = random_collocation
    chi_s, zeta_s = structured_collocation
    assert chi_r.shape == (c["norb"], M)
    assert chi_s.shape == (c["norb"], M)
    err_r = np.linalg.norm(reconstruct_thc(chi_r, zeta_r) - c["eri"])
    err_s = np.linalg.norm(reconstruct_thc(chi_s, zeta_s) - c["eri"])
    assert err_r < 1e-9, err_r
    assert err_s < 1e-9, err_s


def test_G2_fit_thc_fails_the_matched_reconstruction_precondition(
    lih_system, nonlinear_collocation
):
    """RECORDED FINDING, not a bug: at the SAME matched rank M=21, fit_thc's own fast defaults do
    NOT reach the <1e-6 reconstruction fidelity the check's precondition asks for -- the honest
    boundary this spec exists to pin. Corroborated by a larger, non-CI-gated sweep (12 seeds x
    8000 LM evaluations each) that also never got below a ~3.6e-3 residual -- see the PR/handoff;
    not re-run here to keep the gate fast."""
    c = lih_system
    X, Z = nonlinear_collocation
    err = np.linalg.norm(reconstruct_thc(X, Z) - c["eri"])
    assert err > 1e-3, (
        "fit_thc reached the <1e-6 reconstruction precondition -- SPEC_thc_collocation.md's "
        f"honest-boundary claim (G2) is now wrong and must be revised. err={err}"
    )


def test_G3_lambda_formula_sanity_anchors(lih_system, df_baseline, structured_collocation):
    """thc_lambda on the deterministic structured collocation is positive, finite, and within 2x
    of df_lambda -- a sanity bound (NOT claimed exactly equal; that exact-equality claim belongs to
    thc_from_df / SPEC_thc_lambda G2, a different, larger-rank structured THC)."""
    lam_df, _ = df_baseline
    chi_s, zeta_s = structured_collocation
    lam_s = thc_lambda(chi_s, zeta_s, lih_system["h1"])
    assert np.isfinite(lam_s) and lam_s > 0, lam_s
    assert 0.5 * lam_df < lam_s < 2.0 * lam_df, (lam_s, lam_df)


def test_G4_both_kill_directions_evaluated_explicitly(
    lih_system, df_baseline, random_collocation, nonlinear_collocation
):
    """THE COMPARISON chem-5oj asks for. Both directions are computed and recorded; neither is
    forced to a predetermined verdict:

    - Kill A ("penalty runs deeper than unoptimized points"): nonlinear lambda NOT >=5x below
      random lambda. MEASURED: it IS >=5x below (by ~60-68x across seeds 0-4 checked informally;
      seed=0's exact ratio asserted below) -- Kill A does NOT fire.
    - Kill B ("SPEC_thc_lambda G4 must be revised"): nonlinear lambda beats df_lambda outright.
      MEASURED at seed=0: it does (ratio ~0.94). But read this together with G2: fit_thc's
      reconstruction error at this seed is ~0.017 on ||eri||~2.67 (~0.6% relative) -- NOT the
      near-machine-precision reconstruction G4's original claim was about. A poorly-reconstructing
      operator can have an arbitrarily small 1-norm (the zero operator has lambda=0 and infinite
      error), so "beats df_lambda while failing G2" is reported as a genuine, reproducible number,
      not as proof the 62x penalty is solved -- see SPEC_thc_collocation.md Sec 8 R1.

    SPEC_thc_lambda's own G4 (about the EXACT random linear-LS THC, tensor_hypercontraction) is
    UNCHANGED and re-confirmed here: random collocation is still >>1x above df_lambda (~65x on
    this system, consistent with the ~60-62x previously recorded on H2O/N2).
    """
    lam_df, _ = df_baseline
    chi_r, zeta_r = random_collocation
    lam_r = thc_lambda(chi_r, zeta_r, lih_system["h1"])
    X, Z = nonlinear_collocation
    lam_nl = thc_lambda(X, Z, lih_system["h1"])

    # SPEC_thc_lambda G4 (exact random collocation vs DF) is unchanged: still far above df_lambda.
    assert lam_r > 10.0 * lam_df, (lam_r, lam_df)

    ratio_to_random = lam_nl / lam_r
    beats_df = lam_nl < lam_df

    # Kill A check: does NOT fire -- nonlinear collocation is comfortably >=5x below random.
    assert ratio_to_random < 0.2, ("Kill A fires: nonlinear collocation is not >=5x below random",
                                    ratio_to_random, lam_nl, lam_r)
    # Kill B, recorded as measured at the pinned seed (not asserted as a universal law -- a 5-seed
    # spot check during development found 4/5 seeds beat df_lambda and 1/5 did not).
    assert beats_df, ("expected seed=0's fit to beat df_lambda at this pinned seed", lam_nl, lam_df)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
