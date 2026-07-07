"""
Acceptance gates G1-G4 for specs/SPEC_cost_advisor.md.

Claim: `cost_advisor.advise` turns the precision-cost bridge into a per-molecule verdict -- near-term
certified Krylov vs FT-QPE, whichever is cheaper at target accuracy eps under a common cost model
whose only free parameter is rho = cost_per_query/cost_per_shot. The verdict flips exactly at the
crossover eps*(rho) ~ 1/rho, the flip-rho at fixed eps is N/Q, and a rho-sweep gives a cross-validated
robust recommendation.

Real molecular 1-norms (small CAS). PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np
from pyscf import ao2mo, gto, mcscf, scf

from cost_advisor import advise, advise_from_integrals, robust_over_rho
from precision_cost import ft_queries, near_term_shots


def _cas(atom, norb, nelec):
    mf = scf.RHF(gto.M(atom=atom, basis="sto-3g", verbose=0)).run()
    c = mcscf.CASCI(mf, norb, nelec)
    c.kernel()
    na = (nelec + nelec % 2) // 2
    return (np.asarray(c.get_h1eff()[0]), ao2mo.restore(1, c.get_h2eff(), norb),
            float(c.get_h1eff()[1]), (na, nelec - na), norb)


# representative real 1-norms (N2 CAS(6,6), SCDF-shifted lambda_DF ~ 4.0 < lambda_meas ~ 22.8)
_LAM_MEAS, _LAM_DF = 22.844, 4.00


def test_G1_verdict_flips_at_eps_star():
    """At fixed rho the verdict is FT just below eps* and near-term just above, and `cheaper` always
    matches the direct cost comparison."""
    rho = 1e4
    v = advise(_LAM_MEAS, _LAM_DF, 1e-3, rho=rho)
    assert advise(_LAM_MEAS, _LAM_DF, 0.5 * v.eps_star, rho=rho).cheaper == "FT"
    assert advise(_LAM_MEAS, _LAM_DF, 2.0 * v.eps_star, rho=rho).cheaper == "near-term"
    for eps in (1e-2, 1e-3, 1e-4, 1e-5):
        w = advise(_LAM_MEAS, _LAM_DF, eps, rho=rho)
        direct = "FT" if rho * ft_queries(_LAM_DF, eps) < near_term_shots(_LAM_MEAS, eps) else "near-term"
        assert w.cheaper == direct, (eps, w.cheaper, direct)


def test_G2_eps_star_and_flip_rho():
    """DEFINITION OF DONE: eps*(rho) ~ 1/rho (a decade of rho shrinks eps* tenfold), and the
    verdict-flip rho at fixed eps equals N(eps)/Q(eps)."""
    eps = 1e-3
    e1 = advise(_LAM_MEAS, _LAM_DF, eps, rho=1e3).eps_star
    e2 = advise(_LAM_MEAS, _LAM_DF, eps, rho=1e4).eps_star
    assert abs(e1 / e2 - 10.0) < 1e-6, (e1, e2)                       # eps* ~ 1/rho
    v = advise(_LAM_MEAS, _LAM_DF, eps, rho=1.0)
    assert abs(v.rho_flip - near_term_shots(_LAM_MEAS, eps) / ft_queries(_LAM_DF, eps)) < 1e-6
    # exactly at rho = rho_flip the two costs are equal (boundary)
    at = advise(_LAM_MEAS, _LAM_DF, eps, rho=v.rho_flip)
    assert abs(at.rho * at.queries - at.shots) < 1e-3 * at.shots


def test_G3_cross_validated_robustness():
    """At chemical accuracy, FT is robustly cheaper across rho in [1, 1e6] unless rho exceeds the
    (large) flip threshold; `robust_over_rho` reports the correct invariance."""
    eps = 1.6e-3
    rho_flip = advise(_LAM_MEAS, _LAM_DF, eps).rho_flip
    assert rho_flip > 1e5                                              # FT wins to ~2e5x per-query cost
    verdict, robust = robust_over_rho(_LAM_MEAS, _LAM_DF, eps, 1.0, 1e5)
    assert robust and verdict == "FT", (verdict, robust)              # robust below the flip
    # a range that straddles the flip is reported as mixed
    _, robust2 = robust_over_rho(_LAM_MEAS, _LAM_DF, eps, 0.1 * rho_flip, 10.0 * rho_flip)
    assert robust2 is False


def test_G4_integrals_path_and_rho_exposed():
    """`advise_from_integrals` (SCDF shift) agrees with `advise` on the same lambdas, and the Verdict
    exposes rho so the unknown is never hidden."""
    h1, eri, ec, ne, no = _cas("N 0 0 0; N 0 0 1.10", 6, 6)
    v = advise_from_integrals(h1, eri, ec, ne, no, 1.6e-3, rho=1e3, shift=True)
    assert abs(v.lam_df - 4.0) < 0.2                                   # SCDF-shifted N2 lambda_DF
    assert v.cheaper == advise(v.lam_meas, v.lam_df, v.eps, rho=1e3).cheaper
    assert v.rho == 1e3                                                # rho is explicit on the verdict
