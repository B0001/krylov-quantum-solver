"""
Acceptance gates G1-G4 for specs/SPEC_precision_cost.md.

Claim: certifying the ground energy to precision eps costs the near-term arc N ~ (z lambda_meas/eps)^2
shot-measurements (standard limit, exponent -2, from certified_noise) while FT-QPE costs
Q ~ pi lambda_DF/(2 eps) queries (Heisenberg limit, exponent -1), so the resource ratio ~ 1/eps and
FT wins the EXPONENT. But the FT CONSTANT is not free: raw lambda_DF can exceed lambda_meas (N2), and
only the symmetry shift (scdf_lambda) makes lambda_DF < lambda_meas for every molecule, by a margin
that grows with size.

Real molecular 1-norms (small CAS). PySCF/qiskit, no block2; `make gates` runs it in its own
process.
"""
import numpy as np
from pyscf import ao2mo, gto, mcscf, scf

from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals
from precision_cost import (
    crossover_epsilon,
    ft_queries,
    measurement_lambda,
    near_term_shots,
    qubitization_lambda,
    resource_ratio,
)

_EPS = np.array([1e-2, 1e-3, 1e-4, 1e-5])


def _cas(atom, norb, nelec):
    mf = scf.RHF(gto.M(atom=atom, basis="sto-3g", verbose=0)).run()
    c = mcscf.CASCI(mf, norb, nelec)
    c.kernel()
    h1, ec = c.get_h1eff()
    eri = ao2mo.restore(1, c.get_h2eff(), norb)
    na = (nelec + nelec % 2) // 2
    return np.asarray(h1), np.asarray(eri), float(ec), (na, nelec - na), norb


_MOLS = {"H2": ("H 0 0 0; H 0 0 0.74", 2, 2),
         "H2O": ("O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467", 3, 4),
         "N2": ("N 0 0 0; N 0 0 1.10", 6, 6)}


def _lambdas(atom, norb, nel):
    h1, eri, ec, ne, no = _cas(atom, norb, nel)
    mh = build_hamiltonian_from_integrals(h1, eri, ne, ec)
    return (measurement_lambda(mh),
            qubitization_lambda(h1, eri, no),
            qubitization_lambda(h1, eri, no, nelec=ne, shift=True))


def _slope(cost_fn):
    return float(np.polyfit(np.log(_EPS), np.log([cost_fn(e) for e in _EPS]), 1)[0])


def test_G1_standard_vs_heisenberg_exponents():
    """Near-term shots scale as 1/eps^2 (slope -2), FT queries as 1/eps (slope -1), and the
    resource ratio as 1/eps (slope -1) -- the standard-vs-Heisenberg exponent gap."""
    lm, _, ls = _lambdas(*_MOLS["N2"])
    assert abs(_slope(lambda e: near_term_shots(lm, e)) - (-2.0)) < 1e-6
    assert abs(_slope(lambda e: ft_queries(ls, e)) - (-1.0)) < 1e-6
    assert abs(_slope(lambda e: resource_ratio(lm, ls, e)) - (-1.0)) < 1e-6


def test_G2_ft_win_is_the_exponent_not_the_raw_constant():
    """DEFINITION OF DONE (the finding): the RAW qubitization lambda_DF does NOT uniformly beat the
    measurement 1-norm -- for N2 CAS(6,6) lambda_DF_raw > lambda_meas -- so double factorization
    alone does not shrink the FT constant below measurement."""
    lm_n2, lr_n2, _ = _lambdas(*_MOLS["N2"])
    assert lr_n2 > lm_n2, (lr_n2, lm_n2)                       # raw DF exceeds measurement for N2


def test_G3_symmetry_shift_earns_the_constant_growing_with_size():
    """The symmetry shift (scdf_lambda) drops lambda_DF below lambda_meas for EVERY molecule, by a
    ratio that grows with system size (H2 < H2O < N2) -- the shift is load-bearing for the FT
    constant advantage, not a nicety."""
    ratios = {}
    for name, (atom, norb, nel) in _MOLS.items():
        lm, _, ls = _lambdas(atom, norb, nel)
        assert ls < lm, (name, ls, lm)                        # shifted DF beats measurement
        ratios[name] = lm / ls
    assert ratios["H2"] < ratios["H2O"] < ratios["N2"], ratios  # advantage grows with size
    assert ratios["N2"] > 5.0, ratios["N2"]                    # ~5.7x for N2


def test_G4_crossover_exists_and_ratio_diverges():
    """FT wins asymptotically: the resource-count ratio exceeds 1 (near-term costlier) at chemical
    accuracy and diverges as eps -> 0; and a common-cost crossover eps* exists (finite, positive)
    below which FT is cheaper -- its LOCATION set by the per-query cost, its existence robust."""
    lm, _, ls = _lambdas(*_MOLS["N2"])
    assert resource_ratio(lm, ls, 1.6e-3) > 1e3                # near-term needs >1000x the count
    assert resource_ratio(lm, ls, 1e-5) > resource_ratio(lm, ls, 1e-3)   # diverges as eps->0
    for cpq in (1.0, 1e3, 1e6):                               # any per-query cost -> a finite eps*
        eps_star = crossover_epsilon(lm, ls, cost_per_query=cpq)
        assert np.isfinite(eps_star) and eps_star > 0.0, cpq
        # below eps* FT is cheaper in the common model, above it near-term is
        below = 0.5 * eps_star
        assert near_term_shots(lm, below) > cpq * ft_queries(ls, below), (cpq, eps_star)
