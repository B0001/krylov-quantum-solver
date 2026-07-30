"""Gates for specs/SPEC_lambda_h2_bridge.md.

Two claims: (a) `SPEC_lambda_meas_identity`'s identity fix never reached
`certified_noise.hamiltonian_one_norms`, and (b) the near-term/FT bridge charges only lambda_H
although the Temple bracket it represents needs <H^2> as well.

The naive "the undercount grows with system size" reading is FALSE and G5 pins the counterexample.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from certified_noise import hamiltonian_one_norms
from cost_advisor import advise
from hybrid_quantum_solver.molecular_hamiltonian import (
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
)
from precision_cost import measurement_lambda, qubitization_lambda

EPS = 1.6e-3                      # chemical accuracy, the repo's standard advisor target

# The SPEC_certified_noise case list (G3 must not silently narrow it).
_CERTIFIED_NOISE_CASES = {
    "H2": "H 0 0 0; H 0 0 0.74",
    "H4": "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0",
}

# Homogeneous family for the growth claim: H_n / STO-3G, full space.
_HN_FAMILY = [
    ("H2", "H 0 0 0; H 0 0 0.74"),
    ("H4", "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
    ("H6", "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0; H 0 0 4.0; H 0 0 5.0"),
]


def _mh(atom):
    return build_molecular_hamiltonian(atom, basis="sto-3g")


def _legacy_one_norms(mh):
    """The pre-fix expression from certified_noise.py:47, verbatim -- the G2 archaeology oracle."""
    op = mh.qubit_hamiltonian
    return (float(np.abs(op.coeffs).sum()),
            float(np.abs((op @ op).simplify().coeffs).sum()))


def _cas(atom, norb, nelec):
    """CASCI active-space integrals, same helper shape as tests/test_cost_advisor_spec.py."""
    mf = scf.RHF(gto.M(atom=atom, basis="sto-3g", verbose=0)).run()
    c = mcscf.CASCI(mf, norb, nelec)
    c.kernel()
    na = (nelec + nelec % 2) // 2
    return (np.asarray(c.get_h1eff()[0]), ao2mo.restore(1, c.get_h2eff(), norb),
            float(c.get_h1eff()[1]), (na, nelec - na), norb)


def _ratio_from_cas(atom, norb, nelec):
    h1, eri, ec, ne, _ = _cas(atom, norb, nelec)
    mh = build_hamiltonian_from_integrals(h1, eri, ne, ec)
    lam_h, lam_h2 = hamiltonian_one_norms(mh)
    return lam_h, lam_h2, mh.qubit_hamiltonian.num_qubits


# --- G1: the two modules must agree on lambda_H (DEFINITION OF DONE) -----------------------------

@pytest.mark.parametrize("atom", list(_CERTIFIED_NOISE_CASES.values()))
def test_G1_one_norms_match_the_measurement_lambda_convention(atom):
    """THE DEFECT. `precision_cost.measurement_lambda` excludes the identity term (a constant of
    zero variance costs zero shots -- SPEC_lambda_meas_identity). `hamiltonian_one_norms` was never
    updated, so the two modules disagreed by the identity fraction.
    """
    mh = _mh(atom)
    lam_h, _ = hamiltonian_one_norms(mh)
    assert lam_h == pytest.approx(measurement_lambda(mh), rel=0, abs=1e-12)


def test_G1_identity_fraction_is_material():
    """Not a rounding difference: the identity carried ~22-30% of the 1-norm."""
    for atom in _CERTIFIED_NOISE_CASES.values():
        mh = _mh(atom)
        excl, _ = hamiltonian_one_norms(mh)
        incl, _ = hamiltonian_one_norms(mh, include_identity=True)
        assert 0.15 < 1.0 - excl / incl < 0.40, (atom, excl, incl)


# --- G2: archaeology -- the old inflated values stay reachable ------------------------------------

@pytest.mark.parametrize("atom", list(_CERTIFIED_NOISE_CASES.values()))
def test_G2_include_identity_reproduces_the_pre_fix_values(atom):
    mh = _mh(atom)
    assert hamiltonian_one_norms(mh, include_identity=True) == pytest.approx(
        _legacy_one_norms(mh), rel=1e-12)


# --- G3: SPEC_certified_noise's recorded boundary survives the fix --------------------------------

@pytest.mark.parametrize("name,atom", list(_CERTIFIED_NOISE_CASES.items()))
def test_G3_H2_side_is_still_the_noise_expensive_one(name, atom):
    """SPEC_certified_noise G4 recorded lambda_H2 > lambda_H ("the Temple lower bound is the
    noise-expensive side"). The identity fix must not overturn that finding -- it strengthens it
    (H4 ratio 6.3 -> 7.09).
    """
    lam_h, lam_h2 = hamiltonian_one_norms(_mh(atom))
    assert lam_h2 > lam_h, (name, lam_h, lam_h2)


# --- G4: growth WITHIN a homogeneous family -------------------------------------------------------

def test_G4_ratio_grows_across_the_Hn_family():
    """H_n/STO-3G, full space: the undercount ratio strictly increases with n."""
    ratios = []
    for _, atom in _HN_FAMILY:
        lam_h, lam_h2 = hamiltonian_one_norms(_mh(atom))
        ratios.append(lam_h2 / lam_h)
    assert all(a < b for a, b in zip(ratios, ratios[1:])), ratios
    assert ratios[-1] / ratios[0] > 5.0, ratios      # 1.74 -> 16.69 measured


# --- G5: the boundary that KILLS the naive "grows with size" law ----------------------------------

def test_G5_ratio_is_not_a_function_of_qubit_count():
    """THE KILLED FORM. The backlog entry predicted growth "with size". Across heterogeneous active
    spaces that is false: LiH CAS(2,3) at 6 qubits sits BELOW H2 CAS(2,2) at 4 qubits, so the
    undercount cannot be predicted from system size alone and SPEC_precision_cost's "margin grows
    with size" headline is confounded, not cleanly cancelled.
    """
    lam_h_h2, lam_h2_h2, nq_h2 = _ratio_from_cas("H 0 0 0; H 0 0 0.74", 2, 2)
    lam_h_lih, lam_h2_lih, nq_lih = _ratio_from_cas("Li 0 0 0; H 0 0 1.6", 3, 2)
    r_h2, r_lih = lam_h2_h2 / lam_h_h2, lam_h2_lih / lam_h_lih
    assert nq_lih > nq_h2, (nq_h2, nq_lih)           # LiH really is the bigger register
    assert r_lih < r_h2, (r_h2, r_lih)               # ...yet the smaller ratio


# --- G6: the undercount is material, and moves verdicts one way only ------------------------------

def test_G6_pricing_H2_moves_verdicts_and_only_toward_FT():
    """With REAL SCDF-shifted lambda_DF, charging the <H^2> moment moves cost_advisor verdicts, and
    every move runs near-term -> FT.

    The inflation lambda_H + lambda_H2 is a labelled PROXY (SPEC R1), not a derivation: the gate
    asserts direction and existence of movement, never a crossover value.
    """
    cases = [("H2", "H 0 0 0; H 0 0 0.74", 2, 2),
             ("LiH", "Li 0 0 0; H 0 0 1.6", 3, 2),
             ("N2", "N 0 0 0; N 0 0 1.1", 6, 6)]
    moves = 0
    for name, atom, norb, nelec in cases:
        h1, eri, ec, ne, no = _cas(atom, norb, nelec)
        mh = build_hamiltonian_from_integrals(h1, eri, ne, ec)
        lam_h, lam_h2 = hamiltonian_one_norms(mh)
        lam_df = qubitization_lambda(h1, eri, no, nelec=ne, shift=True)
        for rho in (1e0, 1e2, 1e3, 1e4, 1e5, 1e6):
            before = advise(lam_h, lam_df, EPS, rho=rho).cheaper
            after = advise(lam_h + lam_h2, lam_df, EPS, rho=rho).cheaper
            if before != after:
                moves += 1
                # One-way: pricing a cost the near-term side actually pays can only favour FT.
                assert (before, after) == ("near-term", "FT"), (name, rho, before, after)
    assert moves > 0, "no verdict moved -- the undercount would be immaterial to the advisor"
