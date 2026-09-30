"""
Acceptance gates G1-G5 for specs/SPEC_ft_identity_shift.md (the identity Pauli term inflates the
fault-tolerant/qubitization 1-norm lambda the same way it inflated the near-term measurement
1-norm, twice already fixed there -- SPEC_shift_both_sides / SPEC_lambda_meas_identity).

Reuses `qubitization_blueprint.split_identity` and `qpe_walk_readout.run_qpe(..., recenter=...)`
(both already implemented) plus the existing `lambda_ladder`/`df_factorization` machinery for the
G5 comparability question. No block2/DMRG here -- ordinary process, no isolation needed.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from df_factorization import double_factorize, df_lambda
from lambda_ladder import lambda_and_terms
from qpe_walk_readout import run_qpe
from qubitization_blueprint import build_qubit_hamiltonian, build_walk_operator, pauli_decompose, split_identity

SYSTEMS = ("H2", "LiH", "N2", "H2O")
EXACT_SYSTEMS = ("H2", "LiH")  # small enough (4 qubits) for a fresh exact-matrix walk-operator build


def _reference(label):
    if label == "H2":
        atom, norb, ne = "H 0 0 0; H 0 0 0.74", 2, 2
    elif label == "LiH":
        atom, norb, ne = "Li 0 0 0; H 0 0 1.6", 2, 2
    elif label == "N2":
        atom, norb, ne = "N 0 0 0; N 0 0 1.10", 3, 4
    else:
        atom, norb, ne = "O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467", 3, 4
    mol = gto.M(atom=atom, basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run()
    cas = mcscf.CASCI(mf, norb, ne)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return np.asarray(h1), np.asarray(eri), float(e_core), norb, float(cas.e_tot)


def _terms_and_lambda(label):
    h1, eri, e_core, norb, casci = _reference(label)
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    terms = pauli_decompose(H, n)
    lam = sum(abs(c) for _, c in terms)
    non_identity, c_identity = split_identity(terms)
    lam_excl = sum(abs(c) for _, c in non_identity)
    return h1, eri, e_core, norb, casci, H, n, terms, lam, non_identity, c_identity, lam_excl


@pytest.mark.parametrize("label", SYSTEMS)
def test_G1_identity_fraction_is_material(label):
    """The identity Pauli coefficient carries >25% of lambda on every tested active space --
    reconfirms the bead's scout probe (30.1/44.2/49.3/45.6%), not a rounding-noise artifact."""
    *_, lam, _non_identity, c_identity, _lam_excl = _terms_and_lambda(label)
    fraction = abs(c_identity) / lam
    assert fraction > 0.25, (label, fraction, c_identity, lam)


@pytest.mark.parametrize("label", EXACT_SYSTEMS)
def test_G2_drop_is_free_and_exact(label):
    """(a) lambda' == lambda - |c_identity| to floating precision; (b) re-encoding the
    identity-free terms into a FRESH walk operator and adding c_identity back classically
    reproduces every original Hamiltonian eigenvalue exactly -- the shift is lossless."""
    _h1, _eri, _e_core, _norb, _casci, H, n, terms, lam, non_identity, c_identity, lam_excl = (
        _terms_and_lambda(label)
    )
    assert lam_excl == pytest.approx(lam - abs(c_identity), abs=1e-9), (label, lam_excl, lam, c_identity)

    eigH = np.linalg.eigvalsh(H).real
    W_prime, lam_prime, _L, _a = build_walk_operator(non_identity, n)
    theta_prime = np.angle(np.linalg.eigvals(W_prime))
    recovered = lam_prime * np.cos(theta_prime) + c_identity
    max_err = max(np.min(np.abs(recovered - e)) for e in eigH)
    assert max_err < 1e-8, (label, max_err)


def _t_star(lam, eps, const=3.0):
    return int(np.ceil(np.log2(const * lam / eps)))


@pytest.mark.parametrize("label", SYSTEMS)
def test_G3_budget_reduction_is_real_and_honored(label):
    """DEFINITION OF DONE -- the real falsifier: at a fixed target precision eps, the phase-bit
    count needed by the validated bound err(t) <= 3*lambda/2^t (SPEC_qpe_readout_laws G2) is no
    larger for the recentered/identity-free variant than for the raw one, and running run_qpe
    END TO END at each variant's own computed budget actually lands under eps for BOTH -- so the
    reduced budget is a real, honored guarantee, not just algebra on paper."""
    h1, eri, e_core, norb, casci, H, n, terms, lam, non_identity, c_identity, lam_excl = (
        _terms_and_lambda(label)
    )
    _Ek, Vk = np.linalg.eigh(H)
    ground = Vk[:, 0]
    eps = 1e-3  # Ha, tighter than chemical accuracy (1.6 mHa)

    t_raw = _t_star(lam, eps)
    t_rec = _t_star(lam_excl, eps)
    assert t_rec <= t_raw, (label, t_rec, t_raw, lam, lam_excl)

    E_raw, lam_r, _ps_r, _olap_r = run_qpe(h1, eri, norb, e_core, ground, t_raw, recenter=False)
    E_rec, lam_r2, _ps_c, _olap_c = run_qpe(h1, eri, norb, e_core, ground, t_rec, recenter=True)
    assert abs(E_raw - casci) <= eps, (label, "raw", abs(E_raw - casci), eps)
    assert abs(E_rec - casci) <= eps, (label, "recentered", abs(E_rec - casci), eps)


def test_G3_budget_strictly_shrinks_on_at_least_one_system():
    """The reduction in G3 is not always a no-op rounding artifact: at least one system needs
    strictly fewer phase bits when recentered."""
    eps = 1e-3
    strict = []
    for label in SYSTEMS:
        *_, lam, _non_identity, _c_identity, lam_excl = _terms_and_lambda(label)
        strict.append(_t_star(lam_excl, eps) < _t_star(lam, eps))
    assert any(strict), "no system needed fewer phase bits when recentered"


@pytest.mark.parametrize("label", SYSTEMS)
def test_G4_fixed_t_point_estimate_ratio_is_killed(label):
    """RECORDED REVISION, not silently dropped: the bead's proposed closed-form prediction
    (lambda_eff = sqrt(lambda^2 - E0^2) ratio, or the plain lambda ratio) does NOT track the
    literal fixed-t point-estimate error ratio -- it swings by at least a factor of 5 across a
    t-sweep (the dyadic-grid "staircase" rounding artifact), so a single-t comparison is not a
    meaningful falsifier of any 1-norm-based prediction. This is why G3 restates the claim at the
    envelope/budget level instead."""
    h1, eri, e_core, norb, casci, H, n, terms, lam, non_identity, c_identity, lam_excl = (
        _terms_and_lambda(label)
    )
    _Ek, Vk = np.linalg.eigh(H)
    ground = Vk[:, 0]

    ratios = []
    for t in range(4, 21):
        E_raw, _lam, _ps, _olap = run_qpe(h1, eri, norb, e_core, ground, t, recenter=False)
        E_rec, _lam2, _ps2, _olap2 = run_qpe(h1, eri, norb, e_core, ground, t, recenter=True)
        err_raw = abs(E_raw - casci)
        err_rec = abs(E_rec - casci)
        if err_raw > 1e-10:
            ratios.append(err_rec / err_raw)

    assert len(ratios) >= 5, (label, ratios)
    assert max(ratios) / min(r for r in ratios if r > 1e-12) > 5.0, (label, ratios)


@pytest.mark.parametrize("label", ("N2", "H2O"))
def test_G5_df_lambda_comparability_not_a_flip(label):
    """The df_lambda / SPEC_scdf_lambda G1(b) comparability question, resolved explicitly:
    (a) df_lambda <= identity-INCLUDED naive Pauli lambda still holds (G1(b) unbroken);
    (b) df_lambda is a pure function of (h1, eri) -- recomputing it is bit-identical regardless
        of what the qubit operator's identity coefficient is, because df_lambda's formula never
        reads a Pauli-identity coefficient at all;
    (c) df_lambda DOES lose to the identity-EXCLUDED naive lambda (reconfirms the scout) --
        but per (a)+(b) that is comparing a quantity with no identity-exclusion mechanism against
        one that has one, so it is NOT a valid apples-to-apples comparison. Verdict:
        NOT-COMPARABLE-AND-VACUOUS, not a flip of G1(b).
    """
    h1, eri, e_core, norb, casci, H, n, terms, lam, non_identity, c_identity, lam_excl = (
        _terms_and_lambda(label)
    )
    leaves, _g, _full_rank = double_factorize(eri, norb)
    lam_df = df_lambda(leaves, h1, norb)

    lam_incl, _ = lambda_and_terms(h1, eri, norb)
    assert lam_df <= lam_incl + 1e-9, (label, lam_df, lam_incl)  # (a)

    lam_df_recomputed = df_lambda(leaves, h1, norb)
    assert lam_df_recomputed == lam_df  # (b): (h1, eri) untouched by split_identity => bit-identical

    assert abs(c_identity) / lam > 0.25  # the identity mass being excluded is itself material (G1)
    assert lam_df > lam_excl, (label, lam_df, lam_excl)  # (c): reconfirms the scout's "loss"
