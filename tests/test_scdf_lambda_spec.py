"""
Acceptance gates G1-G4 for specs/SPEC_scdf_lambda.md (factorization-native lambda + symmetry shift).

Test-first: ``df_factorization.df_lambda`` and ``df_factorization.symmetry_shift`` do not exist yet,
so this file is RED until the spec is implemented.

The native double-factorization 1-norm is
    lambda_DF = 1/4 * sum_t |g_t| (sum_k |w^(t)_k|)^2  +  ||h1||_nuc
computed directly from the `double_factorize` leaves -- NO brute-force Pauli enumeration, so it
scales past ~4 orbitals (the gap `lambda_ladder.py` names in its docstring). The number-operator
symmetry shift H -> H + (a1 N_e + a2)(N_e - n_e) leaves the spectrum invariant (a posteriori
constant correction) while lowering lambda.

Uses only pyscf (no block2), small CASes so the Pauli oracle is computable.
"""
import numpy as np
import pytest
from pyscf import gto, scf, mcscf, ao2mo, fci

from df_factorization import double_factorize, df_lambda, symmetry_shift  # RED until implemented
from lambda_ladder import lambda_and_terms


def _active_space(atom, norb, nelec, basis="sto-3g"):
    """Return (h1, eri, e_core, (na, nb), norb, casci_energy) for a small CAS."""
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    cas = mcscf.CASCI(mf, norb, nelec)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    na = nb = nelec // 2
    return np.asarray(h1), np.asarray(eri), float(e_core), (na, nb), norb, float(cas.e_tot)


def _n2_small():
    return _active_space("N 0 0 0; N 0 0 1.10", norb=3, nelec=4)


def _h2o_small():
    # norb=3 (6 qubits) keeps the brute-force Pauli oracle in G1 cheap.
    return _active_space("O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467", norb=3, nelec=4)


def _lambda_recompute(leaves, h1):
    """Independent in-test recomputation of the bare DF 1-norm (different code path)."""
    two_body = 0.25 * sum(abs(gt) * (np.abs(w).sum()) ** 2 for gt, w, _ in leaves)
    one_body = float(np.linalg.svd(0.5 * (h1 + h1.T), compute_uv=False).sum())
    return two_body + one_body


def test_G1_native_lambda_formula():
    """(a) df_lambda matches an independent recompute; (b) df_lambda <= naive Pauli LCU lambda."""
    for h1, eri, e_core, nelec, norb, _ in (_n2_small(), _h2o_small()):
        leaves, _, _ = double_factorize(eri, norb)            # full rank
        lam = df_lambda(leaves, h1, norb)
        assert lam == pytest.approx(_lambda_recompute(leaves, h1), abs=1e-9)
        lam_pauli, _ = lambda_and_terms(h1, eri, norb)
        # G1(b): DF is no looser than the trivial Pauli LCU (revisable per spec if a tiny CAS fails)
        assert lam <= lam_pauli + 1e-9, (lam, lam_pauli)


def test_G2_shift_preserves_spectrum():
    """FCI(shifted) + correction == FCI(unshifted) to < 1e-8 Ha (the shift must not move E)."""
    for h1, eri, e_core, nelec, norb, casci in (_n2_small(), _h2o_small()):
        e_unshifted, _ = fci.direct_spin1.kernel(h1, eri, norb, nelec, ecore=e_core)
        h1_s, eri_s, e_shift, _ = symmetry_shift(h1, eri, norb, nelec)
        e_shifted, _ = fci.direct_spin1.kernel(h1_s, eri_s, norb, nelec, ecore=e_core)
        assert abs((e_shifted + e_shift) - e_unshifted) < 1e-8, (e_shifted + e_shift, e_unshifted)
        assert abs(e_unshifted - casci) < 1e-8                # sanity: matches CASCI


def test_G3_shift_lowers_lambda():
    """DEFINITION OF DONE: the symmetry shift drops lambda_DF by >= 20% on N2 CAS(6,6)."""
    h1, eri, e_core, nelec, norb, _ = _active_space("N 0 0 0; N 0 0 1.10", norb=6, nelec=6)
    leaves, _, _ = double_factorize(eri, norb)
    lam_before = df_lambda(leaves, h1, norb)

    h1_s, eri_s, _, _ = symmetry_shift(h1, eri, norb, nelec)
    leaves_s, _, _ = double_factorize(eri_s, norb)
    lam_after = df_lambda(leaves_s, h1_s, norb)

    reduction = 1.0 - lam_after / lam_before
    assert reduction >= 0.20, f"lambda reduction {reduction:.1%} (before={lam_before:.4f}, " \
                              f"after={lam_after:.4f}) -- revise threshold + record if smaller"


def test_G4_lambda_monotone_in_rank():
    """Truncating DF rank upward, df_lambda increases toward the full-rank value."""
    h1, eri, e_core, nelec, norb, _ = _h2o_small()
    _, _, full_rank = double_factorize(eri, norb)
    lams = []
    for R in range(1, full_rank + 1):
        leaves, _, _ = double_factorize(eri, norb, rank=R)
        lams.append(df_lambda(leaves, h1, norb))
    assert all(lams[i + 1] >= lams[i] - 1e-9 for i in range(len(lams) - 1)), lams
