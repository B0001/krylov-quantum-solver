"""
Acceptance gates G1-G4 for specs/SPEC_thc_lambda.md (tensor hypercontraction + qubitization lambda).

Test-first: ``thc_factorization`` does not exist yet, so this file is RED until the spec is
implemented. THC writes (pq|rs) ~ sum_{mu,nu} chi_p^mu chi_q^mu zeta_{mu,nu} chi_r^nu chi_s^nu. We
validate (G1) exact reconstruction + FCI at THC rank norb(norb+1)/2, (G2) the THC 1-norm against the
already-validated df_lambda on the DF-derived structured THC, (G3) the rank compression vs DF-THC,
and (G4) the lambda formula's scale invariance plus the recorded boundary that unoptimized-collocation
THC does NOT beat DF's lambda.

PySCF/NumPy only (no block2), so it runs in the non-DMRG process group -- but `make gates` runs every
test_*_spec.py in its own process anyway.
"""
import numpy as np
import pytest

from df_factorization import double_factorize, df_lambda
from thc_factorization import (
    reconstruct_thc,
    tensor_hypercontraction,
    thc_from_df,
    thc_lambda,
    thc_rank,
)


def _cas(atom):
    """Active-space integrals + FCI for a small full-space molecule (STO-3G)."""
    from pyscf import ao2mo, gto, mcscf, scf
    mol = gto.M(atom=atom, basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    norb = mol.nao_nr()
    na = nb = mol.nelectron // 2
    cas = mcscf.CASCI(mf, norb, (na, nb))
    cas.verbose = 0
    cas.kernel()
    h1, ecore = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return dict(h1=h1, eri=eri, norb=norb, nelec=(na, nb), ecore=ecore, e_fci=cas.e_tot)


H2O = "O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467"
N2 = "N 0 0 0; N 0 0 1.1"


def test_G1_exact_reconstruction_and_fci():
    """DEFINITION OF DONE: linear-LS THC at M=norb(norb+1)/2 reconstructs ERIs exactly, keeps FCI."""
    from pyscf import fci
    for atom in (H2O, N2):
        c = _cas(atom)
        chi, zeta = tensor_hypercontraction(c["eri"], c["norb"])
        assert chi.shape == (c["norb"], thc_rank(c["norb"]))
        eri_thc = reconstruct_thc(chi, zeta)
        assert np.linalg.norm(eri_thc - c["eri"]) < 1e-9, np.linalg.norm(eri_thc - c["eri"])
        e_thc, _ = fci.direct_spin1.kernel(c["h1"], eri_thc, c["norb"], c["nelec"])
        assert abs((e_thc + c["ecore"]) - c["e_fci"]) < 1e-6, (e_thc + c["ecore"], c["e_fci"])


def test_G2_lambda_matches_df_on_structured_thc():
    """The THC 1-norm equals df_lambda on the DF-derived THC -- anchors the formula to vetted code."""
    for atom in (H2O, N2):
        c = _cas(atom)
        chi_df, zeta_df = thc_from_df(c["eri"], c["norb"])
        lam_thc = thc_lambda(chi_df, zeta_df, c["h1"])
        leaves, _, _ = double_factorize(c["eri"], c["norb"])
        lam_df = df_lambda(leaves, c["h1"], c["norb"])
        assert abs(lam_thc - lam_df) < 1e-9, (lam_thc, lam_df)


def test_G3_rank_compression():
    """THC rank norb(norb+1)/2 is strictly below the DF-derived THC rank norb*(DF full rank)."""
    for atom in (H2O, N2):
        c = _cas(atom)
        _, _, full = double_factorize(c["eri"], c["norb"])
        assert thc_rank(c["norb"]) < c["norb"] * full, (thc_rank(c["norb"]), c["norb"] * full)


def test_G4_lambda_scale_invariance_and_boundary():
    """lambda is invariant under chi rescaling; and the recorded finding: naive THC does not beat DF."""
    c = _cas(H2O)
    chi, zeta = tensor_hypercontraction(c["eri"], c["norb"])
    lam = thc_lambda(chi, zeta, c["h1"])
    # scale invariance: chi -> 2 chi, zeta -> zeta/16 leaves the operator (and lambda) unchanged
    lam_scaled = thc_lambda(2.0 * chi, zeta / 16.0, c["h1"])
    assert abs(lam_scaled - lam) / lam < 1e-6, (lam_scaled, lam)
    # recorded boundary: unoptimized collocation gives a much LARGER lambda than DF (~60x here)
    leaves, _, _ = double_factorize(c["eri"], c["norb"])
    lam_df = df_lambda(leaves, c["h1"], c["norb"])
    assert lam > lam_df, (lam, lam_df)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
