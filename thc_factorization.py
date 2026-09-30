#!/usr/bin/env python3
"""
Tensor hypercontraction (THC) of the ERI tensor -- a more compact (in rank) factorization than
double factorization, and the asymptotically cheaper qubitization route (Lee et al., "Even more
efficient quantum computations of chemistry through tensor hypercontraction", PRX Quantum 2,
030305, 2021).

THC writes the two-electron integrals as a single contraction over a collocation matrix chi
(norb x M) and a symmetric central matrix zeta (M x M):

    (pq|rs) ~ sum_{mu,nu} chi_p^mu chi_q^mu  zeta_{mu,nu}  chi_r^nu chi_s^nu

Each THC index mu carries a rank-1 symmetric "density" v^mu (v^mu)^T with v^mu_p = chi_p^mu, so the
operator sum_{pq} chi_p^mu chi_q^mu E_pq is a (scaled) rotated number operator -- exactly the
structure qubitization needs, and the same structure double factorization produces. In fact DF is a
*special structured THC*: substituting each DF leaf L^t_pq = sum_k w^t_k u^tk_p u^tk_q into
(pq|rs) = sum_t g_t L^t_pq L^t_rs gives a THC with chi columns {u^tk} and zeta block-diagonal
(zeta_{(t,k),(t,l)} = g_t w^t_k w^t_l). ``thc_from_df`` returns exactly that, and ``thc_lambda`` on
it reproduces ``df_factorization.df_lambda`` to machine precision -- the anchor that validates the
THC 1-norm here against already-validated code.

HONEST SCOPE (see specs/SPEC_thc_lambda.md): the linear-least-squares THC below reconstructs the
ERIs *exactly* at THC rank M = norb(norb+1)/2 (the symmetric-pair dimension) and compresses the rank
relative to DF-THC, but with *unoptimized* collocation its qubitization 1-norm is far larger than
DF's. The literature THC lambda advantage requires *optimized* (ISDF / nonlinear) collocation points
-- a nonlinear fit that is deliberately out of scope. What is delivered is a correct THC
factorization, a validated lambda, the rank compression, and a precise statement of where the lambda
advantage does and does not come from.

Same active-space inputs as the DF / SQD / Krylov modules: eri from cas.get_h2eff().
"""

import numpy as np

from df_factorization import double_factorize


def thc_rank(norb):
    """The symmetric orbital-pair dimension norb(norb+1)/2 -- the THC rank at which a full-rank
    collocation reconstructs the ERIs exactly."""
    return norb * (norb + 1) // 2


def thc_from_df(eri, norb, rank=None):
    """Exact *structured* THC read straight off the double factorization (no fitting).

    Returns ``(chi, zeta)`` with ``chi`` shape ``(norb, M)`` (M = sum over kept DF leaves of norb)
    and ``zeta`` the block-diagonal central matrix. Reconstructs the (truncated-rank) ERIs exactly
    and gives ``thc_lambda(chi, zeta, h1) == df_lambda(leaves, h1)``. This is the validation anchor,
    not a compression -- its rank M = norb * (DF rank) is larger than ``tensor_hypercontraction``'s.
    """
    leaves, _, _ = double_factorize(eri, norb, rank=rank)
    cols, diag = [], []
    for gt, w, U in leaves:
        for k in range(norb):
            cols.append(U[:, k])
            diag.append((gt, w[k]))
    chi = np.array(cols).T                                  # (norb, M)
    M = chi.shape[1]
    zeta = np.zeros((M, M))
    # block-diagonal per leaf: zeta_{(t,k),(t,l)} = g_t w_k w_l
    off = 0
    for gt, w, _ in leaves:
        n = len(w)
        zeta[off:off + n, off:off + n] = gt * np.outer(w, w)
        off += n
    return chi, zeta


def zeta_from_collocation(eri, norb, chi):
    """Solve the central matrix ``zeta`` by linear least squares for an explicit collocation
    ``chi`` (shape ``(norb, n_thc)``). Shared by every linear (non-nonlinear-fit) THC collocation
    scheme below: the reconstruction is exact iff ``chi``'s induced pair-design matrix
    ``P_{mu,(pq)} = chi_p^mu chi_q^mu`` has rank ``n_thc`` -- i.e. iff the pair space is spanned.
    """
    n_thc = chi.shape[1]
    V = eri.reshape(norb * norb, norb * norb)
    P = np.einsum("pm,qm->mpq", chi, chi).reshape(n_thc, norb * norb)
    Pinv = np.linalg.pinv(P.T)                              # (n_thc, norb^2)
    zeta = Pinv @ V @ Pinv.T
    return 0.5 * (zeta + zeta.T)


def tensor_hypercontraction(eri, norb, n_thc=None, seed=0):
    """Linear-least-squares THC: random full-rank collocation, central matrix by least squares.

    With ``n_thc >= norb(norb+1)/2`` the symmetric-pair space is spanned and the reconstruction is
    exact (to round-off) -- no nonlinear iteration. Returns ``(chi, zeta)`` with ``chi`` shape
    ``(norb, n_thc)`` and symmetric ``zeta`` shape ``(n_thc, n_thc)``.

    Note: this is a *valid and exact* THC, but the collocation is not optimized for a small 1-norm;
    see the module docstring and the spec for the lambda caveat.
    """
    n_thc = thc_rank(norb) if n_thc is None else int(n_thc)
    rng = np.random.default_rng(seed)
    chi = rng.standard_normal((norb, n_thc)) * 0.5
    zeta = zeta_from_collocation(eri, norb, chi)
    return chi, zeta


def pair_indicator_collocation(eri, norb):
    """Deterministic, non-random, non-optimized THC collocation at ``n_thc = thc_rank(norb)``:
    one collocation index per symmetric orbital pair, ``chi^{(i,i)} = e_i`` and
    ``chi^{(i,j)} = e_i + e_j`` (i<j). Central matrix by the same linear least squares as
    ``tensor_hypercontraction``.

    This is "structured" only in the weak sense of being a fixed combinatorial formula with no RNG
    seed and no search for small lambda -- it is NOT the ISDF/optimized collocation the spec puts
    out of scope. It is provably full column rank (each basis direction e_i e_i^T is directly a
    column; each e_i e_j^T is recovered from the (i,i), (j,j), (i,j) columns by
    ``((e_i+e_j)(e_i+e_j)^T - e_ie_i^T - e_je_j^T)/2``), so reconstruction is exact by construction,
    not by numerical luck -- unlike a smooth 1-parameter deterministic curve (Fourier/Chebyshev/
    Vandermonde collocation), which was tried first and tops out at rank 2*norb-1 (see
    specs/SPEC_thc_collocation.md).
    """
    cols = []
    for i in range(norb):
        e = np.zeros(norb)
        e[i] = 1.0
        cols.append(e)
    for i in range(norb):
        for j in range(i + 1, norb):
            e = np.zeros(norb)
            e[i] = 1.0
            e[j] = 1.0
            cols.append(e)
    chi = np.array(cols).T                                  # (norb, norb(norb+1)/2)
    zeta = zeta_from_collocation(eri, norb, chi)
    return chi, zeta


def reconstruct_thc(chi, zeta):
    """Rebuild the ERI tensor (pq|rs) from its THC factors."""
    return np.einsum("pm,qm,mn,rn,sn->pqrs", chi, chi, zeta, chi, chi, optimize=True)


def thc_lambda(chi, zeta, h1):
    """Qubitization 1-norm of the THC two-body operator plus the one-body nuclear norm.

        lambda_THC = 1/4 * sum_{mu,nu} |zeta_{mu,nu}| ||v^mu||^2 ||v^nu||^2  +  ||h1||_nuc

    Each THC density v^mu (v^mu)^T is a rank-1 symmetric one-body operator with single nonzero
    eigenvalue ||v^mu||^2 = sum_p (chi_p^mu)^2, so sum_k |eigenvalue_k| = ||v^mu||^2. This is the
    same qubitization 1-norm convention as ``df_factorization.df_lambda`` (1/4 |g_t|(sum_k|w_k|)^2
    per leaf), specialized to the rank-1 THC densities; on ``thc_from_df`` output the two agree
    exactly. The value is invariant under rescaling chi_mu -> c chi_mu (zeta absorbs c^4), as a
    1-norm of a fixed operator must be.
    """
    vn = np.sum(np.asarray(chi) ** 2, axis=0)              # ||v^mu||^2 per THC index
    two_body = 0.25 * np.sum(np.abs(zeta) * np.outer(vn, vn))
    h1 = np.asarray(h1, dtype=float)
    one_body = float(np.linalg.svd(0.5 * (h1 + h1.T), compute_uv=False).sum())
    return float(two_body + one_body)


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo
    from df_factorization import df_lambda
    from hybrid_quantum_solver.dmrg_reference import fci_energy

    for atom, name in [("O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467", "H2O"),
                       ("N 0 0 0; N 0 0 1.1", "N2")]:
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

        leaves, _, full = double_factorize(eri, norb)
        chi, zeta = tensor_hypercontraction(eri, norb)
        eri_thc = reconstruct_thc(chi, zeta)
        e_thc = fci_energy(h1, eri_thc, (na, nb), e_core=ecore)
        chi_df, zeta_df = thc_from_df(eri, norb)
        print("=" * 78)
        print(f"{name}: norb={norb}  FCI={cas.e_tot:.6f}")
        print(f"  THC rank M = norb(norb+1)/2 = {thc_rank(norb)}   (DF-THC rank = {norb*full})")
        print(f"  recon err = {np.linalg.norm(eri_thc-eri):.2e}   "
              f"FCI(eri_THC) err = {abs(e_thc-cas.e_tot)*1e3:.2e} mHa")
        print(f"  lambda: thc(DF)= {thc_lambda(chi_df, zeta_df, h1):.3f}  "
              f"df= {df_lambda(leaves, h1, norb):.3f}  (equal: validates the formula)")
        print(f"  lambda: thc(random collocation)= {thc_lambda(chi, zeta, h1):.1f}  "
              f">> df -- unoptimized collocation; the THC lambda win needs ISDF points.")
    print("=" * 78)
