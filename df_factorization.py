#!/usr/bin/env python3
"""
Double factorization (DF) of the ERI tensor -- the highest-leverage knob on the
fault-tolerant T-gate budget.

The two-electron integrals (pq|rs) form a positive semidefinite supermatrix in the
composite index (pq)x(rs). Two nested eigendecompositions expose low-rank structure:

  1. First factorization:  (pq|rs) = sum_t g_t L^{(t)}_pq L^{(t)}_rs
  2. Second factorization:  L^{(t)} = U^{(t)} diag(w^{(t)}) U^{(t)T}   (rotated number operators)

Truncating the rank R (number of retained g_t) replaces the O(N^4) integral list with a
handful of tensor factors at a controlled energy cost. In a qubitization block encoding this
shrinks both the number of LCU terms and the 1-norm lambda, and therefore the walk-step /
T-gate count of FT-QPE. The second factorization is what gives each term an efficient circuit
form (an orbital rotation around a layer of number operators).

Same active-space inputs as the SQD / Krylov / qubitization modules: eri from cas.get_h2eff().
"""

import numpy as np


def double_factorize(eri, norb, rank=None, tol=1e-10):
    """Return (leaves, g, full_rank).

    leaves: list of (g_t, w^{(t)}, U^{(t)}) ordered by |g_t| descending, truncated to `rank`.
    g: all first-factorization eigenvalues (magnitude-sorted).
    full_rank: number of factors above `tol` (the exact rank).
    """
    V = eri.reshape(norb * norb, norb * norb)
    V = 0.5 * (V + V.T)
    g, U = np.linalg.eigh(V)
    order = np.argsort(np.abs(g))[::-1]
    g = g[order]
    U = U[:, order]
    full_rank = int(np.sum(np.abs(g) > tol))
    R = full_rank if rank is None else min(rank, full_rank)
    leaves = []
    for t in range(R):
        L = U[:, t].reshape(norb, norb)
        L = 0.5 * (L + L.T)
        w, Uk = np.linalg.eigh(L)        # second factorization
        leaves.append((float(g[t]), w, Uk))
    return leaves, g, full_rank


def reconstruct_eri(leaves, norb):
    """Rebuild the (truncated) ERI tensor from its double-factorized leaves."""
    eri = np.zeros((norb,) * 4)
    for gt, w, Uk in leaves:
        L = (Uk * w) @ Uk.T
        eri += gt * np.einsum("pq,rs->pqrs", L, L)
    return eri


def df_lambda(leaves, h1, norb=None):
    """Factorization-native double-factorization 1-norm (no brute-force Pauli enumeration).

        lambda_DF = 1/4 * sum_t |g_t| (sum_k |w^(t)_k|)^2  +  ||h1||_nuc

    The two-body term is the qubitization 1-norm of the double-factorized two-body operator
    (von Burg et al. 2021; Lee et al., PRX Quantum 2, 030305, 2021): each leaf contributes a
    rotated layer of number operators whose coefficients are sqrt(|g_t|) * w^(t), so its 1-norm
    contribution is 1/4 |g_t| (sum_k |w_k|)^2. The one-body term is the nuclear norm of h1
    (sum of |eigenvalues|), the optimal 1-norm of a one-body operator. lambda_DF is the
    FT-QPE T-gate cost driver (walk steps ~ O(lambda / epsilon)).

    Unlike ``lambda_ladder.lambda_and_terms`` (exact Pauli LCU, feasible only for <= ~4
    orbitals), this is computed straight from the factorization tensors and scales to dozens of
    orbitals. It is a *different and smaller* 1-norm than the naive Pauli LCU -- that reduction
    is the whole point of double factorization, so the Pauli value is an upper reference, not an
    equality oracle.

    Args:
        leaves: ``double_factorize`` output -- list of (g_t, w^(t), U^(t)).
        h1: one-body integrals, shape (norb, norb).
        norb: unused (accepted for a uniform call signature with the other module functions).
    """
    two_body = 0.25 * sum(abs(gt) * float(np.abs(w).sum()) ** 2 for gt, w, _ in leaves)
    h1 = np.asarray(h1, dtype=float)
    one_body = float(np.linalg.svd(0.5 * (h1 + h1.T), compute_uv=False).sum())
    return two_body + one_body


def _apply_number_shift(h1, eri, norb, n_elec, b1, b2):
    """Apply the number-operator symmetry shift H -> H + (b1 N_e + b2)(N_e - n_e).

    On the fixed-electron-number sector (N_e - n_e)|psi> = 0, so the ground state and its energy
    are unchanged -- only the *integrals* (and hence lambda) move. Returns
    ``(h1_s, eri_s, e_shift)`` where ``FCI(h1_s, eri_s) + e_shift == FCI(h1, eri)``.

    Bookkeeping (PySCF's normal-ordered 2-RDM convention, where the two-body operator is
    1/2 sum eri_pqrs (E_pq E_rs - delta_qr E_ps)):
      * adding 2*b1 to eri[p,p,r,r] realises  b1 (N_e^2 - N_e);
      * the residual one-body shift to reach b1 N_e^2 + (b2 - b1 n_e) N_e is therefore
        h1[p,p] += b2 + b1 (1 - n_e);
      * the leftover scalar -b2 n_e is returned as ``e_shift``.
    """
    h1_s = np.array(h1, dtype=float)
    eri_s = np.array(eri, dtype=float)
    diag = np.arange(norb)
    eri_s[diag[:, None], diag[:, None], diag[None, :], diag[None, :]] += 2.0 * b1
    h1_s[diag, diag] += b2 + b1 * (1.0 - n_elec)
    return h1_s, eri_s, -b2 * float(n_elec)


def symmetry_shift(h1, eri, norb, nelec, b1=None, b2=None):
    """Number-operator (BLISS / SCDF-style) symmetry shift that lowers lambda_DF.

    Replaces H with H + (b1 N_e + b2)(N_e - n_e). N_e commutes with H and equals n_e on the
    target sector, so the spectrum on that sector is preserved exactly (Loaiza & Izmaylov 2023;
    Rocca et al., arXiv:2403.03502; Deka & Zak, arXiv:2412.01338). (b1, b2) are chosen to
    minimise lambda_DF; pass them explicitly to skip the optimisation.

    Returns ``(h1_s, eri_s, e_shift, (b1, b2))`` with
    ``FCI(h1_s, eri_s) + e_shift == FCI(h1, eri)``.

    This implements the number-operator (N_e, N_e^2) shift only -- not the full tensor-optimised
    SCDF cost function. See specs/SPEC_scdf_lambda.md.
    """
    n_elec = int(nelec) if isinstance(nelec, (int, np.integer)) else int(sum(nelec))

    if b1 is None or b2 is None:
        from scipy.optimize import minimize

        def objective(params):
            h1_s, eri_s, _ = _apply_number_shift(h1, eri, norb, n_elec, params[0], params[1])
            leaves, _, _ = double_factorize(eri_s, norb)
            return df_lambda(leaves, h1_s, norb)

        # Initialise the one-body shift at the median eigenvalue of h1 (the nuclear-norm
        # optimum for b1 = 0); let Nelder-Mead refine both b1 and b2.
        median_h1 = float(np.median(np.linalg.eigvalsh(0.5 * (np.asarray(h1) + np.asarray(h1).T))))
        res = minimize(objective, x0=[0.0, -median_h1], method="Nelder-Mead",
                       options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 400})
        b1, b2 = float(res.x[0]), float(res.x[1])

    h1_s, eri_s, e_shift = _apply_number_shift(h1, eri, norb, n_elec, b1, b2)
    return h1_s, eri_s, e_shift, (b1, b2)


def rank_for_accuracy(h1, eri, norb, nelec, e_core, casci_energy, target_mHa=1.6):
    """Smallest retained rank whose FCI energy is within target_mHa of CASCI."""
    from pyscf import fci
    _, _, full_rank = double_factorize(eri, norb)
    for R in range(1, full_rank + 1):
        leaves, _, _ = double_factorize(eri, norb, rank=R)
        e_R, _ = fci.direct_spin1.kernel(h1, reconstruct_eri(leaves, norb), norb, nelec)
        if abs((e_R + e_core) - casci_energy) * 1e3 <= target_mHa:
            return R, full_rank
    return full_rank, full_rank


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo, fci

    mol = gto.M(atom="O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    norb = mol.nao_nr()
    na = nb = mol.nelectron // 2
    cas = mcscf.CASCI(mf, norb, (na, nb))
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)

    _, g, full_rank = double_factorize(eri, norb)
    print("=" * 72)
    print(f"H2O STO-3G  norb={norb}  nelec=({na},{nb})  FCI={cas.e_tot:.8f} Ha")
    print(f"ERI integrals O(N^4) = {norb**4};  exact factorization rank = {full_rank}")
    print(f"top factor magnitudes |g_t|: {np.round(np.abs(g[:8]), 4)}")
    print("-" * 72)
    print(f"{'rank R':>7} {'recon_err':>12} {'energy_err_mHa':>16}")
    for R in range(1, full_rank + 1):
        leaves, _, _ = double_factorize(eri, norb, rank=R)
        eri_R = reconstruct_eri(leaves, norb)
        e_R, _ = fci.direct_spin1.kernel(h1, eri_R, norb, (na, nb))
        derr = abs((e_R + e_core) - cas.e_tot) * 1e3
        print(f"{R:>7} {np.linalg.norm(eri_R - eri):>12.2e} {derr:>16.4f}")
    R_ca, _ = rank_for_accuracy(h1, eri, norb, (na, nb), e_core, cas.e_tot)
    print("-" * 72)
    print(f"Chemical accuracy (<1.6 mHa) at rank {R_ca} of {full_rank}: "
          f"{norb**4} integrals -> {R_ca} tensor factors. That compression is the lambda / T-gate win.")
    print("=" * 72)
