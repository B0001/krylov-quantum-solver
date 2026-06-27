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
