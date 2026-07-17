#!/usr/bin/env python3
"""
The lambda ladder: how Hamiltonian factorization changes the qubitization 1-norm lambda,
which sets the FT-QPE T-gate budget (walk steps ~ O(lambda / epsilon)).

Three representations of the same active-space Hamiltonian, compared on lambda, the LCU term
count, and the energy error from truncation:

  1. NAIVE   -- full Pauli LCU. lambda = sum of |Pauli coefficients|. Exact, but no structure.
  2. DF       -- double factorization, truncated to rank R. Trades accuracy for a smaller
                 lambda and far fewer factors (the H2O example: 2401 integrals -> ~21 factors).
  3. THC      -- tensor hypercontraction, (pq|rs) ~ sum_uv X_pu X_qu Z_uv X_rv X_sv.

HONEST CAVEAT (this tool will show it): THC's lambda advantage over DF is ASYMPTOTIC -- it
comes from THC needing O(N) factors where DF needs O(N^2), so it only appears for large
active spaces (dozens of orbitals). On a small CAS, a generic THC fit is comparable to or
denser than DF; do not expect a small-system win. Likewise, the brute-force Pauli lambda used
here is exact but only computable for <= ~4 orbitals; large-N lambda must come from the
factorization-native 1-norm formulas (a careful next step, not done here).

Requires qubitization_blueprint.py and df_factorization.py in the same directory.
"""

import os
import sys
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qubitization_blueprint import build_qubit_hamiltonian, pauli_decompose
from df_factorization import double_factorize, reconstruct_eri


def lambda_and_terms(h1, eri, norb):
    """Exact qubitization 1-norm (brute-force Pauli). Feasible only for small norb (<= ~4)."""
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    terms = pauli_decompose(H, n)
    return sum(abs(c) for _, c in terms), len(terms)


def fci_energy_error(h1, eri, norb, nelec, e_core, casci_energy):
    from pyscf import fci
    e, _ = fci.direct_spin1.kernel(h1, eri, norb, nelec)
    return abs((e + e_core) - casci_energy) * 1e3


def fit_thc(eri, norb, M, restarts=4, seed=0):
    """Least-squares THC fit: (pq|rs) ~ sum_uv X_pu X_qu Z_uv X_rv X_sv, symmetric Z."""
    rng = np.random.default_rng(seed)
    target = eri.reshape(-1)
    iu = np.triu_indices(M)
    best = None
    for _ in range(restarts):
        X0 = rng.standard_normal((norb, M)) * 0.3
        Z0 = rng.standard_normal((M, M)) * 0.1
        Z0 = 0.5 * (Z0 + Z0.T)
        p0 = np.concatenate([X0.ravel(), Z0[iu]])

        def unpack(p):
            X = p[:norb * M].reshape(norb, M)
            Zf = np.zeros((M, M))
            Zf[iu] = p[norb * M:]
            Zf = Zf + Zf.T - np.diag(np.diag(Zf))
            return X, Zf

        def resid(p):
            X, Z = unpack(p)
            return np.einsum("pu,qu,uv,rv,sv->pqrs", X, X, Z, X, X).reshape(-1) - target

        sol = least_squares(resid, p0, method="lm", max_nfev=4000)
        if best is None or sol.cost < best[0]:
            X, Z = unpack(sol.x)
            best = (sol.cost, np.einsum("pu,qu,uv,rv,sv->pqrs", X, X, Z, X, X))
    return best[1]


def lambda_ladder(h1, eri, norb, nelec, e_core, casci_energy, thc_ranks=range(2, 7)):
    """Print the naive / DF / THC comparison for one active space."""
    lam0, L0 = lambda_and_terms(h1, eri, norb)
    print("=" * 74)
    print(f"active space norb={norb} nelec={nelec}  CASCI={casci_energy:.8f}  N^4={norb**4}")
    print(f"[NAIVE] full Pauli LCU: lambda={lam0:.4f}  terms={L0}  (exact)")
    print("-" * 74)

    _, _, full_rank = double_factorize(eri, norb)
    print(f"[DF] exact rank={full_rank}")
    print(f"{'rank':>5} {'lambda':>9} {'terms':>7} {'err_mHa':>12}")
    for R in range(1, full_rank + 1):
        leaves, _, _ = double_factorize(eri, norb, rank=R)
        eriR = reconstruct_eri(leaves, norb)
        lamR, LR = lambda_and_terms(h1, eriR, norb)
        print(f"{R:>5} {lamR:>9.4f} {LR:>7} "
              f"{fci_energy_error(h1, eriR, norb, nelec, e_core, casci_energy):>12.4f}")
    print("-" * 74)

    print("[THC] X X Z X X fit (advantage is asymptotic; small-CAS win not expected)")
    print(f"{'rankM':>5} {'lambda':>9} {'terms':>7} {'err_mHa':>12}")
    for M in thc_ranks:
        eriT = fit_thc(eri, norb, M)
        lamT, LT = lambda_and_terms(h1, eriT, norb)
        print(f"{M:>5} {lamT:>9.4f} {LT:>7} "
              f"{fci_energy_error(h1, eriT, norb, nelec, e_core, casci_energy):>12.4f}")
    print("=" * 74)


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo

    mol = gto.M(atom="N 0 0 0; N 0 0 1.10", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    norb, ne = 3, 4
    cas = mcscf.CASCI(mf, norb, ne)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    lambda_ladder(h1, eri, norb, (ne // 2, ne // 2), e_core, cas.e_tot)
    print("Reading: DF lowers lambda at an accuracy cost (rank-truncation tradeoff). THC matches")
    print("DF here and is denser -- its O(N)-factor advantage needs a large active space to show.")
