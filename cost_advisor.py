#!/usr/bin/env python3
"""
Cost advisor -- a per-molecule near-term-vs-FT verdict with a parametrized crossover.

Turns the precision-cost scaling law (`precision_cost`) into an engineering artifact: for a target
accuracy eps on a Hamiltonian, return which is cheaper -- near-term certified Krylov (N(eps) =
(z lambda_meas/eps)^2 shot-measurements, standard limit) or FT-QPE (Q(eps) = pi lambda_DF/(2 eps)
queries, Heisenberg limit) -- under a common cost model whose only free parameter is the honest
unknown rho = cost_per_query / cost_per_shot. Reports the crossover eps*(rho) and, cross-validated
over a rho-range, whether the verdict is robust to that unknown.

This is the repo's `validate_and_cost` / `cross_check` instinct applied to the certified-vs-FT
choice. The exponent gap (FT ~ 1/eps vs near-term ~ 1/eps^2) is robust; the crossover LOCATION moves
with rho, so rho is never hidden -- it is an explicit input and the recommendation states its
rho-sensitivity.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from precision_cost import (
    crossover_epsilon,
    ft_queries,
    measurement_lambda,
    near_term_shots,
    qubitization_lambda,
)


@dataclass
class Verdict:
    """A costed near-term-vs-FT recommendation at (eps, rho). Costs are resource counts."""
    eps: float
    rho: float              # per-query / per-shot cost ratio (the honest unknown)
    lam_meas: float
    lam_df: float
    shots: float            # near-term cost N(eps)
    queries: float          # FT query count Q(eps)
    eps_star: float         # crossover precision at this rho (below eps* -> FT cheaper)
    cheaper: str            # "FT" | "near-term" at (eps, rho)
    rho_flip: float         # rho at which the verdict flips at this eps (= N/Q)


def advise(lam_meas: float, lam_df: float, eps: float, rho: float = 1.0, z: float = 2.0) -> Verdict:
    """Verdict at target accuracy ``eps`` and cost ratio ``rho``. FT cost = rho * queries."""
    shots = near_term_shots(lam_meas, eps, z)
    queries = ft_queries(lam_df, eps)
    cheaper = "FT" if rho * queries < shots else "near-term"
    return Verdict(eps=eps, rho=rho, lam_meas=lam_meas, lam_df=lam_df, shots=shots, queries=queries,
                   eps_star=crossover_epsilon(lam_meas, lam_df, cost_per_query=rho, z=z),
                   cheaper=cheaper, rho_flip=shots / queries)


def advise_from_integrals(h1, eri, e_core: float, nelec, norb: int, eps: float, rho: float = 1.0,
                          shift: bool = True, z: float = 2.0) -> Verdict:
    """Verdict from active-space integrals. ``shift`` applies the SCDF symmetry shift to lambda_DF."""
    from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals

    mh = build_hamiltonian_from_integrals(h1, eri, nelec, e_core)
    lam_meas = measurement_lambda(mh)
    lam_df = qubitization_lambda(h1, eri, norb, nelec=nelec, shift=shift)
    return advise(lam_meas, lam_df, eps, rho=rho, z=z)


def robust_over_rho(lam_meas: float, lam_df: float, eps: float, rho_lo: float, rho_hi: float,
                    z: float = 2.0):
    """(verdict, robust): is the cheaper method invariant across rho in [rho_lo, rho_hi]? The flip
    is monotone in rho (FT loses as rho grows), so checking the endpoints suffices."""
    lo = advise(lam_meas, lam_df, eps, rho=rho_lo, z=z).cheaper
    hi = advise(lam_meas, lam_df, eps, rho=rho_hi, z=z).cheaper
    return (lo, True) if lo == hi else ("mixed", False)


if __name__ == "__main__":
    from pyscf import ao2mo, gto, mcscf, scf

    def cas(atom, norb, nelec):
        mf = scf.RHF(gto.M(atom=atom, basis="sto-3g", verbose=0)).run()
        c = mcscf.CASCI(mf, norb, nelec)
        c.kernel()
        h1, ec = c.get_h1eff()
        return np.asarray(h1), ao2mo.restore(1, c.get_h2eff(), norb), float(ec), \
            ((nelec + nelec % 2) // 2, nelec - (nelec + nelec % 2) // 2), norb

    mols = {"H2": ("H 0 0 0; H 0 0 0.74", 2, 2),
            "N2(6,6)": ("N 0 0 0; N 0 0 1.10", 6, 6)}
    eps = 1.6e-3
    print(f"Cost-advisor verdict at chemical accuracy (eps = {eps} Ha):")
    print(f"{'mol':9s} | {'lam_meas':>8} | {'lam_DFshift':>11} | {'rho* (flip)':>11} | "
          f"robust FT over rho in [1, 1e6]?")
    for name, (atom, norb, nel) in mols.items():
        h1, eri, ec, ne, no = cas(atom, norb, nel)
        v = advise_from_integrals(h1, eri, ec, ne, no, eps, rho=1.0, shift=True)
        verdict, robust = robust_over_rho(v.lam_meas, v.lam_df, eps, 1.0, 1e6)
        print(f"{name:9s} | {v.lam_meas:8.2f} | {v.lam_df:11.2f} | {v.rho_flip:11.2e} | "
              f"{verdict} (robust={robust})")
    print("\nrho = per-query / per-shot cost (the honest unknown); FT wins until rho exceeds rho*.")
