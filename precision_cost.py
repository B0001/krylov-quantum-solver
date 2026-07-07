#!/usr/bin/env python3
"""
The precision-cost crossover -- near-term certified (1/eps^2) vs fault-tolerant QPE (1/eps).

This bridges the repo's two halves. The certified arc (`temple_bounds` ... `certified_noise`) is a
near-term, sampling method: to certify the ground energy to precision eps it needs
N(eps) = (z lambda_meas / eps)^2 shot-measurements -- the STANDARD quantum limit, 1/eps^2, whose
exponent `certified_noise` measured directly (shot-cost slope -2). Fault-tolerant qubitization/QPE
needs Q(eps) = pi lambda_DF / (2 eps) queries -- the HEISENBERG limit, 1/eps. So the resource ratio
R = N/Q ~ 1/eps: FT wins asymptotically at high precision, by an exponent, not a constant.

THE FINDINGS (specs/SPEC_precision_cost.md), on real molecular 1-norms (H2, H2O, N2 CAS(6,6)):
  * EXPONENTS: near-term shots slope -2, FT queries slope -1, ratio slope -1 (the standard-vs-
    Heisenberg gap). At chemical accuracy (1.6 mHa) R ~ 4e3 - 3e4 -- FT is far cheaper in resource
    COUNT, and the gap widens as eps -> 0.
  * THE FT WIN IS THE EXPONENT, NOT THE CONSTANT. The RAW qubitization lambda_DF does NOT uniformly
    beat the measurement 1-norm lambda_meas: for N2 CAS(6,6) lambda_DF = 24.94 > lambda_meas = 22.84.
    Double factorization alone does not shrink the FT constant below measurement.
  * THE SYMMETRY SHIFT EARNS THE CONSTANT. The number-operator shift (`scdf_lambda`) drops
    lambda_DF to 0.97 / 1.83 / 4.00 -- below lambda_meas for EVERY molecule, by a margin that GROWS
    with system size (2.8x H2 -> 5.7x N2). So FT wins the exponent unconditionally and the constant
    only after the shift; the shift is not a nicety but load-bearing for the FT resource advantage.

HONEST SCOPE: comparing shot-measurements to qubitization queries needs a common cost model -- the
per-query T-gate cost of the block encoding is NOT computed here (that is `ft_resource_estimator`,
openfermion, the chem-ft env). The unit-INDEPENDENT claim is the exponent gap (R ~ 1/eps); the
crossover LOCATION is parametrized by the per-query/per-shot cost ratio. The near-term 1/eps^2 uses
`certified_noise`'s i.i.d. shot model (z*lambda/sqrt(N)); the FT 1/eps uses the standard QPE query
count. Reproduction-adjacent on the scaling laws (both are textbook limits); the composition into a
per-molecule certified-vs-FT crossover on real, shifted lambdas is the contribution.
"""
from __future__ import annotations

import numpy as np

from df_factorization import df_lambda, double_factorize, symmetry_shift
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian


def measurement_lambda(mh: MolecularHamiltonian) -> float:
    """lambda_meas = sum |Pauli coefficients| of the qubit Hamiltonian (the shot-noise 1-norm)."""
    return float(np.abs(mh.qubit_hamiltonian.coeffs).sum())


def qubitization_lambda(h1, eri, norb: int, nelec=None, shift: bool = False) -> float:
    """lambda_DF: the double-factorized qubitization 1-norm; symmetry-shifted (SCDF) if ``shift``."""
    if shift:
        if nelec is None:
            raise ValueError("nelec required for the symmetry shift")
        h1, eri, _, _ = symmetry_shift(h1, eri, norb, nelec)
    return float(df_lambda(double_factorize(eri, norb)[0], h1, norb))


def near_term_shots(lam_meas: float, eps: float, z: float = 2.0) -> float:
    """Shot-measurements to certify the energy to +/- eps (standard limit, from certified_noise):
    N = (z lambda_meas / eps)^2."""
    return (z * lam_meas / eps) ** 2


def ft_queries(lam_df: float, eps: float) -> float:
    """Qubitization/QPE queries to resolve the energy to eps (Heisenberg limit): Q = pi lambda_DF /
    (2 eps)."""
    return np.pi * lam_df / (2.0 * eps)


def resource_ratio(lam_meas: float, lam_df: float, eps: float, z: float = 2.0) -> float:
    """N(eps) / Q(eps) -- unitless resource COUNT ratio (near-term / FT). ~ 1/eps."""
    return near_term_shots(lam_meas, eps, z) / ft_queries(lam_df, eps)


def crossover_epsilon(lam_meas: float, lam_df: float, cost_per_query: float = 1.0,
                      z: float = 2.0) -> float:
    """Precision eps* at which the common-model costs cross: below eps* FT is cheaper. With FT cost
    = Q * cost_per_query and near-term cost = N, eps* = 2 z^2 lambda_meas^2 / (pi cost_per_query
    lambda_DF)."""
    return 2.0 * z ** 2 * lam_meas ** 2 / (np.pi * cost_per_query * lam_df)


if __name__ == "__main__":
    from pyscf import ao2mo, gto, mcscf, scf

    from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals

    def cas(atom, norb, nelec):
        mf = scf.RHF(gto.M(atom=atom, basis="sto-3g", verbose=0)).run()
        c = mcscf.CASCI(mf, norb, nelec)
        c.kernel()
        h1, ec = c.get_h1eff()
        eri = ao2mo.restore(1, c.get_h2eff(), norb)
        na = (nelec + nelec % 2) // 2
        return np.asarray(h1), np.asarray(eri), float(ec), (na, nelec - na), norb

    mols = {"H2": ("H 0 0 0; H 0 0 0.74", 2, 2),
            "H2O(4,3)": ("O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467", 3, 4),
            "N2(6,6)": ("N 0 0 0; N 0 0 1.10", 6, 6)}
    print("Precision-cost crossover: near-term certified (1/eps^2) vs FT-QPE (1/eps)")
    print(f"{'mol':9s} | {'lam_meas':>8} | {'lam_DFraw':>9} | {'lam_DFshift':>11} | "
          f"{'R @1.6mHa':>10} | raw<meas | shift<meas")
    for name, (atom, norb, nel) in mols.items():
        h1, eri, ec, ne, no = cas(atom, norb, nel)
        mh = build_hamiltonian_from_integrals(h1, eri, ne, ec)
        lm = measurement_lambda(mh)
        lr = qubitization_lambda(h1, eri, no)
        ls = qubitization_lambda(h1, eri, no, nelec=ne, shift=True)
        R = resource_ratio(lm, ls, 1.6e-3)
        print(f"{name:9s} | {lm:8.2f} | {lr:9.2f} | {ls:11.2f} | {R:10.2e} | "
              f"{str(lr < lm):>8} | {str(ls < lm):>9}")
    print("\nFT wins the EXPONENT unconditionally (R ~ 1/eps); raw lambda_DF can exceed lambda_meas")
    print("(N2), so the symmetry shift is load-bearing for the FT constant. See SPEC_precision_cost.")
