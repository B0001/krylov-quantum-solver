#!/usr/bin/env python3
"""
Nb3X8 bilayer-cluster charge gaps: exact diagonalization vs the Hubbard-I approximation.

The breathing-mode kagome van der Waals materials Nb3X8 (X = F, Cl, Br, I) downfold, per bilayer, to
a *generalized Hubbard dimer* of two inter-layer-dimerized trimer molecular orbitals (Aretz, Grytsiuk,
Liu, ..., van Loon, Rosner, arXiv:2501.10320). The source paper solves this cluster with cluster
dynamical mean-field theory in the Hubbard-I approximation to obtain the correlated gap. That cluster
is only two orbitals (four spin-orbitals), so it is *exactly diagonalizable* -- and the exact charge
gap is a number the paper did not report.

This study computes, from the paper's own ab-initio cRPA parameters (Table I, LT bulk):
  * the **exact** cluster charge gap  Delta = E(N+1) + E(N-1) - 2 E(N)  at half-filling (N=2), and
  * the **Hubbard-I** gap for the same cluster (the atomic self-energy embedded in the dimer
    dispersion),
and quantifies the Hubbard-I error across the correlation-tuned family.

Finding (see specs/SPEC_nb3x8_gaps.md): Hubbard-I is essentially exact for the strongly-correlated
members (Nb3F8, Nb3Cl8) but **underestimates the charge gap of the weakly-correlated Nb3I8 by ~29%
(~244 meV)**; the error grows monotonically as U0/|t| falls -- exactly where an atomic-limit
self-energy must fail. Both methods agree (-> U0) in the atomic limit t -> 0, which validates the
machinery.

HONEST SCOPE: this is the *isolated* two-orbital cluster. The paper's cluster-DMFT embeds it in a
self-consistent bath, so this quantifies the impurity-solver (Hubbard-I) error *on the cluster*, not
the solid's true gap. It is not a claim that the paper's material gaps are wrong by 29% -- it flags
which cluster-level approximation is least reliable. Density-density interactions only (the paper
reports non-density-density terms of only a few meV). Gaps are in meV.
"""
from __future__ import annotations

import numpy as np

from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals, fixed_filling_energy

# LT bulk ab-initio downfolded parameters (meV) -- Table I of arXiv:2501.10320.
# U0: on-site (per trimer-MO) Hubbard U; t: strong inter-layer hopping t_s_perp;
# Us: inter-site (inter-layer) density-density Coulomb U_s_perp.
NB3X8_LT_BULK = {
    "Nb3F8":  dict(U0=2590.5, t=-4.9,   Us=714.6),
    "Nb3Cl8": dict(U0=1451.4, t=-136.0, Us=400.1),
    "Nb3Br8": dict(U0=1186.6, t=-169.4, Us=342.0),
    "Nb3I8":  dict(U0=787.0,  t=-218.2, Us=258.5),
}


def dimer_cluster_integrals(U0: float, t: float, Us: float) -> ModelIntegrals:
    """The generalized-Hubbard-dimer cluster: two trimer orbitals, on-site ``U0``, inter-layer
    hopping ``t``, inter-site density-density ``Us``. Half-filled (2 electrons)."""
    h1 = np.array([[0.0, t], [t, 0.0]])
    eri = np.zeros((2, 2, 2, 2))
    eri[0, 0, 0, 0] = eri[1, 1, 1, 1] = U0        # on-site Hubbard
    eri[0, 0, 1, 1] = eri[1, 1, 0, 0] = Us        # inter-site density-density
    return ModelIntegrals(h1=h1, eri=eri, e_core=0.0, nelec=(1, 1), norb=2)


def exact_charge_gap(U0: float, t: float, Us: float) -> float:
    """Exact charge (Mott) gap ``E(3) + E(1) - 2 E(2)`` of the cluster by full diagonalization."""
    cluster = dimer_cluster_integrals(U0, t, Us)

    def E(n):
        return fixed_filling_energy(ModelIntegrals(cluster.h1, cluster.eri, 0.0,
                                                   (n - n // 2, n // 2), 2))

    return E(3) + E(1) - 2 * E(2)


def hubbard_i_gap(U0: float, t: float, Us: float) -> float:
    """Hubbard-I charge gap of the dimer cluster.

    The atomic self-energy of a half-filled correlated orbital, ``Sigma_at(w) = w - 2x(x-U0)/(2x-U0)``
    with ``x = w - h`` (h the inter-site Hartree shift), embedded in the bonding/anti-bonding dimer
    dispersion gives poles solving ``w - Sigma_at(w) = +/- t``, i.e. the quadratics
    ``2x^2 - (2U0 +/- 2t) x +/- t U0 = 0``. The gap is the spacing between the two occupied and two
    unoccupied poles (h-independent). Exact in the atomic limit ``t -> 0`` (gap -> U0).
    """
    poles = []
    for s in (+1.0, -1.0):
        poles += list(np.roots([2.0, -(2 * U0 + 2 * s * t), s * t * U0]))
    poles = np.sort(np.real(poles))
    return float(poles[2] - poles[1])


if __name__ == "__main__":
    print("Nb3X8 bilayer-cluster charge gaps (meV), from arXiv:2501.10320 cRPA parameters")
    print(f"{'compound':8} {'U0/|t|':>7} {'exact':>8} {'Hubbard-I':>10} {'HubI error':>11} {'%':>7}")
    for name, p in NB3X8_LT_BULK.items():
        ge, gh = exact_charge_gap(**p), hubbard_i_gap(**p)
        print(f"{name:8} {p['U0']/abs(p['t']):7.1f} {ge:8.1f} {gh:10.1f} "
              f"{gh - ge:11.1f} {(gh - ge) / ge * 100:6.1f}%")
    print("\nFinding: Hubbard-I is exact for the strongly-correlated F/Cl but underestimates the "
          "weakly-\ncorrelated Nb3I8 gap by ~29% -- the atomic self-energy fails where hybridization "
          "competes with U.")
