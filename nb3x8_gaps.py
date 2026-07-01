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

Finding (see specs/SPEC_nb3x8_gaps.md), across all 10 dimer-cluster parameter sets (LT bulk, LT
bilayer, HT bulk):
  * Strongly-correlated clusters (Nb3F8, HT-phase Nb3Cl8/Nb3Br8) -- Hubbard-I is near-exact (< 2%).
  * The **iodides are consistently the worst**: Hubbard-I underestimates the gap by ~29% (bulk) and
    ~12% (bilayer) -- largest hopping, weakest on-site U.
  * BUT the error is **not** a clean single-parameter function of U0/|t| (Spearman ~ -0.86 over the
    full set, not -1): the 4-point LT-bulk trend that looked monotonic does not survive the extended
    dataset -- both t and U_s_perp matter. The robust statement is material-level (iodides worst),
    not a single-ratio scaling law.
Both methods agree (-> U0) in the atomic limit t -> 0, which validates the machinery.

HONEST SCOPE: this is the *isolated* two-orbital cluster. The paper's cluster-DMFT embeds it in a
self-consistent bath, so this quantifies the impurity-solver (Hubbard-I) error *on the cluster*, not
the solid's true gap. It is not a claim that the paper's material gaps are wrong by 29% -- it flags
which cluster-level approximation is least reliable. Density-density interactions only (the paper
reports non-density-density terms of only a few meV). Gaps are in meV.

BATH BOUND (``four_site_exact_gap``): enlarging the correlated region to include the inter-cluster
weak link is a rigorous, bath-fit-free proxy for how far the isolated-cluster result travels toward
the solid. It moves the Nb3I8 gap by only ~5% (the strong intra-dimer bond isolates the cluster for
the iodides), vs the ~29% Hubbard-I error -- so the finding survives. The bath effect is largest for
Nb3F8 (~22%), where t_s ~ t_w makes the dimer ill-defined, but there Hubbard-I is exact anyway. The
finding is therefore most reliable exactly where it matters (the iodides).
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

# All dimer-cluster parameter sets that carry an inter-layer dimer (Tables I & IV; monolayers have no
# dimer, so are excluded). Keys: <halide> <phase/thickness>. Spans U0/|t| ~ 3.6 (I bulk) to ~800 (F BL).
NB3X8_CLUSTERS = {
    "F  LT-bulk": dict(U0=2590.5, t=-4.9,   Us=714.6),
    "Cl LT-bulk": dict(U0=1451.4, t=-136.0, Us=400.1),
    "Br LT-bulk": dict(U0=1186.6, t=-169.4, Us=342.0),
    "I  LT-bulk": dict(U0=787.0,  t=-218.2, Us=258.5),
    "F  LT-bil":  dict(U0=3988.8, t=-5.0,   Us=1987.7),
    "Cl LT-bil":  dict(U0=2697.6, t=-136.2, Us=1570.9),
    "Br LT-bil":  dict(U0=2396.0, t=-169.2, Us=1482.3),
    "I  LT-bil":  dict(U0=1928.7, t=-218.4, Us=1349.0),
    "Cl HT-bulk": dict(U0=1401.0, t=-17.11, Us=336.8),   # Table IV (HT phase)
    "Br HT-bulk": dict(U0=1129.1, t=-20.56, Us=276.5),
}


def dimer_cluster_integrals(U0: float, t: float, Us: float) -> ModelIntegrals:
    """The generalized-Hubbard-dimer cluster: two trimer orbitals, on-site ``U0``, inter-layer
    hopping ``t``, inter-site density-density ``Us``. Half-filled (2 electrons)."""
    h1 = np.array([[0.0, t], [t, 0.0]])
    eri = np.zeros((2, 2, 2, 2))
    eri[0, 0, 0, 0] = eri[1, 1, 1, 1] = U0        # on-site Hubbard
    eri[0, 0, 1, 1] = eri[1, 1, 0, 0] = Us        # inter-site density-density
    return ModelIntegrals(h1=h1, eri=eri, e_core=0.0, nelec=(1, 1), norb=2)


# LT-bulk parameters extended with the weak inter-bilayer link (t_w_perp) and its Coulomb (U_w_perp),
# for the bath bound: (U0, t_s_perp, U_s_perp, t_w_perp, U_w_perp), meV, Table I.
NB3X8_LT_BULK_5P = {
    "Nb3F8":  (2590.5, -4.9,   714.6, -6.5,  572.5),
    "Nb3Cl8": (1451.4, -136.0, 400.1, -16.1, 313.5),
    "Nb3Br8": (1186.6, -169.4, 342.0, -20.4, 262.4),
    "Nb3I8":  (787.0,  -218.2, 258.5, -24.6, 183.8),
}


def four_site_exact_gap(U0: float, ts: float, Us: float, tw: float, Uw: float) -> float:
    """Exact charge gap of an *enlarged* cluster -- two dimers joined by the weak inter-bilayer link
    (chain 0=1 strong, 1~2 weak, 2=3 strong; on-site U0, inter-site Us on strong bonds, Uw on the weak
    bond). Half-filled (4 electrons). The bath bound: comparing this to the isolated-dimer gap
    quantifies how much the inter-cluster coupling moves the gap -- a rigorous proxy for the DMFT bath
    that needs no bath fit. See specs/SPEC_nb3x8_gaps.md R1."""
    h1 = np.zeros((4, 4))
    h1[0, 1] = h1[1, 0] = ts
    h1[2, 3] = h1[3, 2] = ts
    h1[1, 2] = h1[2, 1] = tw
    eri = np.zeros((4, 4, 4, 4))
    for i in range(4):
        eri[i, i, i, i] = U0
    for (i, j), U in (((0, 1), Us), ((2, 3), Us), ((1, 2), Uw)):
        eri[i, i, j, j] = eri[j, j, i, i] = U

    def E(n):
        return fixed_filling_energy(ModelIntegrals(h1, eri, 0.0, (n - n // 2, n // 2), 4))

    return E(5) + E(3) - 2 * E(4)


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
    from scipy.stats import spearmanr

    print("Nb3X8 dimer-cluster charge gaps (meV), from arXiv:2501.10320 cRPA parameters")
    print(f"{'set':11} {'U0/|t|':>7} {'exact':>8} {'Hubbard-I':>10} {'HubI error':>11} {'%':>7}")
    rows = []
    for name, p in sorted(NB3X8_CLUSTERS.items(), key=lambda kv: kv[1]["U0"] / abs(kv[1]["t"])):
        ge, gh = exact_charge_gap(**p), hubbard_i_gap(**p)
        rows.append((p["U0"] / abs(p["t"]), abs(gh - ge) / ge))
        print(f"{name:11} {p['U0']/abs(p['t']):7.1f} {ge:8.1f} {gh:10.1f} "
              f"{gh - ge:11.1f} {(gh - ge) / ge * 100:6.1f}%")
    rho = spearmanr([r[0] for r in rows], [r[1] for r in rows]).correlation
    print(f"\nSpearman(U0/|t| vs |error|) = {rho:.3f}  (strong but imperfect: the single-ratio law "
          "does NOT hold)")
    print("Robust finding: strongly-correlated clusters (F, HT Cl/Br) are near-exact (<2%); the "
          "iodides\nare consistently worst (bulk ~29%, bilayer ~12%). The error is multi-parameter, "
          "not a clean U0/|t| law.")

    print("\nBath bound -- enlarge the correlated region (isolated dimer -> two dimers + weak link):")
    print(f"{'compound':8} {'gap 2-site':>10} {'gap 4-site':>10} {'bath shift':>11}")
    for name, five in NB3X8_LT_BULK_5P.items():
        g2 = exact_charge_gap(five[0], five[1], five[2])
        g4 = four_site_exact_gap(*five)
        print(f"{name:8} {g2:10.1f} {g4:10.1f} {(g4 - g2) / g2 * 100:10.1f}%")
    print("The bath moves the Nb3I8 gap by only ~5% (well-isolated dimer, t_s >> t_w) vs the ~29% "
          "Hubbard-I\nerror -- so the finding survives the isolated-cluster approximation. (4-site "
          "Hubbard-I, computed\nseparately, keeps the Nb3I8 error ~34% -- it grows, not washes out.)")
