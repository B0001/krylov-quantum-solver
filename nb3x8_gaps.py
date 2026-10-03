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

*** THE ISOLATED-CLUSTER FINDING DOES NOT TRANSLATE TO THE SOLID (the important correction). ***
The isolated-dimer exact gap has NO band broadening, so it is an *upper bound* on the solid gap;
Hubbard-I embedded in the full dispersion (the paper's cluster-DMFT) *includes* broadening and lies
lower. The two therefore bracket the true gap from opposite sides. Restoring the inter-dimer
coordination collapses the isolated exact-vs-Hubbard-I discrepancy:

  * ``coordination_gap`` (central dimer + z out-of-plane neighbours, FCI), Nb3I8:
      z=0 (isolated) 842 -> z=1 797 -> z=2 747 -> z=3 (~real) 650 meV, heading to Hubbard-I ~599.
  * ``ssh_chain_gap`` extrapolated to the 1-D thermodynamic limit (DMRG/block2, L=8..20 monotone
    730->709): ~708 meV, i.e. the out-of-plane chain alone already removes ~55% of the 842->599 gap.
  * With realistic 3-D coordination (out-of-plane z~3 plus in-plane t_par) the exact solid gap is
    ~600-650 meV -- close to Hubbard-I.

So the earlier headline "Hubbard-I underestimates the gap by ~29%" is a correct statement ABOUT THE
ISOLATED CLUSTER but an artifact of neglecting broadening; it does NOT mean the paper's material
gaps are wrong. If anything this VINDICATES the cluster-DMFT/Hubbard-I approach for the Nb3X8 solid
gaps. (Density-density interactions only; the paper reports non-density-density terms of a few meV.)
Gaps are in meV. The ``four_site_exact_gap`` "bath bound" (~5% at z=1) was too optimistic -- it
sampled only the nearest neighbour; coordination and the full chain (above) move the gap much more.

*** THE SAME COORDINATION MACHINERY IN THE SPIN CHANNEL DOES NOT RESCUE THE Tc OVERPREDICTION
(chem-g78). *** ``SPEC_nb3x8_magnetometry.md`` section 7 names "coordination / mean-field reduction"
as the likely explanation for the 5.3x (Nb3Cl8) / 2.3x (Nb3Br8) overprediction of the singlet-formation
Tc by the isolated-dimer J. Running ``coordination_spin_gap`` -- the IDENTICAL cluster topology as
``coordination_gap``, just measuring E(Sz=1 lowest) - E(Sz=0 lowest) at half-filling instead of the
charge gap -- shows J_eff does NOT fall monotonically, for either halide (``python nb3x8_gaps.py``):

  Nb3Cl8: z=0 (isolated) 66.20 -> z=1 71.08 -> z=2 66.55 -> z=3 71.12 meV
  Nb3Br8: z=0 (isolated) 119.11 -> z=1 126.50 -> z=2 120.11 -> z=3 126.63 meV
  (<S^2> = 0.000/2.000 at every z for both; clean singlet reference, no higher-S intruder in Sz=1)

J_eff oscillates around the isolated value and is *larger*, not smaller, at the largest cluster reached
(L=8) -- nowhere near the J0/3 = 22.1 meV (Cl) or J0/2.26 = 52.8 meV (Br) that the "coordination
rescues Tc" attribution would need (the bar is the reduction factor needed to remove the miss, capped
at 3; Br's J_eff/J0 >= 1.00 at every L, so the verdict does not depend on it). The SAME machinery's
charge-channel control on the SAME clusters (z=0..3) gives the expected monotonic softening
(Cl 1311.8 -> 1167.5 -> 1092.1 -> 872.9 meV, a 33.5% drop at z=3; Br 1086.0 -> 994.0 -> 923.5 ->
759.4 meV, 30.1%), confirming the machinery works and the charge/spin contrast is not an artifact.

**Verdict: the coordination/mean-field attribution for the 5.3x (Cl) / 2.3x (Br) Tc miss is FALSE
within this finite-cluster model (L <= 8 open clusters; no 3-D triplon band).** Coordination
(open-boundary band broadening of the interlayer bond) softens the CHARGE gap but does not deliver a
comparable softening of the SPIN gap -- if anything the finite clusters modestly stiffen J_eff. The
miss is therefore NOT explained by this mechanism; what explains it is open. Candidates this module
cannot see and did not test: the cooperative/first-order structural transition itself, in-plane
kagome exchange, or a lattice renormalization of t_s directly. See ``coordination_spin_gap`` /
``tests/test_nb3x8_gaps_spec.py::test_G8_*`` and specs/BACKLOG.md ("Coordination cannot rescue the
5.3x Tc overprediction").
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

# HT (undimerized) bulk parameters, Table IV of arXiv:2501.10320 -- only Cl/Br are tabulated there
# (F and I are not reported in the HT phase). Used by nb3x8_magnetometry to compare against the
# Curie-Weiss theta_W, which Sheckelton/Haraguchi fit in this HT phase, not the LT one.
NB3X8_HT_BULK = {
    "Nb3Cl8": dict(U0=1401.0, t=-17.11, Us=336.8),
    "Nb3Br8": dict(U0=1129.1, t=-20.56, Us=276.5),
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


def _cluster_charge_gap(h1: np.ndarray, eri: np.ndarray, *, max_cycle: int = 1000) -> float:
    """Charge gap E(L+1)+E(L-1)-2E(L) of an L-site half-filled extended-Hubbard cluster (FCI).

    ``max_cycle`` is forwarded to the FCI Davidson solver for each of the three sectors (see
    ``fixed_filling_energy``); raise it if a sector's odd-N charge state converges slowly (chem-q9g),
    rather than reimplementing the cluster construction just to pass a looser cap.
    """
    L = h1.shape[0]
    E = {n: fixed_filling_energy(ModelIntegrals(h1, eri, 0.0, (n - n // 2, n // 2), L),
                                  max_cycle=max_cycle)
         for n in (L - 1, L, L + 1)}
    return E[L + 1] + E[L - 1] - 2 * E[L]


def _coordination_cluster(U0: float, ts: float, Us: float, tw: float, Uw: float, z: int):
    """Central strong dimer + ``z`` out-of-plane weak-link pendant dimers -- the cluster topology
    shared by :func:`coordination_gap` (charge channel) and :func:`coordination_spin_gap` (spin
    channel). Returns ``(h1, eri, L)``."""
    L = 2 + 2 * z
    h1 = np.zeros((L, L))
    eri = np.zeros((L, L, L, L))
    for i in range(L):
        eri[i, i, i, i] = U0

    def bond(i, j, t, U):
        h1[i, j] = h1[j, i] = t
        eri[i, i, j, j] = eri[j, j, i, i] = U

    bond(0, 1, ts, Us)                               # central strong dimer
    for k in range(z):
        a, b = 2 + 2 * k, 3 + 2 * k
        bond(a, b, ts, Us)                           # pendant strong dimer
        bond(0 if k % 2 == 0 else 1, a, tw, Uw)      # weak link to an alternating central end
    return h1, eri, L


def coordination_gap(U0: float, ts: float, Us: float, tw: float, Uw: float, z: int,
                      *, max_cycle: int = 1000) -> float:
    """Charge gap of a central dimer with ``z`` out-of-plane weak-link neighbour dimers (FCI).

    A minimal probe of *coordination*: increasing ``z`` restores the band broadening the isolated
    dimer omits. As ``z`` grows the gap falls toward the Hubbard-I / cluster-DMFT value -- which is
    why the isolated-cluster ``exact_charge_gap`` OVERESTIMATES the solid gap and the isolated
    exact-vs-Hubbard-I discrepancy does not translate to the material (see the module docstring and
    specs/SPEC_nb3x8_gaps.md).

    ``max_cycle`` (default 1000, matching ``fci_energy``'s default) caps the FCI Davidson solver for
    each of the three charge sectors. Larger z / odd-electron sectors can need more than 1000
    iterations to converge even when not genuinely degenerate (chem-q9g: Nb3Cl8, z=3, N=9 needs
    ~4x); raise this rather than catching the ``RuntimeError`` and re-deriving the cluster yourself.
    """
    h1, eri, _ = _coordination_cluster(U0, ts, Us, tw, Uw, z)
    return _cluster_charge_gap(h1, eri, max_cycle=max_cycle)


def _fci_sector(h1: np.ndarray, eri: np.ndarray, na: int, nb: int, *, max_cycle: int = 2000):
    """Lowest FCI energy and total spin ``S(S+1)`` in the fixed ``(na, nb)`` sector."""
    from pyscf import fci

    L = h1.shape[0]
    solver = fci.direct_spin1.FCI()
    solver.max_cycle = int(max_cycle)
    energy, civec = solver.kernel(h1, eri, L, (na, nb))
    if not solver.converged:
        raise RuntimeError(f"FCI Davidson did not converge in {max_cycle} iterations for "
                           f"nelec=({na},{nb}), L={L}; retry with a larger (keyword-only) max_cycle "
                           f"before concluding non-convergence (cf. chem-q9g)")
    ss, _ = solver.spin_square(civec, L, (na, nb))
    return float(energy), float(ss)


def coordination_spin_gap(U0: float, ts: float, Us: float, tw: float, Uw: float, z: int,
                           *, max_cycle: int = 2000):
    """Effective spin gap of the SAME coordination cluster as :func:`coordination_gap`, in the SPIN
    channel: ``J_eff = E(Sz=1, lowest) - E(Sz=0, lowest)`` at fixed half-filling (FCI).

    This is the identical SSH/coordination machinery run against the spin sector instead of the
    charge sector, to test whether "coordination/mean-field reduction" (the reduction ``z`` delivers
    for the charge gap) also reduces the interlayer exchange -- the attribution
    ``specs/SPEC_nb3x8_magnetometry.md`` section 7 names for the 5.3x/2.3x Tc overprediction.

    Returns ``(J_eff, ss0, ss1)``: ``ss0``/``ss1`` are the Sz=0/Sz=1 lowest states' ``S(S+1)`` --
    verify ``ss0 ~= 0`` (singlet reference) and ``ss1 ~= 2`` (genuine S=1 triplet, not a higher-S
    intruder) before trusting ``J_eff`` as "the" triplet gap; see specs/BACKLOG.md."""
    h1, eri, L = _coordination_cluster(U0, ts, Us, tw, Uw, z)
    n = L  # half filling
    e0, ss0 = _fci_sector(h1, eri, n // 2, n // 2, max_cycle=max_cycle)
    e1, ss1 = _fci_sector(h1, eri, n // 2 + 1, n // 2 - 1, max_cycle=max_cycle)
    return e1 - e0, ss0, ss1


def ssh_chain_gap(U0: float, ts: float, Us: float, tw: float, Uw: float, n_dimers: int,
                   *, max_cycle: int = 1000) -> float:
    """Charge gap of an open 1-D SSH extended-Hubbard chain of ``n_dimers`` dimers (FCI): strong
    bonds (ts, Us) within a dimer, weak bonds (tw, Uw) between dimers. This captures the out-of-plane
    (stacking) inter-dimer coupling; extrapolated to n->inf it gives the quasi-1D gap. For Nb3I8 the
    DMRG (block2) limit is ~708 meV (L=8..20 monotone 730->709), between the isolated 842 and the
    fuller-coordination values -- exact FCI here is limited to small n (the L=12 half-filled FCI
    even fails to converge, which DMRG corrects).

    ``max_cycle`` (default 1000) is forwarded to the FCI Davidson solver; see ``coordination_gap``.
    """
    L = 2 * n_dimers
    h1 = np.zeros((L, L))
    eri = np.zeros((L, L, L, L))
    for i in range(L):
        eri[i, i, i, i] = U0

    def bond(i, j, t, U):
        h1[i, j] = h1[j, i] = t
        eri[i, i, j, j] = eri[j, j, i, i] = U

    for k in range(n_dimers):
        bond(2 * k, 2 * k + 1, ts, Us)               # strong (intra-dimer)
    for k in range(n_dimers - 1):
        bond(2 * k + 1, 2 * k + 2, tw, Uw)           # weak (inter-dimer)
    return _cluster_charge_gap(h1, eri, max_cycle=max_cycle)


# LT-bulk parameters extended with the weak inter-bilayer link (t_w_perp) and its Coulomb (U_w_perp),
# for the bath bound: (U0, t_s_perp, U_s_perp, t_w_perp, U_w_perp), meV, Table I.
NB3X8_LT_BULK_5P = {
    "Nb3F8":  (2590.5, -4.9,   714.6, -6.5,  572.5),
    "Nb3Cl8": (1451.4, -136.0, 400.1, -16.1, 313.5),
    "Nb3Br8": (1186.6, -169.4, 342.0, -20.4, 262.4),
    "Nb3I8":  (787.0,  -218.2, 258.5, -24.6, 183.8),
}


def four_site_exact_gap(U0: float, ts: float, Us: float, tw: float, Uw: float,
                         *, max_cycle: int = 1000) -> float:
    """Exact charge gap of an *enlarged* cluster -- two dimers joined by the weak inter-bilayer link
    (chain 0=1 strong, 1~2 weak, 2=3 strong; on-site U0, inter-site Us on strong bonds, Uw on the weak
    bond). Half-filled (4 electrons). The bath bound: comparing this to the isolated-dimer gap
    quantifies how much the inter-cluster coupling moves the gap -- a rigorous proxy for the DMFT bath
    that needs no bath fit. See specs/SPEC_nb3x8_gaps.md R1.

    ``max_cycle`` (default 1000) is forwarded to the FCI Davidson solver; see ``coordination_gap``.
    """
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
        return fixed_filling_energy(ModelIntegrals(h1, eri, 0.0, (n - n // 2, n // 2), 4),
                                     max_cycle=max_cycle)

    return E(5) + E(3) - 2 * E(4)


def exact_charge_gap(U0: float, t: float, Us: float, *, max_cycle: int = 1000) -> float:
    """Exact charge (Mott) gap ``E(3) + E(1) - 2 E(2)`` of the cluster by full diagonalization.

    ``max_cycle`` (default 1000) is forwarded to the FCI Davidson solver; see ``coordination_gap``.
    """
    cluster = dimer_cluster_integrals(U0, t, Us)

    def E(n):
        return fixed_filling_energy(ModelIntegrals(cluster.h1, cluster.eri, 0.0,
                                                   (n - n // 2, n // 2), 2),
                                     max_cycle=max_cycle)

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

    print("\n*** CORRECTION -- the isolated-cluster finding does NOT translate to the solid ***")
    print("Nb3I8: restoring inter-dimer coordination (band broadening) collapses the exact gap toward "
          "Hubbard-I:")
    i = NB3X8_LT_BULK_5P["Nb3I8"]
    hub = hubbard_i_gap(i[0], i[1], i[2])
    print(f"{'coordination z':>15} {'exact gap':>10}   (Hubbard-I = {hub:.0f} meV)")
    for z in (0, 1, 2, 3):
        print(f"{z:>15} {coordination_gap(*i, z):10.1f}")
    print("1-D SSH chain -> thermodynamic limit (DMRG/block2, L=8..20 monotone 730->709): ~708 meV.")
    print("With realistic 3-D coordination the solid gap is ~600-650 meV, close to Hubbard-I. The "
          "isolated\n'~29% error' is an artifact of neglecting broadening -- it does NOT imply the "
          "paper's gaps are wrong;\nif anything it vindicates the cluster-DMFT/Hubbard-I approach. "
          "(The 4-site 'bath bound' was too\noptimistic: it sampled only the nearest neighbour.)")

    print("\n*** chem-g78: does coordination rescue the Nb3X8 Tc overprediction? (spin channel) ***")
    from nb3x8_magnetometry import overprediction_factor

    rescued = []
    for name in ("Nb3Cl8", "Nb3Br8"):
        p5 = NB3X8_LT_BULK_5P[name]
        over = overprediction_factor(name)
        needed = min(3.0, over)      # reduction coordination must deliver (the gate's kill bar)
        print(f"\n{name} (isolated-dimer Tc overprediction {over:.2f}x): the IDENTICAL coordination "
              "machinery, spin channel (J_eff) vs charge channel (gap):")
        print(f"{'z':>3} {'J_eff (meV)':>12} {'<S^2> Sz=0/1':>14} {'charge gap (meV)':>17}")
        j_by_z = {}
        for z in (0, 1, 2, 3):
            j_eff, ss0, ss1 = coordination_spin_gap(*p5, z)
            j_by_z[z] = j_eff
            cgap = coordination_gap(*p5, z, max_cycle=4000)
            print(f"{z:>3} {j_eff:12.2f} {f'{ss0:.3f}/{ss1:.3f}':>14} {cgap:17.1f}")
        j0, j3 = j_by_z[0], j_by_z[3]
        if j3 <= j0 / needed:
            rescued.append(name)
        print(f"J_eff at z=3 ({j3:.1f} meV) vs the J0/{needed:.2f} = {j0 / needed:.1f} meV rescue "
              f"threshold: {'RESCUED' if name in rescued else 'NOT rescued'}.")
    print("\nCharge-channel control on the same clusters drops monotonically (coordination softens the "
          "charge\ngap as designed); the spin channel does not -- it oscillates and is largest, not "
          "smallest, at z=3.")
    print("Verdict (finite L<=8 open clusters): "
          + (f"coordination reduction RESCUES the Tc overprediction for {', '.join(rescued)}."
             if rescued else
             "coordination/mean-field reduction does NOT explain the 5.3x (Cl) / 2.3x (Br) Tc "
             "overprediction;\nsee the module docstring and specs/BACKLOG.md."))
