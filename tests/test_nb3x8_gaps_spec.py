"""
Acceptance gates G1-G8 for specs/SPEC_nb3x8_gaps.md (exact Nb3X8 cluster gaps vs Hubbard-I).

Test-first origin: the Nb3X8 bilayer cluster (a generalized Hubbard dimer, from the arXiv:2501.10320
cRPA parameters) is exactly diagonalizable; we compute the exact charge gap and the Hubbard-I gap and
quantify where the Hubbard-I approximation the source paper uses breaks down. Reference: exact
diagonalization + the t->0 atomic limit (both -> U0).

The gates tell an honest, self-correcting story:
  * G1-G4: facts about the ISOLATED cluster -- Hubbard-I is near-exact for strong correlation but
    underestimates the isolated iodide gap by ~29% (G3); the single-ratio U0/|t| law fails (G4).
  * G5: the one-neighbour "bath bound" -- small (~5%) but MISLEADING (it under-samples coordination).
  * G6 (the correction, and the real definition of done): restoring inter-dimer coordination collapses
    the exact gap toward Hubbard-I, so the isolated ~29% is an artifact that does NOT translate to the
    solid. The thermodynamic-limit story (block2 DMRG ~708 meV -> ~600-650 with 3-D coordination) is
    in nb3x8_gaps.ssh_chain_gap / the module docstring.
  * G7 (chem-q9g): ``coordination_gap``'s default ``max_cycle=1000`` FCI cap is too tight for the
    z=3 Cl/Br clusters (Davidson genuinely needs ~4x more iterations, not a near-degeneracy -- it
    agrees with a much looser cap to high precision); the function must accept a ``max_cycle``
    override rather than forcing callers to reimplement the cluster construction.
  * G8 (chem-g78): the SAME coordination clusters in the SPIN channel (``coordination_spin_gap``,
    J_eff = E(Sz=1) - E(Sz=0)) do NOT show the charge channel's softening, so coordination/mean-field
    reduction cannot be what explains the 5.3x Tc overprediction of the isolated-dimer J.

PySCF FCI / NumPy / SciPy only (no block2); `make gates` runs it in its own process.
"""
from scipy.stats import spearmanr

from nb3x8_gaps import (
    NB3X8_CLUSTERS,
    NB3X8_LT_BULK,
    NB3X8_LT_BULK_5P,
    coordination_gap,
    coordination_spin_gap,
    exact_charge_gap,
    four_site_exact_gap,
    hubbard_i_gap,
)
from odmd_spin import dimer_exchange_analytic


def test_G1_atomic_limit_validation():
    """t -> 0: both the exact and Hubbard-I gaps reduce to the atomic Mott gap U0, every cluster."""
    for p in NB3X8_CLUSTERS.values():
        assert abs(exact_charge_gap(p["U0"], 1e-6, p["Us"]) - p["U0"]) < 1e-3
        assert abs(hubbard_i_gap(p["U0"], 1e-6, p["Us"]) - p["U0"]) < 1e-3


def test_G2_exact_gaps_are_the_new_numbers():
    """The exact charge gaps (meV) the paper never reported; positive/insulating for every cluster."""
    expected = {
        "I  LT-bulk": 842.4, "Br LT-bulk": 1086.0, "Cl LT-bulk": 1311.8, "F  LT-bulk": 2580.8,
        "I  LT-bil": 1960.6, "Br LT-bil": 2281.0, "Cl LT-bil": 2550.0, "F  LT-bil": 3978.9,
        "Cl HT-bulk": 1369.0, "Br HT-bulk": 1091.9,
    }
    for name, gap in expected.items():
        computed = exact_charge_gap(**NB3X8_CLUSTERS[name])
        assert abs(computed - gap) < 1.0, (name, computed, gap)
        assert computed > 0.0


def test_G3_isolated_cluster_iodides_worst():
    """ISOLATED-CLUSTER FACT (superseded for the material by G6): strongly-correlated clusters are
    near-exact (<2%); the iodides are the worst, Hubbard-I underestimating the isolated gap. True for
    the isolated cluster, but G6 shows this does not translate to the solid."""
    def rel_err(name):
        p = NB3X8_CLUSTERS[name]
        return (hubbard_i_gap(**p) - exact_charge_gap(**p)) / exact_charge_gap(**p)

    # strongly-correlated clusters: Hubbard-I near-exact
    for name in ("F  LT-bulk", "F  LT-bil", "Cl HT-bulk", "Br HT-bulk"):
        assert abs(rel_err(name)) < 0.02, (name, rel_err(name))
    # iodides: worst, and underestimated (negative error)
    assert rel_err("I  LT-bulk") < -0.20, rel_err("I  LT-bulk")     # ~ -29%
    assert rel_err("I  LT-bil") < -0.10, rel_err("I  LT-bil")       # ~ -12%
    # the iodide-bulk error is the largest in magnitude over the whole set
    worst = max(NB3X8_CLUSTERS, key=lambda n: abs(rel_err(n)))
    assert worst == "I  LT-bulk", worst


def test_G4_single_parameter_law_is_falsified():
    """THE RECORDED NEGATIVE RESULT: the clean 4-point 'error ~ U0/|t|' law does NOT survive the full
    dataset. The correlation is strong but imperfect (Spearman between -0.7 and -1, not -1), and the
    error is not monotone in U0/|t| -- it depends on t and U_s_perp together, not one ratio."""
    ratios, errs = [], []
    for p in NB3X8_CLUSTERS.values():
        ratios.append(p["U0"] / abs(p["t"]))
        errs.append(abs(hubbard_i_gap(**p) - exact_charge_gap(**p)) / exact_charge_gap(**p))
    rho = spearmanr(ratios, errs).correlation
    assert -1.0 < rho < -0.7, rho                                   # strong but NOT perfect
    order = [e for _, e in sorted(zip(ratios, errs))]              # by increasing U0/|t|
    assert not all(order[i] >= order[i + 1] for i in range(len(order) - 1)), order   # non-monotone


def test_G5_nearest_neighbour_shift_is_small_but_misleading():
    """The nearest-neighbour cluster enlargement (2->4 sites) moves the Nb3I8 gap by only ~5% and is
    largest for Nb3F8 (ill-defined dimer, t_s ~ t_w). These facts hold -- but this ONE-neighbour
    'bath bound' is misleading: it under-samples the coordination (see G6). Kept as a recorded fact
    and a caution, not as evidence the isolated-cluster finding survives."""
    shift = {}
    for name, five in NB3X8_LT_BULK_5P.items():
        g2 = exact_charge_gap(five[0], five[1], five[2])
        shift[name] = abs(four_site_exact_gap(*five) - g2) / g2
    assert shift["Nb3I8"] < 0.06, shift                            # ~5% for one neighbour...
    assert shift["Nb3F8"] > shift["Nb3Cl8"] > shift["Nb3Br8"] > shift["Nb3I8"], shift


def test_G6_coordination_collapses_the_isolated_cluster_error():
    """THE CORRECTION (definition of done, superseding the isolated-cluster headline): restoring
    inter-dimer coordination -- the band broadening the isolated dimer omits and cluster-DMFT keeps --
    collapses the exact Nb3I8 gap monotonically toward the Hubbard-I value. So the isolated ~29%
    'error' is an artifact of the isolated-cluster approximation and does NOT translate to the solid;
    Hubbard-I (embedded in the full dispersion) is close to the true gap."""
    p = NB3X8_LT_BULK_5P["Nb3I8"]
    hub = hubbard_i_gap(p[0], p[1], p[2])
    gaps = [coordination_gap(*p, z) for z in (0, 1, 2, 3)]
    assert all(gaps[k + 1] < gaps[k] for k in range(3)), gaps        # monotone drop with coordination
    disc0 = gaps[0] - hub                                            # isolated exact - Hubbard-I
    disc3 = gaps[3] - hub                                            # after z=3 coordination
    assert disc3 < 0.35 * disc0, (gaps, hub)                         # >65% of the discrepancy closed
    assert gaps[3] < gaps[0], gaps                                   # solid gap far below isolated


def test_G7_coordination_gap_max_cycle_override():
    """chem-q9g: the default max_cycle=1000 is too tight for the z=3 Cl/Br clusters (Davidson does
    not converge in the N=9, (5,4)-electron sector), but the function accepts an override so callers
    aren't stuck reimplementing the cluster. A looser cap converges, and agrees with an even looser
    one to high precision -- confirming this is an iteration-cap issue, not genuine non-convergence
    or near-degeneracy."""
    import pytest

    for name in ("Nb3Cl8", "Nb3Br8"):
        p = NB3X8_LT_BULK_5P[name]
        with pytest.raises(RuntimeError, match="did not converge"):
            coordination_gap(*p, z=3)                                # default cap: fails
        g_loose = coordination_gap(*p, z=3, max_cycle=4000)
        g_looser = coordination_gap(*p, z=3, max_cycle=8000)
        assert abs(g_loose - g_looser) < 1e-6, (name, g_loose, g_looser)   # same answer, not drifting

    # z=0..2 (and Nb3I8/Nb3F8 at z=3) already converge under the default cap -- the override must not
    # change their answer.
    p = NB3X8_LT_BULK_5P["Nb3I8"]
    assert abs(coordination_gap(*p, z=3) - coordination_gap(*p, z=3, max_cycle=4000)) < 1e-9


def test_G8_spin_channel_coordination_does_not_rescue_Tc_overprediction():
    """chem-g78: SPEC_nb3x8_magnetometry.md section 7 names coordination/mean-field reduction as the
    likely explanation for the 5.3x (Nb3Cl8) Tc overprediction. Run the IDENTICAL coordination
    machinery in the SPIN channel (J_eff = E(Sz=1 lowest) - E(Sz=0 lowest) at half filling) instead
    of the charge channel, on Nb3Cl8, at L = 4, 6, 8 (z = 1, 2, 3; capped per the module's L=12
    half-filled FCI convergence warning).

    Machinery anchor: at z=0 (the isolated dimer) coordination_spin_gap must reduce to the exact
    closed-form dimer exchange (SPEC_odmd_spin) -- if this drifts, the FCI Sz-sector construction
    itself is wrong and nothing else here can be trusted.
    """
    cl = NB3X8_LT_BULK_5P["Nb3Cl8"]
    j0_analytic = dimer_exchange_analytic(**NB3X8_LT_BULK["Nb3Cl8"])
    j0, ss0_0, ss1_0 = coordination_spin_gap(*cl, 0)
    assert abs(j0 - j0_analytic) < 1e-6, (j0, j0_analytic)           # machinery anchor

    results = {z: coordination_spin_gap(*cl, z) for z in (0, 1, 2, 3)}
    for z, (j_eff, ss0, ss1) in results.items():
        # <S^2> check (the bead's caveat): Sz=0 lowest is a genuine singlet, Sz=1 lowest a genuine
        # S=1 triplet, not a higher-S intruder -- required before J_eff means "the" triplet gap.
        assert abs(ss0) < 1e-3, (z, ss0)
        assert abs(ss1 - 2.0) < 1e-2, (z, ss1)

    j_values = [results[z][0] for z in (0, 1, 2, 3)]

    # PRIMARY KILL: dies (restoring the published attribution) if coordination eventually delivers
    # J_eff(largest cluster) <= J0/3.
    assert j_values[3] > j0 / 3.0, (j_values, j0)

    # No cluster in L<=8 shows anything close to the needed >=3x reduction (open-boundary end
    # effects at L=4 are a bump, not the start of a trend -- gate the whole range, not just z=3).
    assert min(j_values) > j0 / 3.0, j_values

    # The spin gap does NOT fall monotonically the way the charge gap does (G6): it rises then
    # partially returns, and is *larger* at the largest cluster reached than in the isolated dimer.
    assert j_values[3] >= j0, j_values

    # SECONDARY KILL (control): the identical machinery's charge-channel gap on the SAME Nb3Cl8
    # clusters must still reproduce a real (>=20%) reduction, or the machinery itself is broken and
    # neither the spin nor the charge result can be trusted.
    charge0 = coordination_gap(*cl, 0)
    charge3 = coordination_gap(*cl, 3, max_cycle=2000)
    assert (charge0 - charge3) / charge0 >= 0.20, (charge0, charge3)

    # VERDICT (recorded either way): coordination/mean-field reduction does not explain the Tc
    # overprediction -- both kill criteria failed to trigger the published attribution, and the
    # control confirms the machinery works. See the module docstring and specs/BACKLOG.md.
