"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_gaps.md (exact Nb3X8 cluster gaps vs Hubbard-I).

Test-first origin: the Nb3X8 bilayer cluster (a generalized Hubbard dimer, from the arXiv:2501.10320
cRPA parameters) is exactly diagonalizable; we compute the exact charge gap and the Hubbard-I gap and
quantify where the Hubbard-I approximation the source paper uses breaks down. Reference: exact
diagonalization + the t->0 atomic limit (both -> U0).

The gates encode the finding AS REVISED by the full 10-cluster dataset: the material-level conclusion
(iodides worst; strongly-correlated clusters near-exact) is robust, but the clean single-parameter
U0/|t| law seen in the 4-point LT-bulk subset does NOT survive (Spearman ~ -0.86, not -1) -- that
mismatch is itself a recorded finding.

PySCF FCI / NumPy / SciPy only (no block2); `make gates` runs it in its own process.
"""
from scipy.stats import spearmanr

from nb3x8_gaps import (
    NB3X8_CLUSTERS,
    NB3X8_LT_BULK_5P,
    exact_charge_gap,
    four_site_exact_gap,
    hubbard_i_gap,
)


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


def test_G3_material_level_finding_is_robust():
    """DEFINITION OF DONE (robust across all 10 clusters): strongly-correlated clusters are near-exact
    (<2%); the iodides are consistently the worst, with Hubbard-I underestimating their gap."""
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


def test_G5_bath_bound_finding_survives_cluster_enlargement():
    """R1 (the isolated-cluster caveat) is bounded: enlarging the correlated region to include the
    inter-cluster weak link moves the Nb3I8 gap by only ~5% -- far below the ~29% Hubbard-I error --
    so the finding travels toward the solid. The bath effect is largest for Nb3F8 (ill-defined dimer,
    t_s ~ t_w), but Hubbard-I is exact there anyway."""
    shift = {}
    for name, five in NB3X8_LT_BULK_5P.items():
        g2 = exact_charge_gap(five[0], five[1], five[2])
        shift[name] = abs(four_site_exact_gap(*five) - g2) / g2
    assert shift["Nb3I8"] < 0.06, shift                            # iodide dimer well-isolated (~5%)
    assert shift["Nb3I8"] < abs(  # bath shift far below the Hubbard-I error for the iodide
        (hubbard_i_gap(787.0, -218.2, 258.5) - exact_charge_gap(787.0, -218.2, 258.5))
        / exact_charge_gap(787.0, -218.2, 258.5)) / 4, shift
    # bath effect is largest where the dimer is ill-defined (F: t_s ~ t_w), smallest for the iodide
    assert shift["Nb3F8"] > shift["Nb3Cl8"] > shift["Nb3Br8"] > shift["Nb3I8"], shift
