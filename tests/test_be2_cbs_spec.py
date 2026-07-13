"""
Acceptance gates G1-G4 for specs/SPEC_be2_cbs.md.

Claim (BACKLOG.md): core-valence correlation + a cc-pVXZ->CBS extrapolation moves the Be2 well
depth from the frozen-core CAS(4,8)/cc-pVDZ FCI baseline (~305 cm^-1) toward the experimental
929.7 cm^-1 (Merritt, Bondybey & Heaven, Science 2009).

Finding (G4): the composition moves toward experiment (G3) but does NOT reach the original
100 cm^-1 / 0.1 A tolerance -- CASSCF orbital optimization is numerically unstable for Be2 (tried
first, rejected: unconverged / non-monotonic across basis, see specs/SPEC_be2_cbs.md R1), so this
uses the validated fixed-orbital CASCI reference instead, which leaves no orbital relaxation to
help close the gap. G4 pins the actual (Re, De), explicitly outside the original bounds.

PySCF only (no block2); a small R-grid keeps this to ~1-2 minutes (QZ dominates, ~13 s/point).
`make gates` runs this spec in its own process.
"""
import pytest

from be2_cbs import casci_nevpt2_point, cbs_extrapolate_correlation, quadratic_well, HA2CM

EXPERIMENT_DE = 929.7   # cm^-1, Merritt, Bondybey & Heaven, Science 323, 1671 (2009)
EXPERIMENT_RE = 2.4498  # Angstrom


@pytest.fixture(scope="module")
def points():
    """Every CASCI+NEVPT2 point the gates need, computed once and shared."""
    pts = {}
    for R in (2.45, 6.0, 8.0):
        pts[("ccpvdz", R)] = casci_nevpt2_point(R, "ccpvdz")
    for basis in ("ccpvtz", "ccpvqz"):
        for R in (2.4, 2.45, 2.6, 6.0, 8.0):
            pts[(basis, R)] = casci_nevpt2_point(R, basis)
    return pts


def _cbs_curve(points, Rs):
    """CASCI(QZ) reference + CBS(TZ=3/QZ=4)-extrapolated NEVPT2 correlation, per R."""
    out = {}
    for R in Rs:
        lo, hi = points[("ccpvtz", R)], points[("ccpvqz", R)]
        e_corr_cbs = cbs_extrapolate_correlation(3, lo.e_corr, 4, hi.e_corr)
        out[R] = hi.e_casci + e_corr_cbs
    return out


# --- G1: small-basis curve has no real well; TZ/QZ do; asymptote is usably flat ------------


def test_G1_dz_has_no_well_near_re_but_tz_qz_do(points):
    """The "~305 cm^-1" backlog baseline is a far-R artifact (R~4.5 A), not a near-Re well:
    cc-pVDZ is LESS bound at R=2.45 than at R=8.0. cc-pVTZ/QZ, in contrast, ARE bound there."""
    dz_245 = points[("ccpvdz", 2.45)].e_tot
    dz_80 = points[("ccpvdz", 8.0)].e_tot
    assert (dz_245 - dz_80) * HA2CM > 0, "DZ should be UNBOUND at the physical bond length"

    for basis in ("ccpvtz", "ccpvqz"):
        e_245 = points[(basis, 2.45)].e_tot
        e_80 = points[(basis, 8.0)].e_tot
        assert (e_245 - e_80) * HA2CM < 0, f"{basis} should be BOUND at R=2.45"


def test_G1_asymptote_is_flat_enough(points):
    """R=6.0 and R=8.0 agree to within an order of magnitude of the well depth at TZ/QZ --
    flat enough to define De from, even though not fully R->infinity converged."""
    for basis in ("ccpvtz", "ccpvqz"):
        e_60 = points[(basis, 6.0)].e_tot
        e_80 = points[(basis, 8.0)].e_tot
        assert abs(e_60 - e_80) * HA2CM < 150.0


# --- G2: CBS direction is sane ----------------------------------------------------------------


def test_G2_correlation_energy_grows_with_basis(points):
    """More basis functions recover more NEVPT2 dynamic correlation at the physical bond length
    -- the expected direction for the CBS extrapolation (DZ is excluded per G1's finding)."""
    e_corr_tz = points[("ccpvtz", 2.45)].e_corr
    e_corr_qz = points[("ccpvqz", 2.45)].e_corr
    assert abs(e_corr_tz) < abs(e_corr_qz)


# --- G3: the composition moves toward experiment (the literal backlog claim) ------------------


def test_G3_cbs_well_is_closer_to_experiment_than_the_baseline(points):
    """|De_cbs - 929.7| < |De_baseline - 929.7|: CBS(TZ/QZ)+CV genuinely moves the well depth
    toward experiment relative to the frozen-core CAS(4,8)/cc-pVDZ FCI baseline this backlog item
    names (reproduced here from study_be2.py's own grid via the validated fci_energy)."""
    from pyscf import ao2mo, gto, mcscf, scf

    from hybrid_quantum_solver.dmrg_reference import fci_energy

    bond_lengths = [2.0, 2.2, 2.45, 2.7, 3.0, 3.5, 4.5, 6.0]  # study_be2.py's own grid
    rows = []
    for R in bond_lengths:
        mol = gto.M(atom=f"Be 0 0 0; Be 0 0 {R}", basis="ccpvdz", spin=0, verbose=0)
        mf = scf.RHF(mol).run()
        cas = mcscf.CASCI(mf, 8, 4)
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), 8)
        nelecas = (int(cas.nelecas[0]), int(cas.nelecas[1]))
        rows.append((R, fci_energy(h1, eri, nelecas, e_core)))
    asy = rows[-1][1]
    _, e_min = min(rows, key=lambda t: t[1])
    de_baseline = (asy - e_min) * HA2CM
    assert abs(de_baseline - 304.8) < 5.0, de_baseline  # pin the reproduced backlog number

    cbs = _cbs_curve(points, [2.4, 2.45, 2.6, 8.0])
    _, de_cbs = quadratic_well([2.4, 2.45, 2.6], [cbs[2.4], cbs[2.45], cbs[2.6]], 8.0, cbs[8.0])

    assert abs(de_cbs - EXPERIMENT_DE) < abs(de_baseline - EXPERIMENT_DE)


# --- G4: the honest boundary -- the original backlog gate does NOT hold ------------------------


def test_G4_cbs_well_falls_short_of_the_original_gate_by_a_pinned_amount(points):
    """The original backlog tolerance (|De-930|<100, |Re-2.45|<0.1) is NOT met. This pins the
    actual measured well as a regression so a future change is caught whichever way it moves --
    silently regressing further, or silently landing inside the original tolerance without a
    recorded reason (which would mean this finding needs revisiting, not deleting)."""
    cbs = _cbs_curve(points, [2.4, 2.45, 2.6, 8.0])
    Re, De = quadratic_well([2.4, 2.45, 2.6], [cbs[2.4], cbs[2.45], cbs[2.6]], 8.0, cbs[8.0])

    assert 2.50 < Re < 2.70, Re
    assert 400.0 < De < 550.0, De

    # explicitly outside the original backlog tolerance -- the falsification, on the record
    assert not (abs(De - EXPERIMENT_DE) < 100.0 and abs(Re - EXPERIMENT_RE) < 0.1)
