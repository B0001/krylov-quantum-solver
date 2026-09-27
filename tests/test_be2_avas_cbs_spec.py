"""
Acceptance gates G1-G4 for specs/SPEC_be2_avas_cbs.md.

Claim: replacing the CAS(4,8)-on-canonical-HF-virtuals active space (SPEC_be2_cbs.md, found by
bead chem-mom's G5 to drift with basis: outermost active orbital rms spread 7.1 -> 8.3 -> 9.4
bohr TZ->QZ->5Z) with an AVAS(Be 2s/2p, minao='ano')-selected CAS(4,8) removes that confound, and
retests chem-mom's pre-registered CBS(QZ/5Z) vs CBS(TZ/QZ) consistency bar (~20 cm^-1) on the
now-basis-consistent active space.

Finding (G3/G4): G1/G2 confirm the active-space confound is gone (same ncas/nelecas everywhere,
per-orbital spread agrees to <0.01 bohr TZ vs 5Z, not 2.3 bohr). That alone does NOT rescue the
CBS attribution: the QZ/5Z shift shrinks from chem-mom's +1224 cm^-1 to +530 cm^-1 but is still
>> the 20 cm^-1 bar (G3), and the single-basis wells are still non-monotone in X -- though the
decomposition shows cc-pVTZ, not a basis-inconsistent active space, is the outlier (G4).

PySCF only (no block2); cheap R-grid ({2.4, 2.45, 2.6} + {8.0}) keeps this to ~1 min (5Z point
~10 s dominates). `make gates` runs this spec in its own process.
"""
import numpy as np
import pytest

from be2_avas_cbs import (
    active_orbital_diagnostics,
    avas_active_space,
    avas_casci_nevpt2_point,
)
from be2_cbs import BASIS_CARDINAL, HA2CM, cbs_extrapolate_correlation, quadratic_well
from pyscf import gto, scf

RS = (2.4, 2.45, 2.6, 8.0)
BASES = ("ccpvtz", "ccpvqz", "ccpv5z")
CBS_CONFIRM_THRESHOLD_CM1 = 20.0  # pre-registered by bead chem-mom (SPEC_be2_cbs.md G5)


@pytest.fixture(scope="module")
def points():
    """Every CASCI+NEVPT2 point the gates need, computed once and shared."""
    return {(basis, R): avas_casci_nevpt2_point(R, basis) for basis in BASES for R in RS}


def _cbs_curve(points, Rs, lo_basis, hi_basis):
    out = {}
    for R in Rs:
        lo, hi = points[(lo_basis, R)], points[(hi_basis, R)]
        e_corr_cbs = cbs_extrapolate_correlation(BASIS_CARDINAL[lo_basis], lo.e_corr,
                                                 BASIS_CARDINAL[hi_basis], hi.e_corr)
        out[R] = hi.e_casci + e_corr_cbs
    return out


def _well(points, basis, Rs=RS):
    curve = {R: points[(basis, R)].e_tot for R in Rs}
    return quadratic_well([Rs[0], Rs[1], Rs[2]], [curve[Rs[0]], curve[Rs[1]], curve[Rs[2]]],
                          Rs[3], curve[Rs[3]])


def _cbs_well(points, lo, hi, Rs=RS):
    curve = _cbs_curve(points, Rs, lo, hi)
    return quadratic_well([Rs[0], Rs[1], Rs[2]], [curve[Rs[0]], curve[Rs[1]], curve[Rs[2]]],
                          Rs[3], curve[Rs[3]])


# --- G1: the active space has the expected shape everywhere it's used -------------------------


def test_G1_avas_active_space_is_4_electrons_in_8_orbitals_everywhere():
    """avas_active_space must return ncas=8, nelecas=4 at every (basis, R) checked -- a mismatch
    is the literal failure mode chem-mom diagnosed (a silently different-shaped active space) and
    is a hard RuntimeError from the production code itself (checked here, not just asserted)."""
    for basis in BASES:
        for R in RS:
            mol = gto.M(atom=f"Be 0 0 0; Be 0 0 {R}", basis=basis, spin=0, verbose=0)
            mf = scf.RHF(mol).run()
            ncas, nelecas, mo = avas_active_space(mf)  # raises RuntimeError on mismatch
            assert ncas == 8
            assert nelecas == 4
            assert mo.shape[1] == mf.mo_coeff.shape[1]


def test_G1_minao_minao_default_would_have_dropped_2p(caplog):
    """Documents WHY minao='ano' is required: pyscf's AVAS default (minao='minao') gives Be no
    occupied-or-virtual 2p shell to project onto (Be's minimal ground-state basis is 1s^2 2s^2),
    so it silently selects only the 2s orbitals (ncas=2), not the intended CAS(4,8)."""
    from pyscf.mcscf import avas

    mol = gto.M(atom="Be 0 0 0; Be 0 0 2.5", basis="ccpvtz", spin=0, verbose=0)
    mf = scf.RHF(mol).run()
    ncas, nelecas, _ = avas.avas(mf, ("Be 2s", "Be 2p"), minao="minao", canonicalize=False)
    # ncas=2 (only the two Be 2s AOs are found -- no 2p shell in Be's minimal ground-state basis),
    # nelecas=4 (both fully doubly occupied: 0 active orbitals come from the virtual space) --
    # a degenerate, correlation-free "active space", not the intended CAS(4,8).
    assert (ncas, nelecas) == (2, 4), "if this changes, minao='minao' may now be usable for Be"


# --- G2: the active-space character is stable across basis (the actual fix) --------------------


def test_G2_active_orbital_spread_agrees_tz_vs_5z():
    """Per-orbital rms spread (sorted -- degenerate p-pairs can swap order) must agree between
    cc-pVTZ and cc-pV5Z to within 1.0 bohr at R=2.5 (near Re) and R=8.0 (asymptote) -- vs. the
    canonical-orbital baseline's 2.3 bohr (7.1->9.4) drift on the outermost orbital
    (results/be2_cbs_5z/diagnose.txt). This is the gate that would fail if AVAS were not actually
    basis-consistent."""
    for R in (2.5, 8.0):
        _, spread_tz = active_orbital_diagnostics("ccpvtz", R)
        _, spread_5z = active_orbital_diagnostics("ccpv5z", R)
        diff = np.abs(np.sort(spread_tz) - np.sort(spread_5z))
        assert diff.max() < 1.0, (R, spread_tz, spread_5z)


# --- G3: CBS attribution, retested on the fixed active space (chem-mom's own bar) ---------------


def test_G3_cbs_attribution_still_not_confirmed_but_residual_shrinks(points):
    """Pre-registered bar (bead chem-mom, SPEC_be2_cbs.md G5): |De_CBS(QZ/5Z) - De_CBS(TZ/QZ)|
    <= 20 cm^-1 would CONFIRM the method attribution. MEASURED here: shift ~+530 cm^-1 -- still
    NOT CONFIRMED, but well under chem-mom's canonical-orbital +1224 cm^-1 (fixing the active-space
    confound helped, did not solve it). Pins both wells and the shift as a regression."""
    _, de_tzqz = _cbs_well(points, "ccpvtz", "ccpvqz")
    _, de_qz5z = _cbs_well(points, "ccpvqz", "ccpv5z")
    shift = de_qz5z - de_tzqz

    assert abs(shift) > CBS_CONFIRM_THRESHOLD_CM1  # the pre-registered confirmation test FAILS
    assert 0.0 < de_tzqz < 200.0, de_tzqz
    assert 450.0 < de_qz5z < 750.0, de_qz5z
    assert 400.0 < shift < 700.0, shift

    # regression bound on the improvement itself: this composition's shift must stay well below
    # chem-mom's canonical-orbital +1224 cm^-1, or the "fix" claim needs revisiting
    assert shift < 0.7 * 1224.1


# --- G4: single-basis De monotone in X, or the failure mode identified -------------------------


def test_G4_single_basis_wells_nonmonotone_but_tz_is_the_identified_outlier(points):
    """MEASURED: TZ/QZ/5Z single-basis wells (~731/305/446 cm^-1) are still non-monotone in X.
    Decomposition (De = De_casci_only + De_from_nevpt2) shows QZ and 5Z agree with each other far
    better than either agrees with TZ on BOTH pieces -- consistent with cc-pVTZ being pre-
    asymptotic for this CAS(4,8)+NEVPT2 recipe, not with the active space still being basis-
    inconsistent (which G1/G2 directly rule out). Pins the wells and the ordering claim."""
    de = {b: _well(points, b)[1] for b in BASES}
    assert 600.0 < de["ccpvtz"] < 850.0, de
    assert 200.0 < de["ccpvqz"] < 400.0, de
    assert 350.0 < de["ccpv5z"] < 550.0, de
    assert not (de["ccpvtz"] < de["ccpvqz"] < de["ccpv5z"])  # non-monotone, on the record

    win = [2.4, 2.45, 2.6]

    def de_casci_only(basis):
        Re, _ = _well(points, basis)
        pc = np.poly1d(np.polyfit(win, [points[(basis, R)].e_casci for R in win], 2))
        return (points[(basis, 8.0)].e_casci - pc(Re)) * HA2CM

    de_cas = {b: de_casci_only(b) for b in BASES}
    de_corr = {b: de[b] - de_cas[b] for b in BASES}

    # QZ/5Z agree with each other much better than either does with TZ, on both pieces
    assert abs(de_cas["ccpvqz"] - de_cas["ccpv5z"]) < abs(de_cas["ccpvtz"] - de_cas["ccpvqz"])
    assert abs(de_corr["ccpvqz"] - de_corr["ccpv5z"]) < abs(de_corr["ccpvtz"] - de_corr["ccpvqz"])
