"""Regression gate for the ethylene torsion geometry (bead chem-cve).

The generator used to compute the torsion angle in radians and then throw it
away, emitting one identical planar geometry for all seven "torsion" points --
so the 90-degree diradical the sweep exists to probe was never computed, and
the family's seven results were seven copies of one number.

Not named test_*_spec.py on purpose: that glob is the physics spec gates the
Makefile runs one-per-process for block2 isolation. This needs no isolation.
"""
import math

import numpy as np
import pytest

from run_grand_sweep import (
    C2H4_HCC,
    C2H4_RCC,
    C2H4_RCH,
    ETHYLENE_TORSION_ANGLES,
    ethylene_geometry,
)

HARTREE2KCAL = 627.5094740631


def parse(geom: str):
    atoms = []
    for chunk in geom.split(";"):
        sym, *xyz = chunk.split()
        atoms.append((sym, np.array([float(v) for v in xyz])))
    return atoms


def angle_deg(v, w):
    c = np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w))
    return math.degrees(math.acos(np.clip(c, -1.0, 1.0)))


def dihedral_deg(p0, p1, p2, p3):
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m = b2 / np.linalg.norm(b2)
    return abs(math.degrees(math.atan2(np.dot(np.cross(n1, n2), m), np.dot(n1, n2))))


def test_every_torsion_point_is_a_distinct_geometry():
    """The exact defect: 7 angles must not collapse to 1 molecule string."""
    geoms = [ethylene_geometry(d) for d in ETHYLENE_TORSION_ANGLES]
    assert len(set(geoms)) == len(ETHYLENE_TORSION_ANGLES)


@pytest.mark.parametrize("deg", ETHYLENE_TORSION_ANGLES)
def test_torsion_rotates_only_the_dihedral(deg):
    """Twisting must preserve every bond length and valence angle."""
    atoms = parse(ethylene_geometry(deg))
    assert [a[0] for a in atoms] == ["C", "C", "H", "H", "H", "H"]
    c1, c2 = atoms[0][1], atoms[1][1]
    hydrogens = [a[1] for a in atoms[2:]]

    assert np.linalg.norm(c2 - c1) == pytest.approx(C2H4_RCC, abs=1e-4)
    for i, h in enumerate(hydrogens):
        own, far = (c1, c2) if i < 2 else (c2, c1)
        assert np.linalg.norm(h - own) == pytest.approx(C2H4_RCH, abs=1e-3)
        assert angle_deg(h - own, far - own) == pytest.approx(C2H4_HCC, abs=1e-2)


@pytest.mark.parametrize("deg", ETHYLENE_TORSION_ANGLES)
def test_dihedral_equals_the_requested_angle(deg):
    atoms = parse(ethylene_geometry(deg))
    c1, c2 = atoms[0][1], atoms[1][1]
    h_on_c1, h_on_c2 = atoms[2][1], atoms[4][1]
    assert dihedral_deg(h_on_c1, c1, c2, h_on_c2) == pytest.approx(deg, abs=1e-2)


def test_planar_at_zero_and_perpendicular_at_ninety():
    flat = parse(ethylene_geometry(0))
    assert max(abs(a[1][2]) for a in flat) == pytest.approx(0.0, abs=1e-6)
    twisted = parse(ethylene_geometry(90))
    # At 90 deg the C2 methylene leaves the xy-plane entirely.
    assert abs(twisted[4][1][1]) == pytest.approx(0.0, abs=1e-6)
    assert abs(twisted[4][1][2]) > 0.9


def test_torsion_barrier_is_physical():
    """The falsifiable part: breaking the pi bond must cost energy.

    CASCI(2,2)/sto-3g on this rigid scan gives ~86 kcal/mol against an
    experimental barrier near 65. The overestimate is expected and honest --
    the scan does not relax the twisted structure (real ethylene pyramidalizes
    and lengthens C-C at the transition state) and the basis is minimal. The
    gate is the sign, the monotonicity and the order of magnitude, not the
    number: those are what the identical-geometry bug destroyed.
    """
    pytest.importorskip("pyscf")
    from pyscf import gto, mcscf, scf

    energies = {}
    for deg in ETHYLENE_TORSION_ANGLES:
        mol = gto.M(atom=ethylene_geometry(deg), basis="sto-3g", spin=0, verbose=0)
        mc = mcscf.CASCI(scf.RHF(mol).run(), 2, 2)
        mc.kernel()
        energies[deg] = mc.e_tot

    assert len(set(energies.values())) == len(ETHYLENE_TORSION_ANGLES)

    ordered = [energies[d] for d in sorted(energies)]
    assert all(a < b for a, b in zip(ordered, ordered[1:])), energies

    barrier = (energies[90] - energies[0]) * HARTREE2KCAL
    assert 50.0 < barrier < 110.0, f"barrier {barrier:.1f} kcal/mol is unphysical"


def test_recovery_path_shares_one_generator():
    """compile_from_gcs must not re-derive the geometry; it was wrong twice."""
    pytest.importorskip("google.cloud.storage")
    from compile_from_gcs import get_standard_geometries

    cfg = get_standard_geometries()["ethylene_torsion"]
    assert cfg["geom_fn"] is ethylene_geometry
    assert list(cfg["coords"]) == list(ETHYLENE_TORSION_ANGLES)
