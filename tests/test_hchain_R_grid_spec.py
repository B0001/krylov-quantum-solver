"""
Acceptance gates G1-G2 for specs/SPEC_hchain_R_grid.md (H10 EOS: published 10-point R grid).

Pure PySCF FCI -- no block2/DMRG, so unlike test_hchain_tdl_spec.py this needs no process
isolation and can run alongside the rest of the suite.
"""
import csv
import os

from benchmark_hchain_tdl import integrals as driver_integrals
from hybrid_quantum_solver.dmrg_reference import fci_energy

REFERENCE_CSV = os.path.join(os.path.dirname(__file__), "..", "specs", "hchain_R_grid_motta_fci.csv")
N = 10
G1_TOL = 1e-4   # bead chem-4y9 acceptance criterion, verbatim
G2_TOL = 1e-8   # unitary rotation of the full space: FCI must be invariant to orbital choice


def _load_reference():
    with open(REFERENCE_CSV) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 10, "expected all 10 Simons Collaboration grid points"
    return [(float(r["r_bohr"]), float(r["fci_total_ha_n10"])) for r in rows]


def test_G1_fci_matches_published_grid():
    for r_bohr, e_ref in _load_reference():
        h1, eri, ne, ec, _ = driver_integrals(N, R_bohr=r_bohr)
        e_fci = fci_energy(h1, eri, ne, ec)
        assert abs(e_fci - e_ref) < G1_TOL, (r_bohr, e_fci, e_ref)


def test_G2_localized_invariant_across_grid():
    for r_bohr, _ in _load_reference():
        e_can = fci_energy(*driver_integrals(N, R_bohr=r_bohr)[:4])
        e_loc = fci_energy(*driver_integrals(N, R_bohr=r_bohr, localize=True)[:4])
        assert abs(e_loc - e_can) < G2_TOL, (r_bohr, e_can, e_loc)
