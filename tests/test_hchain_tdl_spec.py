"""
Acceptance gates G1-G5 for specs/SPEC_hchain_tdl.md (Hn DMRG bond-dim + thermodynamic-limit study).

block2 (DMRG) is required; the file is skipped without it. It MUST run in its own process
(block2 segfaults if it loads after pyscf/aer in the same interpreter -- see run_in_chem.sh).
"""
import pytest
from pyscf import gto, scf, ao2mo

from hybrid_quantum_solver.dmrg_reference import (
    dmrg_energy_extrapolated,
    thermodynamic_limit_fit,
    fci_energy,
    dmrg_available,
    ExtrapResult,
)

pytestmark = pytest.mark.skipif(not dmrg_available(), reason="block2 not installed")

R_ANG = 1.8 * 0.529177210903   # 1.8 Bohr in Angstrom (~0.9525)


def integrals(n):
    """Full-valence (n,n) minimal-basis integrals for an H_n chain."""
    atom = "; ".join(f"H 0 0 {i * R_ANG:.6f}" for i in range(n))
    mol = gto.M(atom=atom, basis="sto-6g", verbose=0)
    mf = scf.RHF(mol).run()
    mo = mf.mo_coeff
    h1 = mo.T @ mf.get_hcore() @ mo
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mo), n)
    ne = (mol.nelectron // 2, mol.nelectron // 2)
    return h1, eri, ne, float(mol.energy_nuc())


# Well-converged dims: the extrapolation is SOUND (lands near FCI). At these small, FCI-tractable
# sizes the largest-D energy is itself near-exact, so extrapolation validates correctness, not
# improvement -- the improvement only materialises in the under-converged large-n regime (where no
# FCI exists to compare). See specs/SPEC_hchain_tdl.md G1 note.
CONV_DIMS = (80, 160, 300)
# Coarse dims: deliberately under-converged so truncation error is large -- used for the
# monotone-convergence gate (G3).
COARSE_DIMS = (20, 30, 40, 60)
G1_TOL = 2e-4   # extrapolation accuracy at CONV_DIMS for n up to 12 (tighter with larger D)


def test_G1_extrapolation_matches_fci():
    for n in (10, 12):
        h1, eri, ne, ec = integrals(n)
        e_fci = fci_energy(h1, eri, ne, ec)
        res = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=CONV_DIMS)
        assert isinstance(res, ExtrapResult)
        assert res.method == "dweight"
        assert abs(res.energy - e_fci) < G1_TOL, (n, res.energy, e_fci)


def test_G2_variational_window():
    h1, eri, ne, ec = integrals(10)
    e_fci = fci_energy(h1, eri, ne, ec)
    res = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=CONV_DIMS)
    e_dmax = [e for _, _, e in res.per_D][-1]
    assert res.energy <= e_dmax + 1e-9            # extrapolation lowers the largest-D raw energy
    assert res.energy >= e_fci - G1_TOL           # not meaningfully below FCI (mild fit overshoot)


def test_G3_monotone_truncation_convergence():
    h1, eri, ne, ec = integrals(12)
    res = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=COARSE_DIMS)
    Es = [e for _, _, e in res.per_D]
    dws = [dw for _, dw, _ in res.per_D]
    assert all(Es[i + 1] <= Es[i] + 1e-9 for i in range(len(Es) - 1)), Es
    assert all(dws[i + 1] <= dws[i] + 1e-10 for i in range(len(dws) - 1)), dws


def test_G4_thermodynamic_limit_stability():
    ns = [6, 8, 10, 12]
    e_per_atom = []
    for n in ns:
        h1, eri, ne, ec = integrals(n)
        res = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=(40, 60, 100))
        e_per_atom.append(res.energy / n)
    e_inf_all, _ = thermodynamic_limit_fit(ns, e_per_atom)
    e_inf_drop, _ = thermodynamic_limit_fit(ns[:-1], e_per_atom[:-1])
    assert abs(e_inf_all - e_inf_drop) < 1e-3, (e_inf_all, e_inf_drop)


def test_G5_reproducible():
    h1, eri, ne, ec = integrals(10)
    r1 = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=(30, 50), seed=42)
    r2 = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=(30, 50), seed=42)
    assert abs(r1.energy - r2.energy) < 1e-8
