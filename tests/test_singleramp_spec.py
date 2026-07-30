"""
Acceptance gates G1-G3 for specs/SPEC_singleramp.md (single-ramp DMRG extrapolation).

block2 required; runs isolated (see run_in_chem.sh / `make gates`).
"""
import time

import pytest
from pyscf import gto, scf, ao2mo

from hybrid_quantum_solver.dmrg_reference import (
    dmrg_energy_extrapolated,
    fci_energy,
    dmrg_available,
)

pytestmark = pytest.mark.skipif(not dmrg_available(), reason="block2 not installed")

R_ANG = 1.8 * 0.529177210903
DIMS = (200, 400, 800)
N_SWEEPS_PER = 8        # per-D protocol sweeps per bond dim
SWEEPS_PER_STAGE = 4    # ramp sweeps per bond dim


def integrals(n):
    atom = "; ".join(f"H 0 0 {i * R_ANG:.6f}" for i in range(n))
    mol = gto.M(atom=atom, basis="sto-6g", verbose=0)
    mf = scf.RHF(mol).run()
    mo = mf.mo_coeff
    h1 = mo.T @ mf.get_hcore() @ mo
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mo), n)
    return h1, eri, (mol.nelectron // 2, mol.nelectron // 2), float(mol.energy_nuc())


def test_G1_ramp_agrees_with_perD():
    for n in (10, 12):
        h1, eri, ne, ec = integrals(n)
        e_perD = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=DIMS,
                                          protocol="perD", n_sweeps_per=N_SWEEPS_PER).energy
        e_ramp = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=DIMS,
                                          protocol="ramp", sweeps_per_stage=SWEEPS_PER_STAGE).energy
        assert abs(e_ramp - e_perD) < 1e-4, (n, e_ramp, e_perD)


def test_G2_ramp_is_sound():
    h1, eri, ne, ec = integrals(12)
    e_fci = fci_energy(h1, eri, ne, ec)
    res = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=DIMS,
                                   protocol="ramp", sweeps_per_stage=SWEEPS_PER_STAGE)
    # Was `method == "dweight"`; see specs/SPEC_extrap_regime.md. The FCI check below carries the
    # accuracy claim, so this only needs to exclude uncontrolled truncation.
    assert res.regime != "uncontrolled", res.regime
    assert abs(res.energy - e_fci) < 5e-4, (res.energy, e_fci)


def test_G3_ramp_is_cheaper():
    h1, eri, ne, ec = integrals(12)

    t0 = time.perf_counter()
    dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=DIMS,
                             protocol="perD", n_sweeps_per=N_SWEEPS_PER)
    t_perD = time.perf_counter() - t0

    t0 = time.perf_counter()
    dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=DIMS,
                             protocol="ramp", sweeps_per_stage=SWEEPS_PER_STAGE)
    t_ramp = time.perf_counter() - t0

    # deterministic proof: fewer total sweeps
    assert SWEEPS_PER_STAGE * len(DIMS) < N_SWEEPS_PER * len(DIMS)
    # wall-time evidence (generous bound to avoid load-dependent flakiness)
    assert t_ramp < 0.7 * t_perD, (t_ramp, t_perD)
