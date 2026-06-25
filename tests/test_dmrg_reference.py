#!/usr/bin/env python3
"""
Validation for the classical reference path (dmrg_reference.py).

The exact-FCI path (always available) must match PySCF CASCI on the shared integral convention
that the DMRG path also relies on; the DMRG path must fail with a clear ImportError when block2 is
absent, and ``reference_energy(auto)`` must fall back to exact FCI.

Run:  pytest tests/test_dmrg_reference.py -v
"""
import pytest

from hybrid_quantum_solver.dmrg_reference import (
    dmrg_available,
    fci_energy,
    reference_energy,
)


def _lih_cas22():
    from pyscf import ao2mo, gto, mcscf, scf
    mol = gto.M(atom="Li 0 0 0; H 0 0 1.6", basis="sto3g", verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    cas = mcscf.CASCI(mf, 2, 2)
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), 2)
    return h1, eri, cas.nelecas, float(e_core), float(cas.e_tot)


def test_fci_energy_matches_casci():
    h1, eri, n_elec, e_core, casci = _lih_cas22()
    assert abs(fci_energy(h1, eri, n_elec, e_core) - casci) < 1e-6


def test_reference_energy_auto_returns_exact():
    h1, eri, n_elec, e_core, casci = _lih_cas22()
    energy, method = reference_energy(h1, eri, n_elec, e_core, method="auto")
    assert abs(energy - casci) < 1e-6
    assert method == ("dmrg" if dmrg_available() else "fci")


def test_dmrg_energy_requires_block2():
    if dmrg_available():
        pytest.skip("block2 is installed; the missing-dependency path is not exercised")
    from hybrid_quantum_solver.dmrg_reference import dmrg_energy
    h1, eri, n_elec, e_core, _ = _lih_cas22()
    with pytest.raises(ImportError):
        dmrg_energy(h1, eri, n_elec, e_core)
