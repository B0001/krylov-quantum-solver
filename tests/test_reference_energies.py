#!/usr/bin/env python3
"""
Phase 1 validation gate (see REFACTOR_PLAN.md).

The corrected qubit Hamiltonian (hybrid_quantum_solver.molecular_hamiltonian) must
reproduce Full Configuration Interaction (FCI) to <1e-6 Ha on small molecules where
the exact answer is known, cross-checked against an INDEPENDENT PySCF FCI solve.

It also documents the regression: the original AdvancedStochasticCompactor does NOT
reproduce FCI, which is the whole reason for this refactor.

Run:  pytest tests/test_reference_energies.py -v
"""
import numpy as np
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

# Chemical accuracy is ~1.6 mHa; we demand far tighter agreement (exact diagonalisation
# of the same integrals two different ways must agree to numerical precision).
TOL = 1e-6

MOLECULES = {
    #  name : (geometry, basis, expected FCI total energy / Ha)
    "H2":  ("H 0 0 0; H 0 0 0.74",                              "sto3g", -1.137284),
    "H4":  ("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0",         "sto3g", -2.166387),
    "LiH": ("Li 0 0 0; H 0 0 1.6",                             "sto3g", -7.882324),
}


def _pyscf_fci_total(atom, basis):
    """Independent ground-truth: PySCF RHF -> FCI total energy (incl. nuclear repulsion)."""
    from pyscf import gto, scf, fci
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    e_fci, _ = fci.FCI(mf).kernel()
    return float(e_fci)


def _pyscf_rhf_total(atom, basis):
    from pyscf import gto, scf
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol)
    return float(mf.kernel())


@pytest.mark.parametrize("name", list(MOLECULES))
def test_qubit_hamiltonian_reproduces_fci(name):
    """min eig(qubit H) + offset == FCI total, vs both a literal value and a live PySCF FCI."""
    atom, basis, expected = MOLECULES[name]
    mh = build_molecular_hamiltonian(atom=atom, basis=basis)

    e_qubit = mh.ground_state_energy()
    e_fci = _pyscf_fci_total(atom, basis)

    assert abs(e_qubit - e_fci) < TOL, (
        f"{name}: qubit ground state {e_qubit:.8f} != PySCF FCI {e_fci:.8f}")
    assert abs(e_qubit - expected) < 1e-5, (
        f"{name}: qubit ground state {e_qubit:.6f} != expected {expected:.6f}")


@pytest.mark.parametrize("name", list(MOLECULES))
def test_hf_reference_matches_rhf(name):
    """The Hartree-Fock reference state energy must equal the PySCF RHF total energy.

    This guards against the original bug of starting the subspace from the empty
    vacuum |00..0> (which has zero electrons and <0|H|0> = 0).
    """
    atom, basis, _ = MOLECULES[name]
    mh = build_molecular_hamiltonian(atom=atom, basis=basis)
    assert abs(mh.hf_energy - _pyscf_rhf_total(atom, basis)) < TOL


def test_active_space_matches_casci():
    """The active-space path must reproduce PySCF CASCI (needed for the eventual CAS(8,8))."""
    from pyscf import gto, scf, mcscf
    atom, basis = "Li 0 0 0; H 0 0 1.6", "sto3g"
    mh = build_molecular_hamiltonian(atom=atom, basis=basis,
                                     active_electrons=2, active_orbitals=2)
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    cas = mcscf.CASCI(mf, 2, 2)
    cas.kernel()
    assert abs(mh.ground_state_energy() - float(cas.e_tot)) < TOL


def test_integrals_bridge_matches_casci():
    """The CIF/CASCI materials bridge (raw active-space integrals) must match PySCF CASCI."""
    from pyscf import gto, scf, mcscf, ao2mo
    from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals

    mol = gto.M(atom="Li 0 0 0; H 0 0 1.6", basis="sto3g", verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    cas = mcscf.CASCI(mf, 2, 2)
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), 2)

    mh = build_hamiltonian_from_integrals(h1, eri, num_particles=(1, 1), energy_offset=float(e_core))
    assert abs(mh.ground_state_energy() - float(cas.e_tot)) < TOL


def test_original_compactor_does_not_reproduce_fci():
    """Regression marker: the OLD hand-rolled mapping is wrong by ~0.3 Ha on H2.

    If this ever starts passing (i.e. the old compactor suddenly matches FCI), the
    test is stale and should be revisited -- but as shipped, the original pipeline
    cannot reproduce the reference energy.
    """
    from pyscf import gto, scf, mcscf, ao2mo
    from qiskit.quantum_info import SparsePauliOp
    from hybrid_quantum_solver.orchestrate_hybrid_pipeline import AdvancedStochasticCompactor
    from hybrid_quantum_solver.chemistry_gateway import translate_tensors_to_orchestrator

    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto3g", verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    cas = mcscf.CASCI(mf, 2, 2)
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), 2)

    # Build the OLD mapping directly (no orchestrator / legacy sampler needed).
    single_body, two_body = translate_tensors_to_orchestrator(h1, eri)
    compactor = AdvancedStochasticCompactor(n_spin_orbitals=2 * h1.shape[0])
    for p, q, w in single_body:
        compactor.map_one_body_term(p, q, w)
    for p, q, r, s, w in two_body:
        compactor.map_two_body_term(p, q, r, s, w)
    compactor.finalize_and_compile_metrics()

    op = SparsePauliOp([p[::-1] for p in compactor.pauli_strings],
                       np.asarray(compactor.coefficients, dtype=float))
    e_old_total = float(np.linalg.eigvalsh(op.to_matrix())[0]) + float(e_core)

    assert abs(e_old_total - float(cas.e_tot)) > 1e-3, (
        "Original compactor unexpectedly reproduced FCI -- regression marker is stale.")
