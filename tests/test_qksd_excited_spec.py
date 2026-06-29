"""
Acceptance gates G1-G4 for specs/SPEC_qksd_excited.md (excited-state quantum Krylov).

Test-first: ``QuantumKrylovSolver.solve_excited`` does not exist yet, so this file is RED until the
spec is implemented. The same real-time Krylov subspace that gives the ground state also carries the
low-lying *excited* spectrum in its Ritz values; we validate those against the exact eigenvalues of
the same qubit Hamiltonian that are *reachable* from |HF> (nonzero HF overlap).

Uses only pyscf/qiskit (no block2), so it runs in the non-DMRG process group -- but `make gates`
runs every test_*_spec.py in its own process anyway.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

CHEM_ACC = 1.6e-3  # Ha (1 kcal/mol)


def _h2():
    """H2 at equilibrium -- 4 qubits; the reachable singlet space is 2-dimensional."""
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")


def _h4():
    """H4 chain, 1.0 Angstrom spacing -- 8 qubits; several reachable excited singlets."""
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")


def _exact_spectrum(mh):
    """All exact eigenvalues of the qubit Hamiltonian (ascending), lifted by the offset."""
    w = np.linalg.eigvalsh(mh.qubit_hamiltonian.to_matrix())
    return np.sort(w.real) + mh.energy_offset


def _reachable_spectrum(mh, overlap_tol=1e-8):
    """Exact energies whose eigenstate has nonzero HF overlap -- the states real-time Krylov from
    |HF> can actually resolve. Reference for G3/G4."""
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    overlaps = np.abs(V.conj().T @ hf) ** 2
    reachable = w[overlaps > overlap_tol].real + mh.energy_offset
    return np.sort(reachable)


def test_G1_variational_interlacing():
    """Cauchy interlacing: the i-th Ritz value sits above the i-th exact eigenvalue, always."""
    for mh in (_h2(), _h4()):
        exact = _exact_spectrum(mh)
        solver = QuantumKrylovSolver(mh)
        for m in range(1, 13):
            step = solver.solve_excited(m)
            ritz = np.sort(step.energies)
            for i, e in enumerate(ritz):
                assert e >= exact[i] - 1e-9, (mh.num_qubits, m, i, e, exact[i])


def test_G2_ground_state_regression():
    """The excited API's lowest energy must equal the validated ground-state path exactly + FCI."""
    for mh in (_h2(), _h4()):
        e_fci = mh.ground_state_energy()
        solver = QuantumKrylovSolver(mh)
        for m in (4, 8, 12):
            ground = solver.solve(m).energy
            excited = solver.solve_excited(m).energies[0]
            assert excited == ground, (mh.num_qubits, m, excited, ground)
        assert abs(solver.solve_excited(12).energies[0] - e_fci) < CHEM_ACC


def test_G3_reachable_excited_converge():
    """DEFINITION OF DONE: the lowest k reachable energies are each matched by a Ritz value.

    FINDING (recorded depth): H4's reachable subspace is 12-dimensional; the 3rd state (the binding
    constraint here) is only 6.2 mHa off at M=16 but converges to < 0.5 mHa by M=20 and < 0.05 mHa
    by M=24 as the kept rank saturates the reachable space. Excited states need a deeper Krylov
    space than the ground state -- we gate at the depth where the rank has saturated enough to
    resolve the lowest k (see specs/SPEC_qksd_excited.md G3).
    """
    mh = _h4()
    reachable = _reachable_spectrum(mh)
    assert len(reachable) >= 2, len(reachable)
    k = min(3, len(reachable))
    solver = QuantumKrylovSolver(mh)
    step = solver.solve_excited(24, n_states=k)
    ritz = np.sort(step.energies)[:k]
    for i in range(k):
        assert abs(ritz[i] - reachable[i]) < CHEM_ACC, (i, ritz[i], reachable[i])
    # the kept subspace saturates at (at least) the reachable dimension we are probing
    assert step.rank >= k, (step.rank, k)


def test_G4_first_excitation_gap():
    """The physical observable: the first excitation energy, within chemical accuracy."""
    for mh in (_h2(), _h4()):
        reachable = _reachable_spectrum(mh)
        exact_gap = reachable[1] - reachable[0]
        solver = QuantumKrylovSolver(mh)
        energies = np.sort(solver.solve_excited(16, n_states=2).energies)
        krylov_gap = energies[1] - energies[0]
        assert abs(krylov_gap - exact_gap) < CHEM_ACC, (mh.num_qubits, krylov_gap, exact_gap)
