"""
Acceptance gates G1-G4 for specs/SPEC_skqd.md (sample-based Krylov quantum diagonalization).

Test-first: ``hybrid_quantum_solver.skqd`` does not exist yet, so this file is RED until the spec
is implemented. SKQD samples Slater determinants from the real-time-evolved Krylov states
|Psi_k> = e^(-i k dt H)|HF>, unions them into a determinant subspace, and diagonalizes H there.
The ground truth is the validated exact-evolution QuantumKrylovSolver and PySCF FCI.

Uses only pyscf/qiskit (no block2), so it runs in the non-DMRG process group -- but `make gates`
runs every test_*_spec.py in its own process anyway.
"""
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from hybrid_quantum_solver.skqd import SampleKrylovSolver  # noqa: F401  (RED until implemented)

CHEM_ACC = 1.6e-3  # Ha (1 kcal/mol)


def _h4():
    """H4 chain, 1.0 Angstrom spacing -- 8 qubits, FCI-trivial, ground state concentrated."""
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")


def _n2_cas66():
    """N2 near equilibrium in CAS(6,6) -- 12 qubits, the multireference rung."""
    return build_molecular_hamiltonian(
        atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6
    )


def test_G1_variational_floor():
    """A sampled-subspace Rayleigh quotient can never beat FCI; a violation is a mapping bug."""
    for mh in (_h4(), _n2_cas66()):
        e_fci = mh.ground_state_energy()
        for depth in (4, 8):
            for n_shots in (2_000, 20_000):
                step = SampleKrylovSolver(
                    mh, n_shots=n_shots, depth=depth, seed=0
                ).solve()
                assert step.energy >= e_fci - 1e-6, (mh.num_qubits, depth, n_shots,
                                                     step.energy, e_fci)


def test_G2_converges_to_fci():
    """DEFINITION OF DONE: H4 reaches chemical accuracy at depth >= 6 with a large shot budget."""
    mh = _h4()
    e_fci = mh.ground_state_energy()
    step = SampleKrylovSolver(mh, n_shots=50_000, depth=8, seed=1).solve()
    assert abs(step.energy - e_fci) < CHEM_ACC, (step.energy, e_fci)


def test_G3_agrees_with_exact_krylov():
    """At matched dt/depth and high shots, SKQD agrees with exact Krylov within chemical accuracy.

    FINDING (revised gate): SKQD is NOT bounded below by the dense Krylov estimate -- with enough
    samples its determinant subspace spans *more* than the depth-dimensional span{Psi_k}, so it
    reaches essentially full FCI, which lies *below* the finite-M Krylov energy. The real
    variational floor is FCI (checked here), and the two methods agree to < 1.6 mHa.
    """
    mh = _h4()
    e_fci = mh.ground_state_energy()
    exact = QuantumKrylovSolver(mh)
    e_exact = exact.solve(8).energy
    step = SampleKrylovSolver(mh, dt=exact.dt, n_shots=50_000, depth=8, seed=2).solve()
    assert step.energy >= e_fci - 1e-6, (step.energy, e_fci)            # true floor is FCI
    assert abs(step.energy - e_exact) < CHEM_ACC, (step.energy, e_exact)


def test_G4_monotone_in_shots():
    """More shots must not raise the energy beyond the sampling spread (it can only add dets)."""
    mh = _h4()
    schedule = (2_000, 10_000, 50_000)
    steps = SampleKrylovSolver(mh, depth=8, seed=3).convergence(schedule)
    energies = [s.energy for s in steps]
    n_dets = [s.n_dets for s in steps]
    tol = 5e-4  # Ha: allowance for finite-sample noise at fixed depth
    assert all(energies[i + 1] <= energies[i] + tol for i in range(len(energies) - 1)), energies
    assert all(n_dets[i + 1] >= n_dets[i] for i in range(len(n_dets) - 1)), n_dets
