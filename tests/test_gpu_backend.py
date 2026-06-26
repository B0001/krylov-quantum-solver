"""
CPU validation of the GPU backend's numerics.

The ``device="gpu"`` path differs from the validated CPU path in exactly two things: the arrays
live on the GPU (CuPy) and the real-time step uses ``expm_multiply_taylor`` instead of SciPy's
``expm_multiply`` (cupyx has no expm_multiply). CuPy is a drop-in array backend, so the only new
*numerics* to validate is ``expm_multiply_taylor`` -- which is backend-agnostic and therefore
fully testable on CPU here. The actual GPU execution must still be validated on an NVIDIA node.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import (
    QuantumKrylovSolver,
    expm_multiply_taylor,
    solve_generalized_eig,
)


def test_expm_multiply_taylor_matches_scipy_dense_hermitian():
    rng = np.random.default_rng(0)
    n = 64
    A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    H = sp.csr_matrix(A + A.conj().T)                      # Hermitian
    v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    t = 0.37
    lam = float(np.sum(np.abs(np.linalg.eigvalsh(H.toarray()))))   # safe ||H||_2 bound

    ref = expm_multiply(-1j * t * H, v)
    got = expm_multiply_taylor(H, v, t, lam)
    assert np.linalg.norm(got - ref) / np.linalg.norm(ref) < 1e-10


def test_expm_multiply_taylor_on_molecular_hamiltonian():
    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    psi0 = np.asarray(mh.hf_state().data, dtype=complex)
    lam = float(np.sum(np.abs(mh.qubit_hamiltonian.coeffs)))      # >= ||H||_2

    ref = expm_multiply(-1j * 0.5 * H, psi0)
    got = expm_multiply_taylor(H, psi0, 0.5, lam)
    assert np.linalg.norm(got - ref) / np.linalg.norm(ref) < 1e-9


def test_taylor_krylov_energy_matches_cpu_solver():
    """The GPU path builds its basis with the Taylor step; on CPU that must reproduce the
    scipy-expm CPU solver's ground-state estimate (this is the GPU path minus the array move)."""
    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    cpu = QuantumKrylovSolver(mh)
    e_cpu = cpu.convergence(8)[-1].energy

    H, psi0, dt = cpu._H, cpu._psi0, cpu.dt
    lam = float(np.sum(np.abs(mh.qubit_hamiltonian.coeffs)))
    basis = [psi0.copy()]
    for _ in range(7):
        basis.append(expm_multiply_taylor(H, basis[-1], dt, lam))

    B = np.array(basis)
    S = B.conj() @ B.T
    Hm = B.conj() @ H.dot(B.T)
    energy, _ = solve_generalized_eig(0.5 * (Hm + Hm.conj().T), 0.5 * (S + S.conj().T))
    assert abs((energy + mh.energy_offset) - e_cpu) < 1e-6
