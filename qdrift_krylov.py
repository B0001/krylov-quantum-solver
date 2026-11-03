#!/usr/bin/env python3
"""
qDRIFT real-time Krylov -- the correct version of the step the original core got wrong.

Real-time Krylov needs |phi_k> = e^{-i k dt H}|phi_HF>. qDRIFT (Campbell, PRL 123, 070503, 2019)
approximates each e^{-i dt H} by L single-Pauli rotations e^{-i tau sgn(c_j) P_j}, term j drawn
with probability |c_j|/lambda, tau = lambda*dt/L, lambda = sum_j |c_j|. It is a random channel:
one realization is NOT e^{-i dt H}, only the average is. So here one realization supplies a whole
basis {|phi_0>, ..., |phi_{M-1}>} (one random walk, states read off every L rotations), and the
subspace matrices are averaged over R independent realizations:

    S = (1/R) sum_r V_r^dag V_r,   H = (1/R) sum_r V_r^dag H V_r.

WHY THIS IS STILL VARIATIONAL. Stack W = [V_1; ...; V_R]/sqrt(R). Then S = W^dag W and
H = W^dag (I_R (x) H) W -- a genuine Gram/Rayleigh pair in an R-fold direct sum whose spectrum is
H's. Every Ritz value is therefore >= the lowest eigenvalue of H over the space the walk reaches
(exact arithmetic; the thresholded solve keeps it numerically). Averaging the matrix elements of
SEPARATELY drawn walks per (i, j) entry would lose this property; this module does not do that.

THE FINDING (SPEC_qdrift_krylov G4). A single Pauli string does not conserve electron number, so
each realization LEAKS out of the N-electron sector (H2O CAS(4,4): ~11-15% of the norm at L=50,
~0.4-0.6% at L=2000; zero for CAS(2,2) where every rotation happens to stay in-sector). The
guaranteed floor is therefore the FOCK-SPACE minimum of H, which coincides with the N-electron
FCI for the systems gated here but not in general. `postselect=True` projects each state onto the
reference's (N_alpha, N_beta) sector (the classical analogue of a particle-number symmetry check)
and restores the N-electron FCI as a rigorous floor.

HONEST SCOPE: exact statevector simulation of the qDRIFT channel, no shot or device noise;
minimal-basis small active spaces; this reproduces a known method correctly -- it is not new and
shows no quantum advantage. Sample-count bound N = ceil(2 lambda^2 t^2 / eps) is Campbell's
diamond-norm bound for the averaged channel, quoted, not derived here.

See specs/SPEC_qdrift_krylov.md.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import solve_generalized_eig


@dataclass
class QDriftKrylovResult:
    energy: float          # lowest Ritz value + energy_offset (total energy, Ha)
    rank: int              # Krylov directions kept by the thresholded solve
    H: np.ndarray          # realization-averaged subspace Hamiltonian (electronic, no offset)
    S: np.ndarray          # realization-averaged overlap
    max_leakage: float     # largest out-of-sector norm^2 of any basis state (before postselect)
    dt: float
    tau: float             # rotation angle per qDRIFT sample
    lam: float             # 1-norm of the Pauli coefficients


def qdrift_sample_count(lam: float, t: float, eps: float) -> int:
    """Campbell's qDRIFT sample count N = ceil(2 lambda^2 t^2 / eps) for channel error eps.

    Scales as lambda^2 t^2 / eps -- NOT lambda^2 / eps^2 (the error the pasted skill template made).
    """
    return int(np.ceil(2.0 * lam**2 * t**2 / eps))


def sector_mask(mh: MolecularHamiltonian) -> np.ndarray:
    """Boolean mask of computational basis states with the reference's (N_alpha, N_beta).

    Qiskit-Nature's Jordan-Wigner ordering puts alpha spin-orbitals on qubits 0..n_orb-1 (low
    bits) and beta on n_orb..2n_orb-1.
    """
    n_orb = mh.num_spatial_orbitals
    idx = np.arange(2**mh.num_qubits)
    popcount = np.vectorize(lambda x: int(x).bit_count())
    n_alpha = popcount(idx & ((1 << n_orb) - 1))
    n_beta = popcount(idx >> n_orb)
    return (n_alpha == mh.num_particles[0]) & (n_beta == mh.num_particles[1])


def check_reference(mh: MolecularHamiltonian, psi0: np.ndarray, tol: float = 1e-8) -> None:
    """Refuse a reference state outside the N-electron sector (e.g. the |00..0> vacuum)."""
    out = float(np.sum(np.abs(psi0[~sector_mask(mh)]) ** 2))
    if out > tol:
        raise ValueError(
            f"reference has {out:.3g} of its norm outside the (N_alpha, N_beta)="
            f"{mh.num_particles} sector -- time evolution conserves N, so a Krylov space built "
            "from it cannot reach the N-electron ground state (use mh.hf_state())."
        )


def qdrift_krylov(
    mh: MolecularHamiltonian,
    krylov_dim: int = 6,
    samples_per_step: int = 400,
    realizations: int = 20,
    dt: float | None = None,
    threshold: float = 1e-6,
    postselect: bool = False,
    seed: int | None = None,
    psi0: np.ndarray | None = None,
) -> QDriftKrylovResult:
    """qDRIFT real-time Krylov ground-state estimate (see module docstring)."""
    op = mh.qubit_hamiltonian
    H_sparse = op.to_matrix(sparse=True)
    psi0 = mh.hf_state().data if psi0 is None else np.asarray(psi0, dtype=complex)
    check_reference(mh, psi0)

    if dt is None:  # pi / spectral width, as QuantumKrylovSolver does (dense: small spaces only)
        w = np.linalg.eigvalsh(op.to_matrix())
        dt = float(np.pi / (w[-1] - w[0]))
    labels = op.paulis
    coeffs = np.asarray(op.coeffs)
    if np.max(np.abs(coeffs.imag)) > 1e-10:
        raise ValueError("Hamiltonian has complex Pauli coefficients; qDRIFT here needs real ones")
    weights = np.abs(coeffs.real)
    lam = float(weights.sum())
    probs = weights / lam
    signs = np.sign(coeffs.real)
    paulis = [SparsePauliOp(labels[j]).to_matrix(sparse=True) for j in range(len(labels))]
    tau = lam * dt / samples_per_step
    cos_t, sin_t = np.cos(tau), np.sin(tau)
    mask = sector_mask(mh)

    rng = np.random.default_rng(seed)
    M = krylov_dim
    H_acc = np.zeros((M, M), dtype=complex)
    S_acc = np.zeros((M, M), dtype=complex)
    max_leak = 0.0
    for _ in range(realizations):
        states = [psi0]
        x = psi0
        for _k in range(1, M):
            for j in rng.choice(len(labels), size=samples_per_step, p=probs):
                x = cos_t * x - 1j * signs[j] * sin_t * (paulis[j] @ x)  # e^{-i tau sgn P}
            states.append(x)
        V = np.array(states).T
        max_leak = max(max_leak, float(np.max(np.sum(np.abs(V[~mask]) ** 2, axis=0))))
        if postselect:
            V = V * mask[:, None]
        S_acc += V.conj().T @ V
        H_acc += V.conj().T @ (H_sparse @ V)
    S = 0.5 * (S_acc + S_acc.conj().T) / realizations
    H = 0.5 * (H_acc + H_acc.conj().T) / realizations
    e, rank = solve_generalized_eig(H, S, threshold)
    return QDriftKrylovResult(mh.total_energy(e), rank, H, S, max_leak, dt, tau, lam)


def sector_ground_energy(mh: MolecularHamiltonian) -> float:
    """Exact (N_alpha, N_beta)-sector ground energy (dense; reference/validation only)."""
    m = sector_mask(mh)
    Hd = mh.qubit_hamiltonian.to_matrix()
    return mh.total_energy(float(np.linalg.eigvalsh(Hd[np.ix_(m, m)])[0]))
