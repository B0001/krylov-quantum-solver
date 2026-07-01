#!/usr/bin/env python3
"""
Quantum imaginary-time evolution (QITE) -- reach the ground state by approximating the non-unitary
imaginary-time step e^{-Delta tau H} with a *unitary* e^{-i Delta tau A} determined from measurements
(Motta et al., "Determining eigenstates and thermal states on a quantum computer using quantum
imaginary time evolution", Nature Physics 16, 205, 2020; arXiv:1901.07653).

Imaginary-time evolution |psi(beta)> = e^{-beta H}|psi> / ||.|| projects onto the ground state:
E(beta) = <psi(beta)|H|psi(beta)> decreases monotonically and variationally to E_0. QITE realises one
step on a quantum computer by finding the Hermitian A = sum_I a_I sigma_I (Paulis over a *domain*)
whose unitary reproduces the normalised imaginary-time step. Matching the imaginary-time derivative
-(H - <H>)|psi> to -i A|psi> by least squares gives the linear system

    S a = b,   S_IJ = Re<psi|sigma_I sigma_J|psi>,   b_I = Im<psi|sigma_I H|psi>,

and the step is |psi> <- e^{-i Delta tau A}|psi| (renormalised). With the *full* operator domain the
unitary reproduces exact imaginary-time evolution (up to the O(Delta tau^2) step error); a *truncated*
(local, low-weight) domain cannot -- QITE's accuracy is set by whether the domain spans the operators
that build the correlation (Motta's locality point). See specs/SPEC_qite.md.

Reuses the validated qubit Hamiltonian and Hartree-Fock reference (``MolecularHamiltonian``); the
reference is exact imaginary-time evolution and FCI (dense). Small systems: the full Pauli domain is
4^n operators, so this is exercised on H2 (n=4).
"""
from __future__ import annotations

from itertools import product

import numpy as np
from scipy.linalg import expm

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian

_P1 = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_operators(n_qubits: int, max_weight: int | None = None) -> np.ndarray:
    """All n-qubit Pauli matrices (optionally restricted to weight <= ``max_weight``).

    Returns an ``(n_ops, 2^n, 2^n)`` array. ``max_weight=None`` gives the full 4^n domain (exact QITE);
    a finite ``max_weight`` gives a truncated (local) domain.
    """
    ops = []
    for combo in product("IXYZ", repeat=n_qubits):
        weight = sum(1 for c in combo if c != "I")
        if max_weight is not None and weight > max_weight:
            continue
        mat = np.array([[1.0 + 0j]])
        for c in combo:
            mat = np.kron(mat, _P1[c])
        ops.append(mat)
    return np.array(ops)


def exact_imaginary_time(mh: MolecularHamiltonian, betas) -> list[float]:
    """Exact imaginary-time energies E(beta) = <psi(beta)|H|psi(beta)> from the HF reference.

    The variational, monotone-decreasing reference QITE must reproduce. Uses sparse
    ``expm_multiply``, so it works at any system size.
    """
    from scipy.sparse.linalg import expm_multiply

    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    psi = np.asarray(mh.hf_state().data, dtype=complex)
    offset = mh.energy_offset
    energies = []
    for beta in betas:
        v = expm_multiply(-beta * H, psi)
        v /= np.linalg.norm(v)
        energies.append(float(np.real(v.conj() @ (H @ v))) + offset)
    return energies


def qite_evolve(mh: MolecularHamiltonian, dtau: float, n_steps: int,
                operators: np.ndarray) -> list[float]:
    """Run QITE for ``n_steps`` of size ``dtau`` over the Pauli ``operators`` domain.

    Returns the total energy after each step. With the full domain this reproduces
    ``exact_imaginary_time`` up to O(dtau^2); a truncated domain plateaus where its operators can no
    longer lower the energy.
    """
    H = mh.qubit_hamiltonian.to_matrix()
    v = np.asarray(mh.hf_state().data, dtype=complex)
    offset = mh.energy_offset
    energies = []
    for _ in range(n_steps):
        vi = np.einsum("oij,j->oi", operators, v)          # sigma_I |psi>
        S = np.real(vi.conj() @ vi.T)                      # Re<psi|sigma_I sigma_J|psi>
        b = np.imag(vi.conj() @ (H @ v))                   # Im<psi|sigma_I H|psi>
        a, *_ = np.linalg.lstsq(S, b, rcond=1e-10)
        A = np.tensordot(a, operators, axes=([0], [0]))    # sum_I a_I sigma_I
        v = expm(-1j * dtau * A) @ v
        v /= np.linalg.norm(v)
        energies.append(float(np.real(v.conj() @ (H @ v))) + offset)
    return energies


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    fci = mh.ground_state_energy()
    print(f"H2: FCI={fci:.6f}  <H>_HF={mh.hf_energy:.6f}")
    full = pauli_operators(mh.num_qubits)
    e_qite = qite_evolve(mh, 0.1, 40, full)
    e_exact = exact_imaginary_time(mh, np.arange(1, 41) * 0.1)
    for k in (0, 9, 19, 39):
        b = (k + 1) * 0.1
        print(f"  beta={b:.1f}: exact ITE={e_exact[k]:+.6f}  QITE={e_qite[k]:+.6f}  "
              f"QITE-FCI={(e_qite[k] - fci) * 1e3:+.3f} mHa")
    print("  domain truncation:")
    for w in (2, 4):
        e = qite_evolve(mh, 0.05, 80, pauli_operators(mh.num_qubits, max_weight=w))
        print(f"    weight<={w}: QITE(beta=4)-FCI = {(e[-1] - fci) * 1e3:+.3f} mHa")
