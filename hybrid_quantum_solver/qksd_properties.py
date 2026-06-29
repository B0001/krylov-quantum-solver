#!/usr/bin/env python3
"""
qksd_properties.py -- molecular properties (dipoles, transition dipoles, oscillator strengths)
from quantum Krylov subspace eigenstates.

The energies are not the only thing the Krylov subspace carries: once the generalized
eigenproblem is solved, each Ritz eigenstate |Psi_m> = sum_i C[i,m] |phi_i> is a genuine
Hilbert-space vector, and the expectation/transition value of any operator O between two of them
is just <Psi_m| O |Psi_n>. This is the classically-checkable face of the QKSD property framework
of Oumarou et al., *Molecular Properties from Quantum Krylov Subspace Diagonalization*,
arXiv:2501.05286 (whose RDM/QSP machinery measures exactly these matrix elements on hardware) --
here we evaluate them directly on the exact-evolution eigenstates from ``QuantumKrylovSolver``.

Ground truth: the same matrix elements between the exact (dense-diagonalized) eigenstates of the
same qubit operators -- the FCI reference for properties, in the same idiom that
``MolecularHamiltonian.ground_state_energy`` is the FCI reference for energies.

Inputs are the ``(energies, states)`` from ``QuantumKrylovSolver.eigenstates`` and the dipole
operators from ``build_dipole_operators`` (as matrices). See specs/SPEC_qksd_properties.md.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np


def property_matrix(states: np.ndarray, operator) -> np.ndarray:
    """Matrix of ``<Psi_m| O |Psi_n>`` over a set of states.

    Args:
        states: ``(k, N)`` array whose row ``m`` is the ket |Psi_m> in the computational basis.
        operator: an ``(N, N)`` operator -- dense ndarray or any scipy.sparse matrix (only
            ``operator @ vector`` is used, so a sparse ``SparsePauliOp.to_matrix(sparse=True)``
            works without densifying).

    Returns:
        ``(k, k)`` complex array ``O[m, n] = <Psi_m| O |Psi_n>``. Hermitian (up to round-off) for
        a Hermitian operator.
    """
    states = np.asarray(states)
    op_states = (operator @ states.T)                # (N, k); sparse matvec stays sparse
    return states.conj() @ np.asarray(op_states)


def transition_dipoles(states: np.ndarray, dipole_ops: Sequence) -> np.ndarray:
    """Stack of per-axis dipole matrices ``mu[axis][m, n] = <Psi_m| mu_axis |Psi_n>``.

    ``dipole_ops`` is the list of Cartesian dipole operators (e.g. from ``build_dipole_operators``,
    as matrices). The diagonal ``mu[:, m, m]`` is the permanent dipole vector of state ``m``; the
    off-diagonal ``mu[:, 0, n]`` is the transition dipole from the ground state to state ``n``.
    """
    return np.stack([property_matrix(states, op) for op in dipole_ops])


def oscillator_strengths(
    energies: Sequence[float], states: np.ndarray, dipole_ops: Sequence
) -> np.ndarray:
    """Dipole oscillator strengths ``f_n`` from the ground state to each state (length gauge).

        f_n = (2/3) (E_n - E_0) |<Psi_0| mu |Psi_n>|^2          (atomic units)

    with the squared transition-dipole magnitude summed over the three Cartesian axes. ``f_0 = 0``
    by construction. Energies in Hartree and dipoles in a.u. make ``f_n`` dimensionless; a
    dipole-forbidden (dark) transition gives ``f_n ~ 0``.
    """
    energies = np.asarray(energies, dtype=float)
    mu = transition_dipoles(states, dipole_ops)      # (3, k, k)
    trans = mu[:, 0, :]                              # (3, k): ground -> n per axis
    strength = np.sum(np.abs(trans) ** 2, axis=0)   # |mu_0n|^2
    return (2.0 / 3.0) * (energies - energies[0]) * strength
