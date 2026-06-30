#!/usr/bin/env python3
"""
Hamiltonian-moment ground-state energy: the PDS energy functional and the connected-moment
expansion (CMX).

A different family from subspace diagonalization: estimate the ground-state energy purely from the
Hamiltonian *moments* mu_n = <phi|H^n|phi> of a fixed reference |phi> (here Hartree-Fock), with no
time evolution and no eigenproblem in Hilbert space.

* **PDS(K)** (Peeters-Devreese-Soldatov; Peng & Kowalski, Quantum 5, 473, 2021, arXiv:2101.08526).
  Solve the K x K linear system ``M X = -Y`` with ``M_ij = <H^(2K-i-j)>`` and ``Y_i = <H^(2K-i)>``,
  then the ground-state energy is the **smallest root** of ``P_K(E) = E^K + sum_i X_i E^(K-i)``.
  PDS(K) is a *variational upper bound*: ``E_0 <= min root <= <H>`` (PDS(1) = <H> exactly). The bound
  tightens with K.

* **CMX(2)** (connected-moment / Horn-Weinstein expansion; Kowalski & Peng, J. Chem. Phys. 153,
  201102, 2020). ``E = I_1 - I_2^2 / I_3`` with the connected moments (cumulants)
  ``I_1 = mu_1``, ``I_2 = mu_2 - mu_1^2``, ``I_3 = mu_3 - 3 mu_1 mu_2 + 2 mu_1^3``. CMX is **not**
  variational and can dip below E_0 -- the honest contrast to PDS (see specs/SPEC_moment_pds.md).

Reuses the validated qubit Hamiltonian and Hartree-Fock reference (``MolecularHamiltonian``); moments
are computed by repeated sparse mat-vec. The reference for every result is FCI (dense diagonalization
of the same operator).
"""
from __future__ import annotations

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian


def hamiltonian_moments(mh: MolecularHamiltonian, max_order: int):
    """Hamiltonian moments ``mu_n = <phi|H^n|phi>`` (n = 0..max_order) of the Hartree-Fock reference.

    Returns ``(moments, offset)`` with ``moments`` in the electronic frame (``mu_0 = 1``) and
    ``offset = energy_offset``; a physical energy estimate adds ``offset`` to the electronic-frame
    result. Moments come from repeated sparse mat-vec (``w_k = H w_{k-1}``), real because H is
    Hermitian.
    """
    if max_order < 1:
        raise ValueError("max_order must be >= 1")
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    phi = np.asarray(mh.hf_state().data, dtype=complex)
    moments = [1.0]
    w = phi.copy()
    for _ in range(max_order):
        w = H.dot(w)
        moments.append(float(np.real(phi.conj() @ w)))
    return np.array(moments), float(mh.energy_offset)


def pds_energy(moments: np.ndarray, order: int, offset: float = 0.0) -> float:
    """PDS(K) variational ground-state energy: the smallest root of ``P_K(E)``.

    Needs moments up to ``2*order - 1``. ``offset`` lifts the electronic-frame result to a physical
    energy. Variational: the returned value is ``>= E_0`` (up to numerical noise).
    """
    K = int(order)
    if len(moments) < 2 * K:
        raise ValueError(f"PDS({K}) needs moments up to H^{2*K-1} ({2*K} entries), got {len(moments)}")
    M = np.array([[moments[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])
    Y = np.array([moments[2 * K - i] for i in range(1, K + 1)])
    X = np.linalg.solve(M, -Y)
    poly = np.concatenate([[1.0], X])              # E^K + X_1 E^(K-1) + ... + X_K
    roots = np.roots(poly)
    real_roots = roots[np.abs(roots.imag) < 1e-6].real
    if real_roots.size == 0:
        raise ValueError(f"PDS({K}) produced no real root (ill-conditioned moments?)")
    return float(real_roots.min()) + offset


def cmx2_energy(moments: np.ndarray, offset: float = 0.0) -> float:
    """CMX(2) connected-moment energy ``I_1 - I_2^2/I_3`` (non-variational; for contrast with PDS)."""
    mu = moments
    i1 = mu[1]
    i2 = mu[2] - mu[1] ** 2
    i3 = mu[3] - 3.0 * mu[1] * mu[2] + 2.0 * mu[1] ** 3
    return float(i1 - i2 ** 2 / i3) + offset


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H2": dict(atom="H 0 0 0; H 0 0 0.74"),
        "H4": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
        "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    }
    for name, kw in cases.items():
        mh = build_molecular_hamiltonian(**kw)
        fci = mh.ground_state_energy()
        mu, off = hamiltonian_moments(mh, 7)
        print(f"{name}: FCI={fci:.6f}  <H>={mu[1] + off:.6f}")
        for K in (1, 2, 3, 4):
            e = pds_energy(mu, K, off)
            print(f"   PDS({K})={e:.6f}  err={(e - fci) * 1e3:+.3f} mHa  variational={e >= fci - 1e-9}")
        e_cmx = cmx2_energy(mu, off)
        print(f"   CMX(2)={e_cmx:.6f}  err={(e_cmx - fci) * 1e3:+.3f} mHa  variational={e_cmx >= fci - 1e-9}")
