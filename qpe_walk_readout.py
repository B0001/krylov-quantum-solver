#!/usr/bin/env python3
"""
Quantum Phase Estimation readout on the qubitization walk operator.

The walk operator W (from qubitization_blueprint) encodes the whole spectrum in its
eigenphases:  eig(W) = e^{+-i theta_k}  with  E_k = lambda * cos(theta_k).  QPE measures
theta_k on a t-bit phase register; decoding gives the energy. This module simulates textbook
QPE on an eigenstate input via the exact Fejer-kernel outcome distribution (no need to form
the giant QPE unitary) and quantifies the two things that actually gate FT chemistry:

  * PRECISION:  energy resolution ~ lambda / 2^t, so t ~ log2(lambda / epsilon) phase bits.
  * STATE PREP: a trial state's ground-state overlap |<g|psi>|^2 sets the QPE success
                probability. For small molecules HF overlap is ~1; the bottleneck is that
                this overlap can decay (potentially exponentially) for large, strongly
                correlated systems -- which is exactly why good references (CASSCF) and
                near-term methods (SQD, Krylov) matter for systems like Nb3.

Requires qubitization_blueprint.py in the same directory.
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qubitization_blueprint import build_qubit_hamiltonian, pauli_decompose


def qpe_distribution(phases, weights, t_bits):
    """Textbook QPE outcome distribution over 2^t bins for eigen-phases (in [0,1)) with weights."""
    Y = np.arange(2 ** t_bits)
    P = np.zeros(2 ** t_bits)
    for phi, w in zip(phases, weights):
        if w < 1e-14:
            continue
        delta = phi - Y / 2 ** t_bits
        num = np.sin(np.pi * 2 ** t_bits * delta) ** 2
        den = (2 ** t_bits) ** 2 * np.sin(np.pi * delta) ** 2
        kernel = np.where(np.abs(np.sin(np.pi * delta)) < 1e-12, 1.0, num / np.maximum(den, 1e-300))
        P += w * kernel
    return P / P.sum()


def hartree_fock_vector(norb, n_occ_spinorbitals, n_qubits):
    """HF determinant: occupy the lowest spin-orbitals. Qubit 0 = MSB, so bit (n-1-p)."""
    idx = sum(1 << (n_qubits - 1 - p) for p in range(n_occ_spinorbitals))
    v = np.zeros(2 ** n_qubits)
    v[idx] = 1.0
    return v


def run_qpe(h1, eri, norb, e_core, trial_vec, t_bits):
    """Estimate the ground energy via QPE on the walk operator for a given trial state.

    Returns (E_total_est, lambda, success_prob, ground_overlap).
    """
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    Ek, Vk = np.linalg.eigh(H)
    lam = sum(abs(c) for _, c in pauli_decompose(H, n))
    phases = np.arccos(np.clip(Ek / lam, -1, 1)) / (2 * np.pi)   # theta_k/(2pi) in [0, 0.5]
    weights = np.abs(Vk.conj().T @ trial_vec) ** 2
    P = qpe_distribution(phases, weights, t_bits)
    y = int(np.argmax(P))
    E_est = lam * np.cos(2 * np.pi * y / 2 ** t_bits)
    gbin = int(round(phases[0] * 2 ** t_bits))
    win = [b % 2 ** t_bits for b in range(gbin - 1, gbin + 2)]
    return E_est + e_core, lam, float(P[win].sum()), float(weights[0])


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo

    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, 2, 2)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), 2)
    H, n = build_qubit_hamiltonian(h1, eri, 2)
    _, Vk = np.linalg.eigh(H)
    ground = Vk[:, 0]
    hf = hartree_fock_vector(2, 2, n)   # 2 electrons in the lowest 2 spin-orbitals

    print("=" * 72)
    print(f"H2 STO-3G CAS(2,2)  CASCI ground = {cas.e_tot:.8f} Ha")
    print("-" * 72)
    print("Trial = EXACT ground state  (success prob -> 1): precision vs phase bits t")
    for t in [4, 6, 8, 10, 12]:
        E_est, lam, ps, _ = run_qpe(h1, eri, 2, e_core, ground, t)
        print(f"  t={t:>2}  E_est={E_est:.6f}  err={abs(E_est - cas.e_tot) * 1e3:8.3f} mHa  p_success={ps:.3f}")
    print("-" * 72)
    print("Trial = HARTREE-FOCK determinant  (realistic state prep)")
    for t in [8, 10, 12]:
        E_est, lam, ps, olap = run_qpe(h1, eri, 2, e_core, hf, t)
        print(f"  t={t:>2}  E_est={E_est:.6f}  err={abs(E_est - cas.e_tot) * 1e3:8.3f} mHa  "
              f"p_success={ps:.3f}  |<g|HF>|^2={olap:.3f}")
    print("-" * 72)
    print(f"lambda={lam:.4f}.  Energy = lambda*cos(2*pi*phi); precision ~ lambda/2^t.")
    print("HF overlap ~0.99 here, so QPE succeeds. The bottleneck bites when overlap")
    print("decays with system size -- the case for strongly correlated targets like Nb3.")
    print("=" * 72)
