#!/usr/bin/env python3
"""
Qubitization blueprint, verified on a real molecular Hamiltonian.

This is the fault-tolerant backbone the near-term methods (SQD, Krylov) do not touch. It
builds the qubitization *walk operator* W for a molecular Hamiltonian and verifies the
spectral guarantee numerically:

    eig(W) = e^{±i arccos(E_k / lambda)}   <=>   E_k = lambda * cos(theta_k),

where  lambda = sum_l |c_l|  is the block-encoding 1-norm of  H = sum_l c_l P_l.

Why this matters for FT design (each maps to an object below):
  * DATA LOADING  -> PREPARE loads the coefficients c_l into ancilla amplitudes
                     sqrt(|c_l|/lambda). This is the "load millions of numbers" bottleneck.
  * T-GATE BUDGET -> FT-QPE applied to W costs O(lambda / epsilon) walk steps. lambda is THE
                     cost driver; shrinking it (double factorization, tensor hypercontraction)
                     is the single biggest lever on T-count.
  * STATE PREP    -> QPE needs a trial state on the SYSTEM register (separate from the
                     block-encoding ancilla) with non-vanishing ground-state overlap. That
                     overlap sets the success probability and is the contested part of any
                     "quantum advantage" claim for chemistry.

Same active-space inputs as the SQD and Krylov modules: (h1, eri, norb).
Exact-matrix construction is for verification on small systems; the *circuit* version
replaces PREPARE/SELECT with gate sequences and never forms these matrices.
"""

import numpy as np
from itertools import product

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1.0, -1.0]).astype(complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


def _kron(mats):
    out = np.array([[1]], complex)
    for m in mats:
        out = np.kron(out, m)
    return out


def jw_annihilation(p, n):
    """a_p under Jordan-Wigner on n qubits: (Z_0..Z_{p-1}) ⊗ (X+iY)/2 ⊗ I..."""
    return _kron([Z] * p + [(X + 1j * Y) / 2] + [I2] * (n - p - 1))


def build_qubit_hamiltonian(h1, eri, norb):
    """Electronic H over 2*norb spin orbitals (block ordering: alpha=2i, beta=2i+1)."""
    nso = 2 * norb
    a = [jw_annihilation(p, nso) for p in range(nso)]
    ad = [op.conj().T for op in a]
    H = np.zeros((2 ** nso, 2 ** nso), complex)
    for p in range(norb):
        for q in range(norb):
            if abs(h1[p, q]) < 1e-12:
                continue
            for s in range(2):
                H += h1[p, q] * (ad[2 * p + s] @ a[2 * q + s])
    for p in range(norb):
        for q in range(norb):
            for r in range(norb):
                for u in range(norb):
                    v = eri[p, q, r, u]
                    if abs(v) < 1e-12:
                        continue
                    for s in range(2):
                        for t in range(2):
                            H += 0.5 * v * (
                                ad[2 * p + s] @ ad[2 * r + t]
                                @ a[2 * u + t] @ a[2 * q + s]
                            )
    return 0.5 * (H + H.conj().T), nso


def pauli_decompose(H, n, tol=1e-10):
    """H = sum_l c_l P_l over all n-qubit Pauli strings (exponential; verification only)."""
    terms = []
    for labels in product("IXYZ", repeat=n):
        c = np.trace(_kron([PAULI[ch] for ch in labels]).conj().T @ H) / (2 ** n)
        if abs(c.imag) > 1e-8:
            raise ValueError("non-Hermitian Pauli coefficient")
        if abs(c.real) > tol:
            terms.append(("".join(labels), float(c.real)))
    return terms


def build_walk_operator(terms, n):
    """LCU block encoding (PREPARE, SELECT) + qubitization walk W = R_anc · SELECT.

    Returns (W, lambda, n_terms, n_ancilla).
    """
    coeffs = np.array([c for _, c in terms])
    signs = np.sign(coeffs)
    absc = np.abs(coeffs)
    lam = float(absc.sum())
    L = len(terms)
    a = int(np.ceil(np.log2(L)))
    A = 2 ** a
    sys = 2 ** n

    # PREPARE: ancilla unitary whose first column is sqrt(|c_l|/lambda) (zero-padded)
    amp = np.zeros(A)
    amp[:L] = np.sqrt(absc / lam)
    M = np.eye(A, dtype=complex)
    M[:, 0] = amp
    prep, _ = np.linalg.qr(M)
    prep[:, 0] = amp  # enforce the loaded column exactly (QR can flip its sign)

    # SELECT: sum_l |l><l| ⊗ (sign_l P_l); unused indices act as identity on the system
    SELECT = np.zeros((A * sys, A * sys), complex)
    for l in range(A):
        proj = np.zeros((A, A), complex)
        proj[l, l] = 1.0
        if l < L:
            Pl = _kron([PAULI[ch] for ch in terms[l][0]]) * signs[l]
        else:
            Pl = np.eye(sys, dtype=complex)
        SELECT += np.kron(proj, Pl)

    prep0 = prep[:, 0]
    R_anc = 2.0 * np.outer(prep0, prep0.conj()) - np.eye(A)
    W = np.kron(R_anc, np.eye(sys, dtype=complex)) @ SELECT
    return W, lam, L, a


def verify_qubitization(h1, eri, norb, e_core=0.0, casci_energy=None):
    """Build W and confirm lambda*cos(theta_k) reproduces the Hamiltonian spectrum."""
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    eigH = np.linalg.eigvalsh(H).real
    terms = pauli_decompose(H, n)
    W, lam, L, a = build_walk_operator(terms, n)

    theta = np.angle(np.linalg.eigvals(W))
    recovered = np.unique(np.round(lam * np.cos(theta), 8))
    true = np.unique(np.round(eigH, 8))
    max_err = max(np.min(np.abs(recovered - e)) for e in true)

    info = {
        "n_qubits_system": n, "n_pauli_terms": L, "n_ancilla": a,
        "lambda_1norm": lam, "spectral_max_error": float(max_err),
        "ground_energy_total": float(eigH.min() + e_core),
    }
    if casci_energy is not None:
        info["casci_match"] = abs(eigH.min() + e_core - casci_energy) < 1e-8
    return info


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

    info = verify_qubitization(h1, eri, 2, e_core=e_core, casci_energy=cas.e_tot)
    print("=" * 70)
    print("Qubitization walk operator for H2 / STO-3G, CAS(2,2)")
    print("-" * 70)
    for k, v in info.items():
        print(f"  {k:>22}: {v}")
    print("-" * 70)
    print("  lambda is the FT-QPE cost driver: walk steps ~ O(lambda / epsilon).")
    print("  Reducing lambda (double factorization) directly cuts the T-gate budget.")
    print("=" * 70)
