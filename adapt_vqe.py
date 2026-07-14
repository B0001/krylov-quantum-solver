#!/usr/bin/env python3
"""
ADAPT-VQE (exact statevector) -- the variational, gradient-driven member of the stack.

This is the principled version of the EfficientSU2 ansatz your hardware gateway used: instead
of a fixed circuit with guessed angles, ADAPT-VQE GROWS the ansatz one operator at a time,
each iteration adding the pool operator with the largest energy gradient, then re-optimizing
all parameters. The result is a compact, problem-tailored ansatz.

    |psi> = exp(theta_m G_m) ... exp(theta_1 G_1) |HF>,   G_k = T_k - T_k^dagger

with T_k the fermionic single/double excitations (anti-Hermitian generators -> unitary gates).
The ansatz subspace lies within the active-space CI space, so it is variational: the energy
approaches CASCI from above and can never fall below it.

Hardware mapping: a real device measures the gradients (commutator expectation values) and the
energy on the QPU; here we compute them exactly via the dense Hamiltonian, making this the
validation oracle for a hardware ADAPT-VQE run. Same active-space interface as the rest of the
stack: (h1, eri, e_core, nelec, norb).

Requires qubitization_blueprint.py (for the JW operators) in the same directory.
"""

import os
import sys
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qubitization_blueprint import build_qubit_hamiltonian, jw_annihilation


def hf_state(na, nb, n_qubits):
    """HF determinant statevector. Qubit index = 2*spatial + spin; qubit 0 = MSB."""
    occ = [2 * i for i in range(na)] + [2 * i + 1 for i in range(nb)]
    idx = sum(1 << (n_qubits - 1 - q) for q in occ)
    v = np.zeros(2 ** n_qubits, dtype=complex)
    v[idx] = 1.0
    return v


def build_pool(norb):
    """Generalized singles + doubles as anti-Hermitian generators G = T - T^dagger."""
    nso = 2 * norb
    a = [jw_annihilation(p, nso) for p in range(nso)]
    ad = [op.conj().T for op in a]
    pool = []
    for p in range(nso):
        for q in range(p):
            G = ad[p] @ a[q]
            G = G - G.conj().T
            if np.linalg.norm(G) > 1e-10:
                pool.append(G)
    seen = set()
    for p in range(nso):
        for q in range(p):
            for r in range(nso):
                for s in range(r):
                    key = tuple(sorted((p, q))) + tuple(sorted((r, s)))
                    if (p, q) <= (r, s) or key in seen:
                        continue
                    seen.add(key)
                    G = ad[p] @ ad[q] @ a[r] @ a[s]
                    G = G - G.conj().T
                    if np.linalg.norm(G) > 1e-10:
                        pool.append(G)
    return pool


def adapt_vqe(H, hf, pool, grad_tol=1e-3, max_ops=40):
    """Run ADAPT-VQE. Returns (psi, params, history) where history is a list of
    (n_ops, electronic_energy, selected_gradient) per growth step."""
    eig = [np.linalg.eigh(1j * G) for G in pool]   # iG is Hermitian -> real eigenvalues
    chosen, params, history = [], [], []
    psi = hf.copy()

    for _ in range(max_ops):
        Hpsi = H @ psi
        grads = [abs(2.0 * np.real(np.vdot(Hpsi, G @ psi))) for G in pool]
        k = int(np.argmax(grads))
        if grads[k] < grad_tol:
            break
        chosen.append(eig[k])
        params.append(0.0)

        def state(theta):
            s = hf.copy()
            for (lam, V), th in zip(chosen, theta):
                s = V @ (np.exp(-1j * th * lam) * (V.conj().T @ s))
            return s

        def energy(theta):
            s = state(theta)
            return float(np.real(np.vdot(s, H @ s)))

        res = minimize(energy, np.array(params), method="BFGS",
                       options={"gtol": 1e-7, "maxiter": 500})
        params = list(res.x)
        psi = state(params)
        history.append((len(chosen), res.fun, grads[k]))

    return psi, params, history


def fixed_order_vqe(H, hf, pool, order, max_ops):
    """The `adapt_vqe` comparison baseline (specs/SPEC_adapt_vqe_compactness.md): grows the SAME
    pool, re-optimizing all parameters the same way (BFGS, warm-started), but in a CALLER-GIVEN
    fixed ``order`` instead of by gradient -- isolates whether gradient-greedy selection buys a
    more compact ansatz, or a fixed order does just as well. Returns a history of
    (n_ops, electronic_energy) pairs."""
    eig = [np.linalg.eigh(1j * pool[k]) for k in order[:max_ops]]
    chosen, params, history = [], [], []

    for lam, V in eig:
        chosen.append((lam, V))
        params.append(0.0)

        def state(theta):
            s = hf.copy()
            for (lam, V), th in zip(chosen, theta):
                s = V @ (np.exp(-1j * th * lam) * (V.conj().T @ s))
            return s

        def energy(theta):
            s = state(theta)
            return float(np.real(np.vdot(s, H @ s)))

        res = minimize(energy, np.array(params), method="BFGS",
                       options={"gtol": 1e-7, "maxiter": 500})
        params = list(res.x)
        history.append((len(chosen), res.fun))

    return history


def adapt_ground_state(h1, eri, e_core, nelec, norb, grad_tol=1e-3, max_ops=40):
    """Total ground-state energy via ADAPT-VQE. total = <psi|H|psi> + e_core."""
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    hf = hf_state(nelec[0], nelec[1], n)
    pool = build_pool(norb)
    psi, params, history = adapt_vqe(H, hf, pool, grad_tol=grad_tol, max_ops=max_ops)
    e_total = float(np.real(np.vdot(psi, H @ psi))) + e_core
    return e_total, {"n_operators": len(params), "pool_size": len(pool),
                     "history": history}


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo

    def reference(atom, norb, ne, spin=0, basis="sto-3g"):
        mol = gto.M(atom=atom, basis=basis, spin=spin)
        mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
        mf.verbose = 0
        mf.kernel()
        na, nb = (ne + spin) // 2, (ne - spin) // 2
        cas = mcscf.CASCI(mf, norb, (na, nb))
        cas.verbose = 0
        cas.kernel()
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        return h1, eri, float(e_core), (na, nb), float(cas.e_tot)

    systems = {
        "H2  CAS(2,2)": ("H 0 0 0; H 0 0 0.74", 2, 2, 0),
        "LiH CAS(2,2)": ("Li 0 0 0; H 0 0 1.6", 2, 2, 0),
        "H4  CAS(4,4)": ("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4, 0),
    }
    print("=" * 70)
    for label, (atom, norb, ne, spin) in systems.items():
        h1, eri, e_core, nelec, casci = reference(atom, norb, ne, spin)
        e_adapt, info = adapt_ground_state(h1, eri, e_core, nelec, norb)
        assert e_adapt >= casci - 1e-6, "VARIATIONAL VIOLATION"
        print(f"[{label}] pool={info['pool_size']}  CASCI={casci:.6f}  "
              f"ADAPT={e_adapt:.6f}  ops={info['n_operators']}  "
              f"Δ={abs(e_adapt - casci) * 1e3:.4f} mHa")
        for n_ops, e_el, g in info["history"][:6]:
            print(f"     +op {n_ops:>2}:  E={e_el + e_core:.6f}  (selected grad {g:.2e})")
    print("=" * 70)
    print("ADAPT-VQE converges to CASCI from above, growing a compact problem-tailored ansatz.")
