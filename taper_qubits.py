#!/usr/bin/env python3
"""
Z2 qubit tapering (Bravyi, Gambetta, Mezzacapo, Temme 2017).

Molecular Hamiltonians carry Z2 symmetries -- Pauli operators that commute with every term
(spin-up parity, spin-down parity, and combinations). Each independent symmetry lets you remove
one qubit without changing the spectrum in the relevant sector. For an active space this is a
free reduction in the qubit (and circuit) cost of every downstream method -- directly useful for
fitting a system like Nb3 onto limited hardware.

Method: find the symmetry group as the GF(2) kernel of the symplectic check matrix; for each
Z-type generator, a Clifford U = (Z_sym + X_pivot)/sqrt(2) maps it to a single-qubit X; a
Hadamard turns that into Z; the qubit is then fixed to the sector eigenvalue (read from the HF
state) and traced out. The result is a smaller Hamiltonian with the same ground energy.

Validated by: tapered ground eigenvalue == CASCI, with the qubit count reduced (H2: 4 -> 1).
Same active-space interface as the rest of the stack: (h1, eri, e_core, nelec, norb).

Requires qubitization_blueprint.py and adapt_vqe.py in the same directory.
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qubitization_blueprint import build_qubit_hamiltonian, pauli_decompose, _kron, PAULI
from adapt_vqe import hf_state

_HAD = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)


def _symplectic(terms, n):
    """Pauli labels -> binary (X, Z) matrices. x=1 for X/Y, z=1 for Z/Y."""
    X = np.zeros((len(terms), n), int)
    Z = np.zeros((len(terms), n), int)
    for i, (lab, _) in enumerate(terms):
        for j, ch in enumerate(lab):
            if ch in "XY":
                X[i, j] = 1
            if ch in "ZY":
                Z[i, j] = 1
    return X, Z


def _gf2_kernel(M):
    """Basis for the null space of M over GF(2)."""
    M = M.copy() % 2
    rows, cols = M.shape
    pivots = {}
    r = 0
    for c in range(cols):
        piv = next((rr for rr in range(r, rows) if M[rr, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for rr in range(rows):
            if rr != r and M[rr, c]:
                M[rr] = (M[rr] + M[r]) % 2
        pivots[c] = r
        r += 1
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fc in free:
        v = np.zeros(cols, int)
        v[fc] = 1
        for c, pr in pivots.items():
            v[c] = M[pr, fc]
        basis.append(v % 2)
    return basis


def _gf2_independent(A):
    """Keep a maximal GF(2)-independent subset of the rows of A."""
    reduced, kept = [], []
    for r in A % 2:
        cur = r.copy()
        for b in reduced:
            piv = int(np.argmax(b))
            if cur[piv]:
                cur = (cur + b) % 2
        if cur.any():
            reduced.append(cur)
            kept.append(r)
    return kept


def find_symmetries(terms, n):
    """Pauli (x|z) operators commuting with every term: kernel of [Z | X] over GF(2).
    Returns list of (x_part, z_part)."""
    X, Z = _symplectic(terms, n)
    ker = _gf2_kernel(np.hstack([Z, X]))
    return [(v[:n], v[n:]) for v in ker if v.any()]


def _drop_qubit(M, v, q, n, keep):
    """Project qubit q onto |keep> and trace it out (M must commute with Z_q)."""
    M4 = M.reshape([2] * n + [2] * n)
    M4 = np.take(M4, keep, axis=q)
    M4 = np.take(M4, keep, axis=(n - 1) + q)   # ket axis shifted left by 1 after first take
    newn = n - 1
    v4 = np.take(v.reshape([2] * n), keep, axis=q)
    return M4.reshape(2 ** newn, 2 ** newn), v4.reshape(2 ** newn)


def taper_hamiltonian(h1, eri, e_core, nelec, norb):
    """Taper the active-space qubit Hamiltonian. Returns a dict with the reduced Hamiltonian,
    qubit counts, removed-qubit list, and the (validated) ground energy."""
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    hf = hf_state(nelec[0], nelec[1], n)
    terms = pauli_decompose(H, n)
    syms = find_symmetries(terms, n)
    z_syms = [z for (x, z) in syms if not x.any() and z.any()]

    H_cur, hf_cur, ncur, removed, qubit_ids = H, hf, n, [], list(range(n))
    if z_syms:
        gens = [r.copy() for r in _gf2_independent(np.array(z_syms))]
        while gens:
            g = gens[0]
            m = len(qubit_ids)
            cand = [q for q in range(m) if g[q] == 1]
            pivot = next((q for q in cand if all(o[q] == 0 for o in gens[1:])), cand[0])
            Zg = _kron([PAULI["Z"] if g[q] else PAULI["I"] for q in range(m)])
            Xp = _kron([PAULI["X"] if q == pivot else PAULI["I"] for q in range(m)])
            U = (Zg + Xp) / np.sqrt(2)                      # Clifford: Z_g <-> X_pivot
            H_cur = U @ H_cur @ U.conj().T
            sector = int(np.round(np.real(np.vdot(hf_cur, Zg @ hf_cur))))   # +-1 from HF
            hf_cur = U @ hf_cur
            Had = _kron([_HAD if q == pivot else PAULI["I"] for q in range(m)])
            H_cur, hf_cur = Had @ H_cur @ Had, Had @ hf_cur
            H_cur, hf_cur = _drop_qubit(H_cur, hf_cur, pivot, m, keep=0 if sector == 1 else 1)
            removed.append(qubit_ids.pop(pivot))
            ncur -= 1
            gens = [np.delete(o, pivot) for o in gens[1:]]

    ground = float(np.linalg.eigvalsh(H_cur)[0].real) + e_core
    return {"H_tapered": H_cur, "n_qubits_original": n, "n_qubits_tapered": ncur,
            "qubits_removed": removed, "ground_energy": ground}


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo

    def reference(atom, norb, ne, basis="sto-3g"):
        mol = gto.M(atom=atom, basis=basis)
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.kernel()
        cas = mcscf.CASCI(mf, norb, ne)
        cas.verbose = 0
        cas.kernel()
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        return h1, eri, float(e_core), (ne // 2, ne // 2), float(cas.e_tot)

    systems = {
        "H2  CAS(2,2)": ("H 0 0 0; H 0 0 0.74", 2, 2),
        "LiH CAS(2,2)": ("Li 0 0 0; H 0 0 1.6", 2, 2),
        "H4  CAS(4,4)": ("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4),
    }
    print("=" * 72)
    for label, (atom, norb, ne) in systems.items():
        h1, eri, e_core, nelec, casci = reference(atom, norb, ne)
        r = taper_hamiltonian(h1, eri, e_core, nelec, norb)
        ok = abs(r["ground_energy"] - casci) < 1e-7
        print(f"[{label}] {r['n_qubits_original']} -> {r['n_qubits_tapered']} qubits "
              f"(removed {len(r['qubits_removed'])})  CASCI={casci:.8f}  "
              f"tapered={r['ground_energy']:.8f}  match={ok}")
        assert ok, "tapering changed the ground energy"
    print("=" * 72)
    print("Z2 tapering removes qubits for free -- same ground energy, smaller register for every")
    print("downstream method (SQD, Krylov, ADAPT, qubitization).")
