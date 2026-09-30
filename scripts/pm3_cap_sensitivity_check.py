#!/usr/bin/env python3
"""
chem-pm3 follow-up: is the eta_bound LP's reported optimum materially LOOSE because the
Fourier-coefficient box (+-1000, see spec_pm3_subspace_eta_bound.py's _COEFF_CAP docstring) is
actively binding at the optimum? The main run's "[warning] coefficient cap active" fired at every
M on every witness (unlike the earlier isolated H4 sanity check, which wasn't cap-active) -- so
this needs to be checked directly on the real data, not inferred from the toy case.

Targets the borderline/kill witness (linear H6 R=1.1, d=3, M=20 -> G-b KILLED, eta_bound(M=24) =
0.0593, just above the 0.05 threshold) at M=24, and re-solves the LP at cap in {1e3, 1e4, 1e5,
1e6} to see whether a looser box changes the answer enough to flip the verdict.

Run: uv run python scripts/pm3_cap_sensitivity_check.py
"""
from __future__ import annotations

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from spec_pm3_subspace_eta_bound import _build_moment_bound_lp, _chain, _REACHABLE_TOL

W = dict(name="linear H6 R=1.1, d=3, M=20", atom=_chain(6, 1.1), d=3, m_floor=20)
M_TEST = 24
CAPS = [1e3, 1e4, 1e5, 1e6]


def main():
    mh = build_molecular_hamiltonian(atom=W["atom"])
    solver = QuantumKrylovSolver(mh)
    offset = mh.energy_offset
    u = np.asarray(mh.hf_state().data, dtype=complex)
    d = W["d"]

    Hd = mh.qubit_hamiltonian.to_matrix()
    wv, V = np.linalg.eigh(Hd)
    reach = np.where(np.abs(V.conj().T @ u) ** 2 > _REACHABLE_TOL)[0]
    reach_e = wv[reach]
    e_d_true = float(reach_e[d])

    energies, states = solver.eigenstates(W["m_floor"], n_states=d + 1)
    from hybrid_quantum_solver.certified_overlap import residual_norm
    Hs = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    centers = [energies[k] - offset for k in range(d + 1)]
    sigmas = [residual_norm(Hs, states[k], centers[k]) for k in range(d + 1)]
    beta_self = centers[d] - sigmas[d]

    print(f"witness: {W['name']}")
    print(f"e_d_true (electronic) = {e_d_true:.6f}, beta_self (electronic) = {beta_self:.6f}")

    solver._ensure_basis(M_TEST)
    basis = solver._basis[:M_TEST]
    moments = np.array([np.vdot(basis[0], basis[k]) for k in range(M_TEST)])
    dt = solver.dt
    e_min, e_max = float(wv.min()), float(wv.max())

    print(f"{'cap':>10} {'eta_bound':>12} {'cap_active':>11}")
    for cap in CAPS:
        bound = _build_moment_bound_lp(moments, e_min, e_max, dt, e_d_true, beta_self,
                                        coeff_cap=cap)
        print(f"{cap:>10.0e} {bound:>12.6f}")


if __name__ == "__main__":
    main()
