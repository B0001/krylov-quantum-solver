#!/usr/bin/env python3
"""
chem-pm3 follow-up: test_pm3_eta_bound_spec.py caught that _build_moment_bound_lp's pointwise
grid constraint (f(E_j) >= indicator(E_j) at n_grid samples, not truly "for all E") lets the LP
find a polynomial that dips below the band indicator BETWEEN grid points -- an under-resolved
grid reports a bound that is too LOW (a false "tighter than真" result), not too high. Confirmed on
a synthetic case: n_grid=4000 underestimated eta by 0.0028 (bound 0.1472 vs true 0.1500);
n_grid=400000 converged to 0.150000. Direction matters here: as the grid refines, the reported
bound only INCREASES toward the true continuous-LP optimum -- so under-resolution biases in favor
of "SURVIVES" (bound <= 0.05), which is exactly the risky direction for chem-pm3's borderline
witnesses (main run used n_grid=12000).

This reruns the LP at the three borderline (M, witness) pairs from the main sweep at n_grid in
{12000 (as run), 60000, 200000} to check whether any G-b verdict flips once the grid is properly
resolved:
  - witness 1 (linear H6 R=1.0), M=12: eta_bound=0.0466 (SURVIVES, margin 0.0034 -- thin)
  - witness 2 (linear H6 R=1.1),  M=24: eta_bound=0.0593 (KILLED already -- refining can only
    raise the reported bound further, so this can only become MORE killed, not less; checked
    anyway for completeness)
  - witness 3 (asym H6),          M=16: eta_bound=0.0443 (SURVIVES, margin 0.0057 -- thin)

Run: uv run python scripts/pm3_grid_convergence_check.py
"""
from __future__ import annotations

import numpy as np

from hybrid_quantum_solver.certified_overlap import residual_norm
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from spec_pm3_subspace_eta_bound import _asym_chain, _build_moment_bound_lp, _chain, _REACHABLE_TOL

CASES = [
    dict(name="linear H6 R=1.0, d=3, M=16 floor -- LP at M=12",
         atom=_chain(6, 1.0), d=3, m_floor=16, m_lp=12),
    dict(name="linear H6 R=1.1, d=3, M=20 floor -- LP at M=24",
         atom=_chain(6, 1.1), d=3, m_floor=20, m_lp=24),
    dict(name="asym H6, d=3, M=6 floor -- LP at M=16",
         atom=_asym_chain([0.9, 0.9, 2.2, 0.9, 0.9]), d=3, m_floor=6, m_lp=16),
]

N_GRIDS = [12000, 60000, 200000]


def run_case(c):
    print("=" * 100)
    print(c["name"])
    mh = build_molecular_hamiltonian(atom=c["atom"])
    solver = QuantumKrylovSolver(mh)
    offset = mh.energy_offset
    u = np.asarray(mh.hf_state().data, dtype=complex)
    d = c["d"]

    Hd = mh.qubit_hamiltonian.to_matrix()
    wv, V = np.linalg.eigh(Hd)
    reach = np.where(np.abs(V.conj().T @ u) ** 2 > _REACHABLE_TOL)[0]
    reach_e = wv[reach]
    e_d_true = float(reach_e[d])

    Hs = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    energies, states = solver.eigenstates(c["m_floor"], n_states=d + 1)
    centers = [energies[k] - offset for k in range(d + 1)]
    sigmas = [residual_norm(Hs, states[k], centers[k]) for k in range(d + 1)]
    beta_self = centers[d] - sigmas[d]

    m_lp = c["m_lp"]
    solver._ensure_basis(m_lp)
    basis = solver._basis[:m_lp]
    moments = np.array([np.vdot(basis[0], basis[k]) for k in range(m_lp)])
    dt = solver.dt
    e_min, e_max = float(wv.min()), float(wv.max())

    print(f"  band = [{e_d_true:.6f}, {beta_self:.6f}) electronic, width={beta_self - e_d_true:.6f} Ha, "
          f"spectrum width = {e_max - e_min:.4f} Ha")
    for n_grid in N_GRIDS:
        bound = _build_moment_bound_lp(moments, e_min, e_max, dt, e_d_true, beta_self,
                                        n_grid=n_grid)
        flag = "n/a" if bound is None else ("YES" if bound <= 0.05 else "no")
        print(f"    n_grid={n_grid:>7}  eta_bound={bound:.6f}  <=0.05? {flag}")


def main():
    for c in CASES:
        run_case(c)


if __name__ == "__main__":
    main()
