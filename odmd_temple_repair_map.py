#!/usr/bin/env python3
"""
Map where ODMD's E_1 estimate repairs SPEC_temple_bracket's oracle-free premise, and where it
silently breaks it (chem-z3h; specs/SPEC_odmd_temple_repair_map.md, gated by
tests/test_odmd_temple_repair_map_spec.py).

`SPEC_temple_bracket`'s self mode feeds eps = theta_1 - sigma_1 (from the SAME Krylov data as the
ground state) into Temple's inequality as a stand-in for E_1; its G4(b) records that this premise
(eps <= E_1) FAILS at M=4. `odmd.py` reads the survival amplitude s_k = <phi_0|U^k|phi_0>, which
is exactly the first row of the solver's own S matrix, and its second DMD eigenphase is another
E_1 estimate. That estimate is free (no extra measurement) only for K <= M: S is Toeplitz,
S_ij = s_{j-i}, so an M-dim solve has measured s_0..s_{M-1} and nothing more. ODMD is
non-variational, so the substitution can repair the premise at one K and break it -- or produce a
genuine CONTAINMENT ESCAPE (lower > E_0, a false certificate) -- at another.

Per system and per (M, K) this records: eps_ODMD(K) vs E_1 (the premise), whether the Temple
bracket still contains E_0 (an ORACLE check -- the only real INVALID verdict), and whether the
oracle-free `EnergyBracket.premise_refuted` (eps > theta_1(M) >= E_1) can tell. The measured map
and verdicts are in the spec, section 10.

HONEST SCOPE: noiseless statevector; reachable-sector E_0/E_1 by dense diagonalisation (validation
only); ODMD runs in the solver's UNCENTERED frame with `solver.dt`, not SPEC_odmd's validated
centered frame, and K < 8 is below SPEC_odmd's validated depths; a finite (M, K) grid on five
STO-3G system variants -- a measured map, not a certificate.
"""
from __future__ import annotations

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from odmd import odmd_spectrum
from reachability import reachable_eigenpairs
from temple_bounds import krylov_bracket

LIH = "Li 0 0 0; H 0 0 1.6"
SYSTEMS = {
    "h2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "h4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "lih_cas": dict(atom=LIH, active_electrons=2, active_orbitals=5),   # certified_gaps' LiH
    "lih": dict(atom=LIH),                                              # temple_bracket's LiH
    "n2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
M_DEPTHS = (2, 4, 6, 8)
K_DEPTHS = (4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28)
TOL = 1e-9   # Ha: SPEC_temple_bracket G1's containment tolerance, for every comparison here


def odmd_eps1(mh, solver: QuantumKrylovSolver, k: int) -> float | None:
    """ODMD's E_1 estimate (total energy) from the first row of the solver's OWN S matrix at
    depth k (extra overlap measurements only when k exceeds the solve's own depth). None if the
    noiseless DMD fit does not resolve a 2nd mode."""
    _, S = solver._subspace_matrices(k)
    energies, _, _ = odmd_spectrum(S[0, :k], solver.dt)
    if len(energies) < 2:
        return None
    return float(energies[1] + mh.energy_offset)


def k_vs_validity_map(key: str) -> dict:
    """One system's raw map: oracle ``e0``/``e1`` (total Ha) and ``n_reachable``; ``self``
    {M: self-mode EnergyBracket}; ``odmd`` {K: eps_ODMD or None}; ``cells`` {(M, K): EnergyBracket
    with eps = eps_ODMD(K)}. One shared solver, so every cell reuses the same Krylov vectors."""
    mh = build_molecular_hamiltonian(**SYSTEMS[key])
    w, _ = reachable_eigenpairs(mh, tol=1e-8)
    solver = QuantumKrylovSolver(mh)
    odmd = {k: odmd_eps1(mh, solver, k) for k in K_DEPTHS}
    return dict(e0=float(w[0]) + mh.energy_offset, e1=float(w[1]) + mh.energy_offset,
                n_reachable=len(w),
                self={m: krylov_bracket(mh, m, solver=solver) for m in M_DEPTHS},
                odmd=odmd,
                cells={(m, k): krylov_bracket(mh, m, eps=eps, solver=solver)
                       for m in M_DEPTHS for k, eps in odmd.items() if eps is not None})


def summarize(mp: dict) -> dict:
    """The recorded form of one system's map (spec section 10); cells are (M, K), energies mHa."""
    e0, e1, cells, self4 = mp["e0"], mp["e1"], mp["cells"], mp["self"][4]
    over = {c for c, br in cells.items() if br.eps > e1 + TOL}
    escapes = {c for c, br in cells.items() if br.lower > e0 + TOL}
    return dict(
        self_m4_mha=(self4.eps - e1) * 1e3,                 # > 0: self-mode premise violated
        self_m4_contains=bool(self4.lower <= e0 + TOL),
        overshoot_k=sorted({k for _, k in over}),
        escapes=sorted(escapes),
        worst_mha=max((br.lower - e0) * 1e3 for br in cells.values()),
        degraded=sorted(over - escapes),
        flagged=sorted(c for c, br in cells.items() if br.premise_refuted),
        repair_k_m4=sorted(k for (m, k), br in cells.items()
                           if m == 4 and self4.eps > e1 + TOL and br.eps <= e1 + TOL
                           and np.isfinite(br.lower)),
    )


if __name__ == "__main__":
    for key in SYSTEMS:
        mp = k_vs_validity_map(key)
        e0, e1 = mp["e0"], mp["e1"]
        print("=" * 96)
        print(f"{key}: E0={e0:.6f} E1={e1:.6f} Ha (reachable gap {(e1 - e0) * 1e3:.3f} mHa, "
              f"{mp['n_reachable']} reachable levels)")
        for m, sb in mp["self"].items():
            print(f"  M={m}: self eps-E1 {(sb.eps - e1) * 1e3:+10.4f} mHa, lower-E0 "
                  f"{(sb.lower - e0) * 1e3:+10.4f} mHa, theta1-E1 {(sb.theta1 - e1) * 1e3:+10.4f}")
            print(f"    {'K':>3} | {'eps-E1 (mHa)':>13} | {'lower-E0 (mHa)':>14} | "
                  f"{'eps-theta1 (mHa)':>16} | verdict")
            for k, eps in mp["odmd"].items():
                br = mp["cells"].get((m, k))
                if br is None:
                    print(f"    {k:3d} | {'no E1 fit':>13} |")
                    continue
                verdict = ("INVALID" if br.lower > e0 + TOL else
                           "degraded" if eps > e1 + TOL else "ok")
                print(f"    {k:3d} | {(eps - e1) * 1e3:+13.4f} | {(br.lower - e0) * 1e3:+14.4f} | "
                      f"{(eps - br.theta1) * 1e3:+16.4f} | {verdict}"
                      f"{' FLAGGED' if br.premise_refuted else ''}{' free' if k <= m else ''}")
        for name, val in summarize(mp).items():
            print(f"  {name}: {val}")
    print("=" * 96)
