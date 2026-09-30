#!/usr/bin/env python3
"""
Map where ODMD's E_1 estimate repairs SPEC_temple_bracket's oracle-free premise, and where it
silently breaks it (chem-z3h; backlog: specs/BACKLOG.md, "Method rungs", "A non-variational
estimator cannot supply the Temple premise").

`SPEC_temple_bracket`'s self-consistent mode feeds `eps = theta_1 - sigma_1` (from the SAME
Krylov data as the ground state) into Temple's inequality as an estimate of E_1; G4(b) records
that this premise (eps <= E_1) FAILS at M=4 on H4/N2. `odmd.py` measures the identical physical
object -- the survival amplitude s_k = <phi_0|U^k|phi_0>, which is exactly the first row of the
solver's own S matrix (``solver._subspace_matrices(K)[1][0, :K]``, K >= the Krylov depth already
built) -- and its second DMD eigenphase is an independent, no-extra-measurement estimate of E_1.
Substituting it for theta_1 - sigma_1 is not free: ODMD is explicitly non-variational (no
eps <= E_1 guarantee at any depth), so the substitution can repair the Temple premise at one K
and violate it -- or worse, produce a genuine bracket CONTAINMENT ESCAPE (lower > E_0, a false
certificate) -- at another.

THE MEASURED FINDING (this script; also pinned in tests/test_odmd_temple_repair_map.py):
- G-a: at M=4, the self-mode premise is violated on H4 and N2 (matches SPEC_temple_bracket G4b)
  but NOT on H2 -- H2's 4-qubit reachable sector is exhausted by M=4 (theta_1 = E_1 exactly to
  float precision), so the "repair" question is vacuous there. The repair is non-vacuous on
  2 of the 3 ODMD reference systems (SPEC_odmd.md: H2, H4 chain, N2 CAS(6,6)), not all 3.
- G-b: the real danger is NOT at deep K, it is at the SHALLOWEST usable K. At K=4 (the minimum
  `_dmd_modes` accepts; only a rank-2 fit -- 2 DMD modes resolved), the ODMD E_1
  estimate overshoots the true E_1 by ~1-1.5 Ha (H4: +950 mHa, N2: +1530 mHa) and this is not
  merely a premise violation: it is a genuine CONTAINMENT ESCAPE -- the resulting Temple lower
  bound exceeds the true E_0 by up to 1.26 mHa (N2, M=2) -- reproduced at M in {2,4,6} (H4) and
  M in {2,4,6,8} (N2). By K=6 (rank >= 3) the escape disappears and does not recur through
  K=28 at M in {2,4}, even though the *premise* (eps <= E_1) flips invalid again at several K in
  that range (H4: K=6-14, K=28; N2: K=6-10, K=18-28) without ever re-escaping containment at
  these M. So "premise violated" and "bracket actually wrong" are NOT the same event here --
  the escape is concentrated at the single most under-resolved depth, not spread across "deep K"
  as an earlier informal backlog scout probe suggested (whose specific eps/E_1 numbers do not
  reproduce against the current H4 geometry in tests/test_temple_bracket_spec.py; this script
  supersedes those numbers with a systematic K-sweep against the oracle E_1).
- K=4 is already below every depth SPEC_odmd.md validates (its gates start at K=8); this is a
  known-out-of-scope depth for ODMD, not a new failure mode of ODMD itself -- the finding is
  that `temple_bounds.krylov_bracket` has no internal check that would catch the escape if K=4
  ODMD output were fed to it as if it were an oracle (SPEC_odmd_uq's recorded blind spot: a
  single signal's own resampling cannot see this kind of model bias -- see specs/SPEC_odmd_uq.md
  G4). The escape is caught here only because the probe has independent oracle access to E_0/E_1
  (dense diagonalisation) -- production code does not.

HONEST SCOPE: noiseless/exact statevector (matches SPEC_temple_bracket); reachable-sector E0/E1
by dense diagonalisation (SPEC_qksd_excited.md convention); 3 systems (H2, H4 chain, N2 CAS(6,6)
-- SPEC_odmd.md's own reference set, not SPEC_temple_bracket's LiH-inclusive set, since LiH is
not one of ODMD's validated systems); base Krylov depths M in {2,4,6,8}; ODMD signal depths K in
{4,5,6,7,8,10,...,28}. This is a measured map, not a certificate: it does not claim K=4 is the
ONLY dangerous depth on every system/geometry, only that it is dangerous on these three.
"""
from __future__ import annotations

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from odmd import odmd_spectrum
from temple_bounds import krylov_bracket

SYSTEMS = {
    "h2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "h4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "n2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
M_DEPTHS = (2, 4, 6, 8)
K_DEPTHS = (4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28)


def reachable_e0_e1(mh) -> tuple[float, float]:
    """REFERENCE ONLY (dense O(2^n)): exact ground/first-excited reachable-sector energies."""
    w_eig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    psi0 = np.asarray(mh.hf_state().data, dtype=complex)
    pops = np.abs(V.conj().T @ psi0) ** 2
    reach = w_eig[pops > 1e-8].real + mh.energy_offset
    return float(reach[0]), float(reach[1])


def odmd_eps1(mh, solver: QuantumKrylovSolver, k: int) -> float | None:
    """ODMD's E_1 estimate (total energy) from the first row of the solver's OWN S matrix at
    depth k -- no extra measurement beyond what the Krylov solve already built. None if the
    noiseless DMD fit does not resolve a 2nd mode (rank < 2, e.g. k=5 on H4 here)."""
    solver._ensure_basis(k)
    _, S = solver._subspace_matrices(k)
    energies, _, _ = odmd_spectrum(S[0, :k], solver.dt)
    if len(energies) < 2:
        return None
    return float(energies[1] + mh.energy_offset)


def k_vs_validity_map(key: str) -> None:
    mh = build_molecular_hamiltonian(**SYSTEMS[key])
    e0, e1 = reachable_e0_e1(mh)
    print("=" * 100)
    print(f"{key}: E0={e0:.6f} E1={e1:.6f} Ha  (reachable gap {(e1 - e0) * 1e3:.3f} mHa)")
    for m in M_DEPTHS:
        solver = QuantumKrylovSolver(mh)
        br_self = krylov_bracket(mh, m, solver=solver)
        print(f"  self-mode (theta1-sigma1) @ M={m:2d}: eps={br_self.eps:.6f}  "
              f"premise_violated={br_self.eps > e1}")
        header = (f"    {'K':>3} | {'eps_odmd (Ha)':>14} | {'vs E1 (mHa)':>12} | "
                  f"{'premise_ok':>10} | {'escape (mHa)':>13} | {'CONTAINMENT':>11}")
        print(header)
        for k in K_DEPTHS:
            eps = odmd_eps1(mh, solver, k)
            if eps is None:
                print(f"    {k:3d} | {'--':>14} | {'--':>12} | {'--':>10} | {'--':>13} | "
                      f"{'no E1 fit':>11}")
                continue
            br = krylov_bracket(mh, m, eps=eps, solver=solver)
            escape_mha = (br.lower - e0) * 1e3
            valid = br.lower <= e0 + 1e-9
            premise_ok = eps <= e1
            print(f"    {k:3d} | {eps:14.6f} | {(eps - e1) * 1e3:+12.4f} | "
                  f"{str(premise_ok):>10} | {escape_mha:+13.4f} | "
                  f"{'VALID' if valid else 'INVALID':>11}")


if __name__ == "__main__":
    for key in SYSTEMS:
        k_vs_validity_map(key)
    print("=" * 100)
