#!/usr/bin/env python3
"""
Certified two-sided energy brackets from the Krylov solve -- Temple/Weinstein lower bounds on the
QKSD ground Ritz eigenstate.

Every energy the validated pipeline produces is a variational UPPER bound (QKSD Ritz values, PDS,
DMRG); nothing certifies from below. Temple's inequality (Proc. R. Soc. A 119, 276, 1928) closes
the gap: for any state with mean th = <H> and variance var = <H^2> - <H>^2, and any eps <= E_1
with th < eps,

    E_0  >=  th - var / (eps - th)          (Temple)
    E_0  >=  th - sqrt(var)                 (Weinstein, when th lies closer to E_0 than E_1)

Applied to the Krylov eigenstate |Psi_0(M)> (whose variance -> 0 as M grows), each solve carries a
certified bracket [E_Temple, E_Ritz] containing the exact reachable-sector ground energy, at the
cost of ONE extra expectation value <Psi_0|H^2|Psi_0> (a sparse matvec here). Modern chemistry
usage: Pollak & Martinazzo, JCTC 15, 1498 (2019); JCP 152, 244110 (2020).

HONEST SCOPE (see specs/SPEC_temple_bracket.md): certification is sector-restricted (E_0/E_1 are
the lowest reachable levels -- the same scope as QKSD itself). The Temple premise eps <= E_1 is
rigorous only with an oracle gap; the oracle-free mode eps = theta_1 - sigma_1 cannot verify its
own premise and is valid here only for M >= 6 (gated boundary). Exact statevector: the hardware
shot cost of <H^2> (a ~lambda^2-sized Pauli expansion) is not modeled. A lower bound of -inf
(eps <= theta_0) is valid but vacuous -- check ``width`` is finite before quoting.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver


@dataclass
class EnergyBracket:
    """A certified bracket at one Krylov dimension. All energies are TOTAL energies in Ha."""
    m: int                  # Krylov dimension
    upper: float            # Ritz value (variational upper bound)
    lower: float            # Temple lower bound (-inf when eps <= theta_0: valid but vacuous)
    weinstein_lower: float  # looser variance-only lower bound
    width: float            # upper - lower (inf when the Temple bound is vacuous)
    variance: float         # <H^2> - <H>^2 of the ground Ritz state (electronic frame, Ha^2)
    eps: float              # gap input fed to Temple (total energy)
    eps_source: str         # "oracle" (caller-supplied E_1) | "self" (theta_1 - sigma_1)


def mean_and_variance(H, psi):
    """(<H>, <H^2> - <H>^2) for a normalized state, from ONE matvec -- the measurement every
    certificate in this repo is built on (Temple, gap, dipole, and their noise variants).

    The variance is clipped at zero: it is mathematically non-negative, but near convergence
    <H^2> - <H>^2 is a difference of two nearly equal O(1) numbers and rounds negative, which
    would make sqrt(var) NaN and silently poison a bound downstream.
    """
    hpsi = H @ psi
    th = float((psi.conj() @ hpsi).real)
    var = float((hpsi.conj() @ hpsi).real - th ** 2)
    return th, max(var, 0.0)


def krylov_bracket(mh: MolecularHamiltonian, m: int, eps: Optional[float] = None,
                   solver: Optional[QuantumKrylovSolver] = None) -> EnergyBracket:
    """Certified [Temple, Ritz] bracket from an M-dimensional Krylov solve.

    ``eps``: exact E_1 as a TOTAL energy (oracle mode -- validation), or None to estimate it from
    the same Krylov data as theta_1 - sigma_1 (self mode; premise unverifiable, see module
    docstring). Pass a shared ``solver`` to reuse the cached Krylov basis across calls.
    """
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    offset = mh.energy_offset
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    energies, states = solver.eigenstates(m, n_states=2)
    th0, var0 = mean_and_variance(H, states[0])
    if eps is None:
        if len(energies) > 1:
            th1, var1 = mean_and_variance(H, states[1])
            eps_e = th1 - np.sqrt(var1)
        else:
            eps_e = -np.inf                              # rank-1 subspace: no E_1 estimate
        eps_source = "self"
    else:
        eps_e = float(eps) - offset
        eps_source = "oracle"
    lower_e = th0 - var0 / (eps_e - th0) if eps_e > th0 else -np.inf
    upper = th0 + offset
    lower = lower_e + offset if np.isfinite(lower_e) else -np.inf
    return EnergyBracket(m=m, upper=upper, lower=lower,
                         weinstein_lower=th0 - np.sqrt(var0) + offset,
                         width=upper - lower, variance=var0,
                         eps=eps_e + offset if np.isfinite(eps_e) else -np.inf,
                         eps_source=eps_source)


def bracket_ladder(mh: MolecularHamiltonian, dims: Sequence[int], eps: Optional[float] = None,
                   solver: Optional[QuantumKrylovSolver] = None) -> List[EnergyBracket]:
    """Brackets at each Krylov dimension in ``dims`` (basis built once and reused)."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    return [krylov_bracket(mh, m, eps=eps, solver=solver) for m in dims]


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
        "N2 CAS(6,6)": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
    }
    for name, spec in cases.items():
        mh = build_molecular_hamiltonian(**spec)
        fci = mh.ground_state_energy()
        solver = QuantumKrylovSolver(mh)
        print("=" * 78)
        print(f"{name}: FCI={fci:.6f} Ha (shown only for comparison -- the bracket is computed "
              "without it)")
        print("   M |  Ritz err (mHa) | self-Temple bracket width (mHa) | FCI inside?")
        for br in bracket_ladder(mh, (4, 6, 8, 12, 16, 20, 24), solver=solver):
            inside = br.lower <= fci <= br.upper
            print(f"  {br.m:2d} | {(br.upper - fci) * 1e3:14.4f} | {br.width * 1e3:30.4f} | "
                  f"{inside}   ({br.eps_source})")
    print("=" * 78)
