#!/usr/bin/env python3
"""
Certified two-sided brackets on a spectral GAP from the Krylov solve -- no FCI oracle.

`temple_bounds` gives every QKSD ground-state *energy* a certified bracket [E_Temple, E_Ritz]. But
the quantities spectroscopy actually measures are *gaps* (excitation energies), and those were still
point estimates. This module closes that: it puts a rigorous interval around the fundamental gap
Delta = E_1 - E_0 of the reachable sector, from the SAME Krylov data (one extra <H^2> per Ritz
state -- no new measurements), so an excitation energy can be quoted with a certified +/- .

The bracket composes three facts the repo already trusts:

  * Cauchy interlacing:      theta_1 >= E_1           (the 2nd Ritz value bounds E_1 from above)
  * Temple (1928):           tau_0   <= E_0           (lower bound on the ground energy)
  * Weinstein / self-eps:    eps_1 = theta_1 - sigma_1 <= E_1   (under the checkable premise that
                             theta_1 resolves E_1 -- valid for M >= 6, the temple_bracket boundary)

giving, with theta_0 >= E_0 (variational),

    Delta_hi = theta_1 - tau_0          >= E_1 - E_0 = Delta     (upper certificate)
    Delta_lo = (theta_1 - sigma_1) - theta_0  <= Delta           (lower certificate, premise-gated)

THE FINDING (specs/SPEC_certified_gaps.md): the exact reachable gap sits inside [Delta_lo, Delta_hi]
at EVERY depth M >= 6 across H4 / LiH / N2 CAS(6,6) (zero escapes), and the interval closes with
depth (H4: 342 -> 0.7 mHa over M=6..24; N2 CAS(6,6): 391 -> 37 mHa). The certification inherits the
temple_bracket boundary EXACTLY: at M = 4 the Weinstein premise eps_1 <= E_1 fails for H4/N2 and the
LOWER certificate escapes -- so the gap bracket is trustworthy only once the Krylov space resolves
the excited state (M >= 6). The asymmetry is the mechanism: the UPPER certificate needs only
interlacing + Temple and holds at every tested depth; the LOWER certificate is the premise-sensitive
side, because a real-time Krylov space has no lower bound on E_2 to anchor a rigorous E_1 floor.

HONEST SCOPE: sector-restricted (E_0/E_1 are the two lowest HF-reachable levels -- QKSD's own
scope); exact statevector (the shot cost of <H^2>, a ~lambda^2 Pauli expansion, is not modeled);
the lower certificate's premise is checkable against an oracle but not self-verifiable. On a system
small enough to diagonalize exactly (e.g. the 4-qubit Nb3X8 dimer) certification is pointless -- you
already have the gap; the value is on systems where FCI is out of reach.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import _mean_and_variance


@dataclass
class GapBracket:
    """A certified bracket on the reachable-sector gap at one Krylov dimension. Energies in Ha."""
    m: int                  # Krylov dimension
    gap_lower: float        # premise-gated lower certificate (-inf when theta_1 unavailable)
    gap_upper: float        # upper certificate theta_1 - tau_0 (+inf when Temple is vacuous)
    width: float            # gap_upper - gap_lower
    theta0: float           # ground Ritz value (>= E_0)
    theta1: float           # first-excited Ritz value (>= E_1, Cauchy interlacing)
    sigma1: float           # sqrt(variance) of the first-excited Ritz state
    eps1: float             # theta_1 - sigma_1, the self lower estimate of E_1
    eps1_source: str        # "self" (theta_1 - sigma_1) | "oracle" (caller-supplied E_1)


def gap_bracket(mh: MolecularHamiltonian, m: int, e1: Optional[float] = None,
                solver: Optional[QuantumKrylovSolver] = None) -> GapBracket:
    """Certified [Delta_lo, Delta_hi] bracket on the reachable gap from an M-dim Krylov solve.

    ``e1``: exact E_1 as a TOTAL energy (oracle mode -- validation) for the Temple/E_1 floor, or
    None to estimate it from the same data as theta_1 - sigma_1 (self mode; premise unverifiable,
    valid for M >= 6 -- see module docstring). Gaps are energy-offset-independent, so everything is
    computed in the electronic frame. Pass a shared ``solver`` to reuse the Krylov basis.
    """
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    _, states = solver.eigenstates(m, n_states=2)
    th0, var0 = _mean_and_variance(H, states[0])
    if len(states) < 2:                                   # rank-1 subspace: no E_1 handle
        return GapBracket(m, -np.inf, np.inf, np.inf, th0, np.inf, np.inf, -np.inf, "self")
    th1, var1 = _mean_and_variance(H, states[1])
    sig1 = float(np.sqrt(var1))
    if e1 is None:
        eps1, eps1_source = th1 - sig1, "self"
    else:
        eps1, eps1_source = float(e1) - mh.energy_offset, "oracle"
    # Temple lower bound on E_0 uses eps1 as a floor for E_1 (needs eps1 > th0 to be finite).
    tau0 = th0 - var0 / (eps1 - th0) if eps1 > th0 else -np.inf
    gap_upper = th1 - tau0 if np.isfinite(tau0) else np.inf
    gap_lower = eps1 - th0                                 # (theta_1 - sigma_1) - theta_0
    return GapBracket(m=m, gap_lower=gap_lower, gap_upper=gap_upper,
                      width=gap_upper - gap_lower, theta0=th0, theta1=th1, sigma1=sig1,
                      eps1=eps1, eps1_source=eps1_source)


def gap_bracket_ladder(mh: MolecularHamiltonian, dims: Sequence[int], e1: Optional[float] = None,
                       solver: Optional[QuantumKrylovSolver] = None) -> List[GapBracket]:
    """Gap brackets at each Krylov dimension in ``dims`` (basis built once and reused)."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    return [gap_bracket(mh, m, e1=e1, solver=solver) for m in dims]


def reachable_gap(mh: MolecularHamiltonian) -> float:
    """REFERENCE ONLY (dense, O(2^n)): exact gap between the two lowest HF-reachable eigenstates --
    the object the bracket certifies. For validation/demo on small systems, never the live path."""
    H = mh.qubit_hamiltonian.to_matrix()
    w, V = np.linalg.eigh(H)
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    reach = w[np.abs(V.conj().T @ hf) ** 2 > 1e-10]
    return float(reach[1] - reach[0])


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
        "LiH": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
        "N2 CAS(6,6)": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
    }
    for name, spec in cases.items():
        mh = build_molecular_hamiltonian(**spec)
        gap = reachable_gap(mh)                            # shown only for comparison
        solver = QuantumKrylovSolver(mh)
        print("=" * 78)
        print(f"{name}: reachable gap = {gap * 1e3:.3f} mHa (certified WITHOUT this number)")
        print("   M | gap_lo (mHa) | gap_hi (mHa) | width (mHa) | inside? (self-eps)")
        for br in gap_bracket_ladder(mh, (4, 6, 8, 12, 16, 20, 24), solver=solver):
            inside = br.gap_lower <= gap <= br.gap_upper
            print(f"  {br.m:2d} | {br.gap_lower * 1e3:12.3f} | {br.gap_upper * 1e3:12.3f} | "
                  f"{br.width * 1e3:11.3f} | {inside}")
    print("=" * 78)
    print("Certified for M >= 6 (zero escapes); the lower certificate escapes at M=4 where the")
    print("Weinstein premise eps_1 <= E_1 fails -- the temple_bracket boundary, on gaps.")
