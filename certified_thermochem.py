#!/usr/bin/env python3
"""
Certified error bars on a RELATIVE energy (reaction / dissociation / stretch) -- no FCI oracle.

The certified arc bounds absolute energies (`temple_bounds`), gaps (`certified_gaps`), and
properties (`certified_dipole`). But chemistry's currency is *relative* energies -- reaction
energies, barriers, dissociation energies -- and those reach experiment most directly. A certified
absolute energy at each of two geometries composes into a certified difference:

    E_A in [tau_A, rho_A]  (Temple lower, Ritz upper),  E_B in [tau_B, rho_B]
    =>  Delta = E_B - E_A  in  [tau_B - rho_A,  rho_B - tau_A]

a rigorous interval on the relative energy, from Krylov data at each endpoint (no FCI).

THE FINDING (specs/SPEC_certified_thermochem.md): on the H4 symmetric stretch (equilibrium ->
stretched) the exact relative energy (8.2255 eV) sits inside the certified interval at every depth
(zero escapes), and the interval is dominated ENTIRELY by the strongly-correlated stretched
endpoint -- the equilibrium bracket closes ~25x faster (width 0.001 vs 0.025 eV at M=6), so the
certified error bar LOCALIZES where the correlation difficulty lives. And the certificate inherits
the temple-bracket premise at that endpoint: at intermediate depth the stretched Temple lower bound
is vacuous, so Delta carries a one-sided (upper-only) certificate there before closing two-sided at
larger M -- the same premise regime `certified_gaps`/`gap_selfcheck` chart.

HONEST SCOPE: sector-restricted ground states at each geometry (QKSD's scope); the lower bound rests
on the temple premise (valid at sufficient depth; the harder endpoint sets it); exact statevector
(the <H^2> shot cost is not modeled); the certificate contains the *in-basis FCI* relative energy --
basis-set error vs experiment is a separate, uncertified quantity. On systems small enough to
diagonalize the certificate is a demonstration, not a need.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import krylov_bracket

HARTREE_TO_EV = 27.211386245988


@dataclass
class RelEnergyBracket:
    """Certified interval on Delta = E_B - E_A (total energies, Ha). Endpoint widths exposed so the
    caller can see which geometry dominates the uncertainty."""
    m: int
    delta: float            # point estimate: Ritz_B - Ritz_A
    delta_lower: float      # tau_B - rho_A  (-inf if the harder endpoint's Temple is vacuous)
    delta_upper: float      # rho_B - tau_A  (+inf if the other endpoint's Temple is vacuous)
    width: float            # delta_upper - delta_lower (inf if either side vacuous)
    width_a: float          # certified bracket width at endpoint A (inf if vacuous)
    width_b: float          # certified bracket width at endpoint B


def certified_relative_energy(mh_a: MolecularHamiltonian, mh_b: MolecularHamiltonian, m: int,
                              solver_a: Optional[QuantumKrylovSolver] = None,
                              solver_b: Optional[QuantumKrylovSolver] = None,
                              e1_a: Optional[float] = None,
                              e1_b: Optional[float] = None) -> RelEnergyBracket:
    """Certified [Delta_lo, Delta_hi] on Delta = E(B) - E(A) from an M-dim Krylov solve at each
    endpoint. Pass exact E_1's (oracle) or None (self mode; premise at sufficient depth)."""
    ba = krylov_bracket(mh_a, m, eps=e1_a, solver=solver_a)
    bb = krylov_bracket(mh_b, m, eps=e1_b, solver=solver_b)
    delta = bb.upper - ba.upper
    delta_lower = bb.lower - ba.upper           # tau_B - rho_A
    delta_upper = bb.upper - ba.lower           # rho_B - tau_A
    return RelEnergyBracket(m=m, delta=delta, delta_lower=delta_lower, delta_upper=delta_upper,
                            width=delta_upper - delta_lower, width_a=ba.width, width_b=bb.width)


def certified_relative_energy_ladder(mh_a: MolecularHamiltonian, mh_b: MolecularHamiltonian,
                                     dims: Sequence[int]) -> List[RelEnergyBracket]:
    """Certified relative-energy bracket at each Krylov dimension (one solver per endpoint)."""
    sa, sb = QuantumKrylovSolver(mh_a), QuantumKrylovSolver(mh_b)
    return [certified_relative_energy(mh_a, mh_b, m, solver_a=sa, solver_b=sb) for m in dims]


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    def chain(r):
        return f"H 0 0 0; H 0 0 {r}; H 0 0 {2 * r}; H 0 0 {3 * r}"

    mh_eq = build_molecular_hamiltonian(atom=chain(0.9))
    mh_st = build_molecular_hamiltonian(atom=chain(2.3))
    de_fci = (mh_st.ground_state_energy() - mh_eq.ground_state_energy()) * HARTREE_TO_EV
    print("=" * 76)
    print(f"H4 symmetric stretch (0.9 -> 2.3 A): certified relative energy (FCI {de_fci:.4f} eV, "
          "shown only for comparison)")
    print("   M | certified Delta (eV)        | width | endpoint widths eq / stretch (eV) | in?")
    for rb in certified_relative_energy_ladder(mh_eq, mh_st, (6, 8, 10, 12, 16, 20)):
        lo, hi = rb.delta_lower * HARTREE_TO_EV, rb.delta_upper * HARTREE_TO_EV
        we = rb.width_a * HARTREE_TO_EV if rb.width_a != float("inf") else float("inf")
        ws = rb.width_b * HARTREE_TO_EV if rb.width_b != float("inf") else float("inf")
        inside = lo - 1e-6 <= de_fci <= hi + 1e-6
        span = f"[{lo:8.3f}, {hi:8.3f}]"
        w = f"{(hi - lo):6.3f}" if rb.width != float("inf") else "  inf "
        print(f"  {rb.m:2d} | {span:>26} | {w} | {we:8.3f} / {ws:8.3f}            | {inside}")
    print("=" * 76)
    print("Certified relative energy = Ritz_B - Ritz_A bracketed by the two Temple/Ritz brackets;")
    print("the strongly-correlated (stretched) endpoint dominates the error bar.")
