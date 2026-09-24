#!/usr/bin/env python3
"""
Interval-dominance screening -- pruning a candidate library by CERTIFIED bracket dominance,
and the two conditions that decide whether it saves anything at all.

The idea (specs/SPEC_interval_dominance_screening.md): solve every candidate cheaply first, and
whenever a candidate's certified LOWER bound exceeds the best certified UPPER bound anywhere in
the library, it provably cannot be the global minimum -- drop it before paying for convergence.
The dominance test itself is sound arithmetic. What v1.0 got wrong is that the SAVING is not a
property of the method:

  * Pruning power is set by ONE ratio -- the library's energy spread over the bracket width at the
    cheap dimension. Measured on H4 chains at M=2: spread/width = 4.45 prunes 6/8 immediately,
    spread/width = 1.00 prunes 0/8. Same loop, same code, same budget. A percentage quoted without
    the library that produced it means nothing, which is why `sweep_metrics` returns the ratio
    alongside the saving.
  * The bound has to be a real bound. Temple's SELF mode (eps = theta_1 - sigma_1, the only mode
    available when you have not already solved the candidate) is NOT rigorous -- measured here at
    up to 0.695 mHa ABOVE the true energy on H4, which is licence to prune the winner. Oracle mode
    (eps = true E_1) never violates. This is the same kill SPEC_excited_state_certification found
    one rung up, and screening inherits it: the cheap certificate a screening loop actually has is
    the unsound one.

On a HOMOGENEOUS library the unsoundness is largely common-mode -- every candidate's bound is
inflated together, so the winner's margin stays negative and no false prune occurs (gated, with
the measured margin, in G1b). That is a measured property of such libraries, NOT a guarantee, and
a heterogeneous library whose members converge at different rates has no such protection.

HONEST SCOPE: exact statevector Krylov (no shot noise, no device noise). "QPU cost" is a COUNT of
Krylov basis vectors -- one controlled time-evolution family each -- not a wall-clock or gate
count; the accounting assumes the basis is extended incrementally so a candidate stopped at M
costs M. Ground-state energies only. STO-3G H2 is excluded by construction: its Krylov space
saturates at M=2, so "cheap vs converged" does not exist there (the same category error
SPEC_excited_state_certification's G2 records).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import krylov_bracket


@dataclass
class Candidate:
    """One library member and its running certified state. ``cost`` counts Krylov basis vectors."""
    name: str
    mh: MolecularHamiltonian
    eps: float | None = None          # oracle E_1 (total Ha); None -> Temple self mode
    active: bool = True
    current_m: int = 0
    lower: float = -np.inf
    upper: float = np.inf
    cost: int = 0
    pruned_at: int | None = None
    _solver: QuantumKrylovSolver | None = field(default=None, repr=False)

    def refine(self, m: int) -> None:
        """Advance this candidate to Krylov dimension ``m`` and update its certified bracket."""
        if self._solver is None:
            self._solver = QuantumKrylovSolver(self.mh)
        b = krylov_bracket(self.mh, m, eps=self.eps, solver=self._solver)
        self.current_m, self.lower, self.upper, self.cost = m, b.lower, b.upper, m


def _dominance_prune(cands: list[Candidate]) -> list[Candidate]:
    """Deactivate every active candidate whose certified lower bound exceeds the best upper bound.

    ``L_j > U_best`` proves j cannot be the global minimum. A -inf lower bound (Temple vacuous)
    never prunes, which is the correct conservative behaviour.
    """
    active = [c for c in cands if c.active]
    if len(active) <= 1:
        return active
    u_best = min(c.upper for c in active)
    for c in active:
        if np.isfinite(c.lower) and c.lower > u_best:
            c.active, c.pruned_at = False, c.current_m
    return [c for c in cands if c.active]


def interval_dominance_sweep(library: list[Candidate], max_m: int = 12, m_step: int = 2,
                             start_m: int = 2) -> tuple[Candidate, dict]:
    """Adaptive screen: refine every survivor, prune by dominance, repeat until one or max_m."""
    pool_curve = []
    for m in range(start_m, max_m + 1, m_step):
        for c in library:
            if c.active:
                c.refine(m)
        survivors = _dominance_prune(library)
        pool_curve.append((m, len(survivors)))
        if len(survivors) <= 1:
            break
    winner = min((c for c in library if c.active), key=lambda c: c.upper)
    return winner, sweep_metrics(library, max_m, pool_curve)


def brute_force_sweep(library: list[Candidate], target_m: int = 12) -> tuple[Candidate, dict]:
    """Reference: take every candidate to ``target_m``, no pruning."""
    for c in library:
        c.refine(target_m)
    winner = min(library, key=lambda c: c.upper)
    return winner, sweep_metrics(library, target_m, [(target_m, len(library))])


def sweep_metrics(library: list[Candidate], max_m: int, pool_curve) -> dict:
    """Cost, saving, and the spread/width ratio that actually governs the saving."""
    spent = sum(c.cost for c in library)
    brute = len(library) * max_m
    widths = np.array([c.upper - c.lower for c in library])
    finite = widths[np.isfinite(widths)]
    uppers = np.array([c.upper for c in library])
    spread = float(uppers.max() - uppers.min())
    med_w = float(np.median(finite)) if finite.size else np.inf
    return dict(basis_vectors=spent, brute_force_basis_vectors=brute,
                saving=1.0 - spent / brute, pool_curve=pool_curve,
                spread_Ha=spread, median_width_Ha=med_w,
                spread_over_width=spread / med_w if med_w > 0 else np.inf,
                pruned=sum(1 for c in library if not c.active))


def h4_library(spacings, oracle: bool = False) -> list[Candidate]:
    """Linear H4 chains at the given spacings -- the gated library.

    STO-3G H2 is deliberately not used: its Krylov space saturates at M = 2 (bracket width 0.00
    mHa), so every candidate is exactly converged at the cheapest dimension and the loop measures
    nothing. H4 is the smallest chain here that is still unconverged at M = 2.
    ``oracle`` supplies the exact E_1 to Temple, the only mode in which the bound is rigorous.
    """
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    out = []
    for a in spacings:
        mh = build_molecular_hamiltonian(
            atom=f"H 0 0 0; H 0 0 {a:.5f}; H 0 0 {2 * a:.5f}; H 0 0 {3 * a:.5f}")
        eps = None
        if oracle:
            w = np.linalg.eigvalsh(mh.qubit_hamiltonian.to_matrix())
            eps = float(w[1]) + mh.energy_offset
        out.append(Candidate(name=f"a={a:.3f}", mh=mh, eps=eps))
    return out


if __name__ == "__main__":
    import csv

    rows = []
    for tag, spacings in (("wide", np.linspace(0.8, 2.2, 8)),
                          ("tight", np.linspace(1.00, 1.07, 8))):
        for mode in ("oracle", "self"):
            lib = h4_library(spacings, oracle=(mode == "oracle"))
            winner, m = interval_dominance_sweep(lib, max_m=12)
            ref = h4_library(spacings, oracle=(mode == "oracle"))
            bwin, bm = brute_force_sweep(ref, target_m=12)
            rows.append(dict(library=tag, mode=mode, winner=winner.name,
                             brute_winner=bwin.name, agree=winner.name == bwin.name,
                             basis_vectors=m["basis_vectors"],
                             brute_basis_vectors=m["brute_force_basis_vectors"],
                             saving=round(m["saving"], 4),
                             spread_mHa=round(m["spread_Ha"] * 1e3, 3),
                             median_width_mHa=round(m["median_width_Ha"] * 1e3, 3),
                             spread_over_width=round(m["spread_over_width"], 3),
                             pool_curve=";".join(f"{k}:{v}" for k, v in m["pool_curve"])))
            print(f"  {tag:5s}/{mode:6s}: winner={winner.name} (brute {bwin.name}) "
                  f"cost {m['basis_vectors']}/{m['brute_force_basis_vectors']} "
                  f"saving={m['saving'] * 100:5.1f}%  spread/width={m['spread_over_width']:.2f}")
    with open("data/screening_loop_performance.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote data/screening_loop_performance.csv")
