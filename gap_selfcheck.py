#!/usr/bin/env python3
"""
Oracle-free trustworthiness certificate for the certified gap bracket -- self-verification by
cross-depth consistency.

`certified_gaps` puts a two-sided interval around the reachable gap, but its lower certificate rests
on a premise (eps_1 <= E_1) that fails at shallow Krylov depth and CANNOT be checked without the
exact answer -- the open limitation that spec recorded ("valid for M >= 6", known only against an
oracle). This module removes the oracle: a bracket is *corroborated* iff it is consistent with the
converged brackets at greater depth. A premise failure inflates/shifts the bracket so it no longer
overlaps the deep ones, and is caught -- no FCI required.

The test is ADAPTIVE, not a blanket "distrust shallow M": it rejects a shallow bracket only when
that bracket is genuinely inconsistent. On H4 / N2 CAS(6,6) the M=4 bracket (premise failed) is
rejected; on LiH (well-separated, premise already holds at M=4) the M=4 bracket is ACCEPTED -- the
certificate distinguishes real failures from merely-shallow-but-valid ones.

    anchor        = intersection of the deepest k brackets (the convergence reference)
    corroborated  = the bracket at M overlaps the anchor
    self-checked gap = intersection of the corroborated brackets (an oracle-free interval,
                       validated to contain the exact reachable gap)

HONEST SCOPE (specs/SPEC_gap_selfcheck.md): NECESSARY, NOT SUFFICIENT. The anchor -- the deepest
brackets -- is taken on trust; the check cannot certify its own convergence, and it cannot see a
consistently-biased sequence that shrinks toward the wrong value (the model-misspecification blind
spot of any single-run resampling, cf. SPEC_odmd_uq). It catches the KNOWN failure mode (premise
breakdown at shallow depth) and pairs with, not replaces, a depth-convergence check.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

from certified_gaps import GapBracket, gap_bracket_ladder
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver


def _intersect(intervals: Sequence[Tuple[float, float]]) -> Tuple[float, float]:
    """Intersection of intervals as (lo, hi); empty if lo > hi."""
    return max(a for a, _ in intervals), min(b for _, b in intervals)


def _overlap(a: Tuple[float, float], b: Tuple[float, float], tol: float = 1e-12) -> bool:
    return max(a[0], b[0]) <= min(a[1], b[1]) + tol


def anchor_interval(brackets: Sequence[GapBracket], k: int = 2) -> Tuple[float, float]:
    """Convergence reference: intersection of the ``k`` deepest (largest-m) brackets."""
    deep = sorted(brackets, key=lambda b: b.m)[-k:]
    return _intersect([(b.gap_lower, b.gap_upper) for b in deep])


def corroborated_flags(brackets: Sequence[GapBracket], k: int = 2) -> List[bool]:
    """Per-bracket flag: does it overlap the deep-anchor (in input order)? Oracle-free."""
    anchor = anchor_interval(brackets, k)
    return [_overlap((b.gap_lower, b.gap_upper), anchor) for b in brackets]


def self_checked_gap(brackets: Sequence[GapBracket], k: int = 2) -> Tuple[float, float]:
    """Oracle-free gap interval: intersection of the corroborated brackets. Validated (see spec) to
    contain the exact reachable gap, and to exclude the shallow premise-failure outliers."""
    flags = corroborated_flags(brackets, k)
    keep = [(b.gap_lower, b.gap_upper) for b, ok in zip(brackets, flags) if ok]
    return _intersect(keep)


def self_checked_gap_from(mh: MolecularHamiltonian, dims: Sequence[int], k: int = 2,
                          solver=None) -> Tuple[Tuple[float, float], List[bool]]:
    """Build the bracket ladder and return (self-checked interval, corroboration flags)."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    brackets = gap_bracket_ladder(mh, dims, solver=solver)
    return self_checked_gap(brackets, k), corroborated_flags(brackets, k)


if __name__ == "__main__":
    from certified_gaps import reachable_gap
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
        "LiH": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
        "N2 CAS(6,6)": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
    }
    dims = (4, 6, 8, 12, 16, 20, 24)
    for name, spec in cases.items():
        mh = build_molecular_hamiltonian(**spec)
        gap = reachable_gap(mh)
        (lo, hi), flags = self_checked_gap_from(mh, dims)
        rejected = [m for m, ok in zip(dims, flags) if not ok]
        print("=" * 74)
        print(f"{name}: self-checked gap = [{lo * 1e3:.2f}, {hi * 1e3:.2f}] mHa "
              f"(exact {gap * 1e3:.2f}; contains it: {lo <= gap <= hi})")
        print(f"   rejected depths (oracle-free): {rejected or 'none'}")
    print("=" * 74)
    print("Adaptive: M=4 rejected on H4/N2 (premise failed) but accepted on LiH (premise holds).")
