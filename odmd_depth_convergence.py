#!/usr/bin/env python3
"""
A depth-convergence check on ODMD's eps(K) -- the guard SPEC_odmd_temple_repair_map.md named and
deferred (chem-z3h -> chem-36m).

That spec's measured (M, K) map showed ODMD's E_1 estimate is NOT monotone in K: on H4 the Temple
premise is repaired at K=16-26 and violated again at K=28 (+0.74 mHa); on N2 it is repaired at
K=12-16 and violated again at K=18-28 (+14.8...+35.4 mHa). "The premise currently holds at this K"
is therefore not a safe reason to trust eps(K) -- the next K can flip it back, with no oracle to
warn you. This module adds the cheap, oracle-free guard that spec's SS7 deferred: require eps(K) to
have STOPPED MOVING over a trailing window of already-measured K (no look-ahead -- only
information available if you were actually running this live) before treating it as a candidate
Temple floor.

HONEST SCOPE: NECESSARY, NOT SUFFICIENT (specs/SPEC_odmd_depth_convergence.md G3). Convergence
means the estimate is stable, not that it is correct -- N2's eps(K) sits +34...+35 mHa above the
true E_1 while stable within ~1 mHa over K=24-28, a consistently-biased fixed point. This is the
same model-misspecification blind spot SPEC_odmd_uq G4 recorded for single-signal resampling: no
amount of watching one signal settle down can see that it settled on the wrong answer. Pairs with,
does not replace, ``temple_bounds.EnergyBracket.premise_refuted`` /
``certified_gaps.GapBracket.premise_refuted`` and an eventual rigorous (Lehmann/Kato) E_1 floor.

Window/tol (3 grid points, 2 mHa) were chosen by inspecting the already-gated
SPEC_odmd_temple_repair_map.md SS10 numbers -- re-derived here against that committed record, not
re-measured blind.
"""
from __future__ import annotations

from typing import Mapping, Optional, Sequence

from odmd_temple_repair_map import K_DEPTHS, odmd_eps1

TOL = 1e-9   # Ha: SPEC_temple_bracket G1 / SPEC_odmd_temple_repair_map's containment tolerance


def depth_converged(odmd: Mapping[int, Optional[float]], k_grid: Sequence[int], k: int,
                    window: int = 3, tol: float = 2e-3) -> bool:
    """True iff the ``window`` deepest grid points at or below ``k`` are all resolved (not
    ``None``) and span no more than ``tol`` Ha. No look-ahead: only K's <= k are used, matching
    what is actually knowable while running an ODMD sweep live (watch the estimate until it stops
    moving, then stop -- never using a future, unmeasured K to judge the present one)."""
    ks = sorted(kk for kk in k_grid if kk <= k)
    if len(ks) < window:
        return False
    chosen = ks[-window:]
    vals = [odmd.get(kk) for kk in chosen]
    if any(v is None for v in vals):
        return False
    return (max(vals) - min(vals)) <= tol


def converged_depths(odmd: Mapping[int, Optional[float]], k_grid: Sequence[int],
                     window: int = 3, tol: float = 2e-3) -> list[int]:
    """Every K in ``k_grid`` at which :func:`depth_converged` holds."""
    return [k for k in k_grid if depth_converged(odmd, k_grid, k, window=window, tol=tol)]


def safe_eps1(mh, solver, k: int, k_grid: Sequence[int] = K_DEPTHS, window: int = 3,
             tol: float = 2e-3) -> Optional[float]:
    """ODMD's E_1 estimate (total energy) at depth ``k``, gated by depth-convergence -- ``None``
    (abstain) unless the trailing window has stabilized. Prefer this over calling
    ``odmd_temple_repair_map.odmd_eps1`` directly when the result will be substituted as the E_1
    floor into ``temple_bounds.krylov_bracket(eps=...)`` or ``certified_gaps.gap_bracket(e1=...)``
    -- both accept an externally supplied floor, so one guard covers either consumer. Computes
    eps(K') for every K' in ``k_grid`` up to ``k`` (the cost the convergence check itself needs);
    pass a shared ``solver`` to reuse the Krylov basis with other calls on the same system.

    A ``True`` result is NOT a correctness proof (SPEC_odmd_depth_convergence.md G3) -- pair the
    returned value with the oracle-free ``premise_refuted`` flag on whichever bracket consumes it.
    """
    odmd = {kk: odmd_eps1(mh, solver, kk) for kk in k_grid if kk <= k}
    if not depth_converged(odmd, [kk for kk in k_grid if kk <= k], k, window=window, tol=tol):
        return None
    return odmd[k]


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
    from odmd_temple_repair_map import SYSTEMS

    for key in ("h4", "n2", "lih_cas", "lih"):
        mh = build_molecular_hamiltonian(**SYSTEMS[key])
        solver = QuantumKrylovSolver(mh)
        odmd = {k: odmd_eps1(mh, solver, k) for k in K_DEPTHS}
        conv = converged_depths(odmd, K_DEPTHS)
        print(f"{key}: converged K = {conv}")
        for k in conv:
            print(f"   K={k:3d}  eps(K) = {odmd[k]:.6f} Ha")
