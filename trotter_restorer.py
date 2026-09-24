#!/usr/bin/env python3
"""
Signal-domain Richardson for Trotter ODMD -- the leak that isn't there, and the window where
extrapolating the SIGNAL beats extrapolating the ENERGY.

specs/SPEC_trotter_floor_restorer.md v1.0 set out to plug a "Trotter variational-floor leak": the
claim that Trotterization shifts the effective spectrum so a Krylov/ODMD solve returns E < E_FCI,
destroying certified bounds. Two things came out of measuring it.

  * THERE IS NO LEAK. Across 4 systems x reps {1,2,4} x 8 window lengths -- 96 configurations --
    the order-2 Suzuki Trotter ODMD energy is ABOVE E_FCI every single time, minimum excess
    +0.26 mHa, typical +1 to +20 mHa. The effective-Hamiltonian bias here is POSITIVE, so the
    premise the spec was built to rescue does not occur (G1).
  * SIGNAL-domain Richardson is not uniformly better than the SCALAR (energy-domain) Richardson
    that `trotter_odmd.richardson_energy` already ships -- it is better in a narrow window and
    then stops working entirely. Extrapolating s_k linearly only cancels the Trotter bias while
    the bias's accumulated PHASE is small, because the bias enters s_k as exp(-i k tau dE). Past
    K*tau*dE ~ 0.15 rad the extrapolant degenerates onto the un-extrapolated fine signal and buys
    nothing. Measured on stretched H2: at K = 6 (0.12 rad) signal Richardson leaves -0.0045 mHa,
    10x better than scalar's -0.0476; at K >= 8 (>= 0.16 rad) it returns exactly the raw reps=4
    value, +1.0267 mHa, while scalar stays flat at -0.0476 for every K.

The refinement parameter is `reps`, NOT dt. v1.0's formula pairs s_k(dt) with s_k(2dt), which sit
at DIFFERENT physical times k*dt and 2k*dt, so the O(dt^2) term cannot cancel. Refining reps at
fixed tau keeps every sample at the same physical time k*tau while halving dt_eff = tau/reps --
the construction `trotter_odmd.build_trotter_odmd_problem(mh, n, reps=r)` already provides.

HONEST SCOPE: statevector-simulated Trotter circuits (no device noise, no native-gate
transpilation), inherited from `trotter_odmd`. Energies are compared to FCI in the centered frame.
The phase-window law is calibrated on these systems, not derived. Nothing here computes a Lehmann
bracket: ODMD works from the scalar signal alone, so it yields eigenphases, not Ritz STATES, and
<u|H^2|u> is not available from s_k -- v1.0's step 4 is out of scope for this path and is recorded
as such rather than faked.
"""
from __future__ import annotations

import numpy as np

from odmd import odmd_energy
from trotter_odmd import build_trotter_odmd_problem, richardson_energy

_PROBLEMS: dict[tuple, object] = {}


def trotter_problems(mh, n: int = 32, reps=(1, 2, 4)) -> dict:
    """Trotter ODMD problems at each rep count, built once per (mh id, n) and cached."""
    key = (id(mh), n, tuple(reps))
    if key not in _PROBLEMS:
        _PROBLEMS[key] = {r: build_trotter_odmd_problem(mh, n=n, reps=r) for r in reps}
    return _PROBLEMS[key]


def signal_richardson(s_coarse, s_fine, step_ratio: float = 2.0, order: int = 2) -> np.ndarray:
    """Richardson-extrapolate two time-ALIGNED Trotter signals to dt_eff -> 0.

    ``s_coarse`` and ``s_fine`` must be sampled at the SAME physical times (refine ``reps`` at
    fixed tau, never dt), with dt_eff differing by ``step_ratio``. Returns
    (w*s_fine - s_coarse)/(w-1), w = step_ratio**order -- the signal-domain analogue of
    :func:`trotter_odmd.richardson_energy`. Valid only while the accumulated phase error is small;
    see :func:`phase_budget`.
    """
    w = float(step_ratio) ** order
    return (w * np.asarray(s_fine) - np.asarray(s_coarse)) / (w - 1.0)


def phase_budget(K: int, tau: float, bias: float) -> float:
    """Accumulated Trotter phase error K*tau*|bias| (radians) over the ODMD window.

    The diagnostic that decides whether :func:`signal_richardson` can work at all: the bias enters
    the signal as exp(-i k tau bias), so linear extrapolation of s_k cancels it only while this
    stays well below 1. Measured boundary ~0.15 rad (gated in tests/test_trotter_restorer_spec.py).
    """
    return float(K * tau * abs(bias))


def restored_energy(problems: dict, K: int, mode: str = "signal_extrap",
                    pair: tuple[int, int] = (2, 4)) -> float:
    """Ground energy (centered frame) from Trotter signals under one mitigation mode.

    ``raw``           -- un-mitigated ODMD on the fine signal.
    ``signal_extrap`` -- Richardson on the complex signal BEFORE diagonalization (v1.0's proposal).
    ``scalar_extrap`` -- Richardson on the two solved energies (the shipped `trotter_odmd` path).
    """
    c, f = pair
    if mode == "raw":
        return odmd_energy(problems[f].s[:K], problems[f].tau)[0]
    if mode == "signal_extrap":
        sx = signal_richardson(problems[c].s[:K], problems[f].s[:K], step_ratio=f / c)
        return odmd_energy(sx, problems[c].tau)[0]
    if mode == "scalar_extrap":
        ec = odmd_energy(problems[c].s[:K], problems[c].tau)[0]
        ef = odmd_energy(problems[f].s[:K], problems[f].tau)[0]
        return richardson_energy(ec, ef, step_ratio=f / c)
    raise ValueError(f"unknown mode {mode!r}")


if __name__ == "__main__":
    import csv

    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    systems = {"H2@0.74": "H 0 0 0; H 0 0 0.74", "H2@1.6": "H 0 0 0; H 0 0 1.6",
               "H2@2.5": "H 0 0 0; H 0 0 2.5",
               "H4@0.9": "H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"}
    rows = []
    for name, atom in systems.items():
        mh = build_molecular_hamiltonian(atom=atom)
        P = trotter_problems(mh, n=32)
        ref, tau = P[1].ref, P[1].tau
        bias2 = odmd_energy(P[2].s, tau)[0] - ref
        print(f"\n=== {name} ===  (bias at reps=2: {bias2 * 1e3:+.4f} mHa)")
        for K in (6, 8, 12, 16, 24, 32):
            e = {m: restored_energy(P, K, m) for m in ("raw", "signal_extrap", "scalar_extrap")}
            pb = phase_budget(K, tau, bias2)
            rows.append(dict(system=name, K=K, phase_budget_rad=round(pb, 4),
                             **{m: round((v - ref) * 1e3, 5) for m, v in e.items()},
                             raw_leak=bool(e["raw"] < ref),
                             extrap_undershoot=bool(min(e["signal_extrap"],
                                                        e["scalar_extrap"]) < ref)))
            print(f"  K={K:2d} phase={pb:6.3f} rad  raw={(e['raw'] - ref) * 1e3:+9.5f}  "
                  f"signal={(e['signal_extrap'] - ref) * 1e3:+9.5f}  "
                  f"scalar={(e['scalar_extrap'] - ref) * 1e3:+9.5f} mHa")
    with open("data/trotter_restorer_bench.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote data/trotter_restorer_bench.csv  raw Trotter below FCI in "
          f"{sum(r['raw_leak'] for r in rows)}/{len(rows)} configurations; the EXTRAPOLANT "
          f"undershoots in {sum(r['extrap_undershoot'] for r in rows)}")
