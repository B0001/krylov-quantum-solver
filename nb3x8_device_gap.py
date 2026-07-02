#!/usr/bin/env python3
"""
Nb3X8 cluster charge gaps through the simulated-hardware pipeline -- the capstone composition.

The repo's two validated threads meet: the Nb3X8 downfolded dimer clusters (nb3x8_gaps.py, exact
ED gaps from arXiv:2501.10320's cRPA parameters) and the device-validated ODMD stack
(SPEC_odmd -> SPEC_trotter_odmd -> SPEC_device_odmd). The charge gap

    Delta = E(N+1) + E(N-1) - 2 E(N),      N = 2 (half filling)

needs only GROUND-state energies in three particle sectors; real-time evolution conserves
particle number, so each sector's HF reference pins its sector, and each energy is one
depolarizing-immune ODMD run on genuinely-Trotterized Hadamard-test circuits. The Trotter bias
does NOT cancel between sectors (reps=1 leaves 12% on Nb3I8 -- the recorded finding; the
near-commuting Nb3F8 sits at 0.1 meV), so Richardson extrapolation across reps is what turns the
circuit data into the exact gap.

HONEST SCOPE (see specs/SPEC_nb3x8_device_gap.md): a study/composition, not a method rung. This
is the ISOLATED-CLUSTER gap (842.44 meV for Nb3I8 LT-bulk) -- an upper bound on the solid gap,
which band broadening moves to ~600-650 meV (SPEC_nb3x8_gaps' corrected conclusion). Simulated
device (Aer depolarizing; no coherent errors/crosstalk/drift); density-density interactions
only; the dimer's charged sectors are nearly free -- a pipeline demonstration at validation
scale, not a hard correlated-electron benchmark. Energies in meV throughout.
"""
from __future__ import annotations

from typing import Dict

from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals, fixed_filling_energy
from nb3x8_gaps import dimer_cluster_integrals
from odmd import build_odmd_problem, odmd_energy
from device_odmd import centered_frame, device_odmd_energy, measure_survival_signal
from trotter_odmd import build_trotter_odmd_problem, richardson_energy

# Charge-gap combination: E(3) + E(1) - 2 E(2).
SECTOR_WEIGHTS = ((1, 1.0), (2, -2.0), (3, 1.0))


def sector_models(U0: float, t: float, Us: float) -> Dict[int, ModelIntegrals]:
    """The N = 1, 2, 3 electron sectors of the generalized Hubbard dimer cluster."""
    base = dimer_cluster_integrals(U0, t, Us)
    return {1: ModelIntegrals(base.h1, base.eri, 0.0, (1, 0), 2),
            2: base,
            3: ModelIntegrals(base.h1, base.eri, 0.0, (2, 1), 2)}


def exact_gap(U0: float, t: float, Us: float) -> float:
    """Sector-FCI charge gap (the reference; pinned against SPEC_nb3x8_gaps' recorded values)."""
    models = sector_models(U0, t, Us)
    return sum(w * fixed_filling_energy(models[n]) for n, w in SECTOR_WEIGHTS)


def statevector_gap(U0: float, t: float, Us: float, n: int = 16) -> float:
    """Charge gap from exact-evolution ODMD signals (one run per sector)."""
    models = sector_models(U0, t, Us)
    gap = 0.0
    for nelec, w in SECTOR_WEIGHTS:
        prob = build_odmd_problem(models[nelec].to_hamiltonian(), n=n)
        gap += w * (odmd_energy(prob.s, prob.tau)[0] + prob.offset)
    return gap


def circuit_gap(U0: float, t: float, Us: float, reps: int, n: int = 16) -> float:
    """Charge gap from the EXACT ground eigenphases of the Trotterized step circuits.

    The Trotter bias of this gap does not cancel between sectors; it follows the order-2 law in
    ``reps`` and is removed by :func:`trotter_odmd.richardson_energy` across a reps pair.
    """
    models = sector_models(U0, t, Us)
    gap = 0.0
    for nelec, w in SECTOR_WEIGHTS:
        tp = build_trotter_odmd_problem(models[nelec].to_hamiltonian(), n=n, reps=reps)
        gap += w * (tp.e_circuit + tp.offset)
    return gap


def device_gap(U0: float, t: float, Us: float, shots: int, noise_model, seed: int,
               trotter_reps: int = 1, n: int = 8) -> float:
    """Charge gap measured through the Aer hardware stack (Hadamard tests, device noise).

    One damping-robust ODMD run per sector; the per-sector Aer seed is
    ``seed*100 + trotter_reps*10 + nelec`` (distinct, reproducible draws). The returned gap still
    carries the Trotter bias of ``trotter_reps`` -- see :func:`device_gap_richardson`.
    """
    models = sector_models(U0, t, Us)
    sigma = 1.0 / float(shots) ** 0.5
    gap = 0.0
    for nelec, w in SECTOR_WEIGHTS:
        mh = models[nelec].to_hamiltonian()
        _, tau, mu = centered_frame(mh)
        s = measure_survival_signal(mh, n, shots=shots, noise_model=noise_model,
                                    seed=seed * 100 + trotter_reps * 10 + nelec,
                                    trotter_reps=trotter_reps)
        gap += w * (device_odmd_energy(s, tau, sigma) + mh.energy_offset + mu)
    return gap


def device_gap_richardson(U0: float, t: float, Us: float, shots: int, noise_model,
                          seed: int, n: int = 8) -> float:
    """Richardson-extrapolated device gap over the reps (1, 2) pair.

    Extrapolation is linear, so it commutes with the gap combination. It pays only while the
    Trotter bias exceeds the device-noise floor (SPEC_trotter_odmd R1 -- gated crossover:
    cx = 1e-4 yes, cx = 3e-4 no, at 32768 shots).
    """
    g1 = device_gap(U0, t, Us, shots, noise_model, seed, trotter_reps=1, n=n)
    g2 = device_gap(U0, t, Us, shots, noise_model, seed, trotter_reps=2, n=n)
    return richardson_energy(g1, g2)


if __name__ == "__main__":
    from hybrid_quantum_solver.noise import build_depolarizing_noise_model
    from nb3x8_gaps import NB3X8_LT_BULK

    print("Nb3X8 LT-bulk dimer-cluster charge gaps [meV] -- the pipeline ladder")
    print("(isolated-cluster gaps: upper bounds on the solid; see SPEC_nb3x8_gaps)")
    print(f"{'material':>8} | {'FCI':>8} | {'statevec-ODMD':>13} | {'circuit reps=1':>14} | "
          f"{'Richardson(1,2)':>15}")
    for name, p in NB3X8_LT_BULK.items():
        ref = exact_gap(**p)
        sv = statevector_gap(**p)
        c1, c2 = circuit_gap(**p, reps=1), circuit_gap(**p, reps=2)
        print(f"{name:>8} | {ref:8.2f} | {sv:13.2f} | {c1:14.2f} | "
              f"{richardson_energy(c1, c2):15.2f}")
    p = NB3X8_LT_BULK["Nb3I8"]
    nm = build_depolarizing_noise_model(1e-5, 1e-4, 1e-4)
    g = device_gap_richardson(**p, shots=32768, noise_model=nm, seed=0)
    print(f"\nNb3I8 through the noisy device (cx=1e-4, 32768 shots, seed 0): "
          f"{g:.1f} meV vs FCI {exact_gap(**p):.2f}")
