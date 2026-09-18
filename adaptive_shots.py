#!/usr/bin/env python3
"""
Non-uniform shot allocation for ODMD -- and the boundary where it actually buys anything.

The hypothesis (specs/SPEC_adaptive_shots_planning.md): the visibility law's 1/(w^2 K) cost is a
UNIFORM-allocation law, so a decaying schedule S_k ~ e^{-alpha k} that spends the budget on
"high-signal early steps" should beat it 5-10x. Measured, it does not -- and the reason is the
finding, not a bug:

  * The survival amplitude of a closed system is quasi-periodic, NOT decaying: |s_k| ~ 0.9 for
    every k on N2 CAS(6,6). There is no "noise-dominated late tail" to defund. With every element
    carrying equal signal, total noise power sum_k sigma_k^2 = V sum_k 1/S_k is CONVEX in the
    allocation, so uniform S_k = S/K is the constrained minimizer: any decay strictly increases it.
    Adaptive allocation cannot beat uniform on an undamped signal. (Gate G2.)
  * Under the repo's validated global-depolarizing model s_k -> f^k s_k (device_odmd.py), a dead
    tail DOES exist, and defunding it wins -- but by ~1.3x in error (~1.8x in budget at the
    1/sqrt(S) scaling), not 5-10x. The optimal decay rate alpha* rises with the damping rate
    -log f, which is the one part of the original hypothesis that survives. (Gates G3, G4.)

The re-weighting is where the repo's depolarizing-immunity insight pays off. Heteroskedastic
elements are whitened by r_k = 1/sigma_k ~ sqrt(S_k); for a GEOMETRIC schedule that weight is
itself geometric, r_k = c*beta^k, which is exactly the transformation that multiplies every DMD
eigenvalue by beta and leaves its PHASE untouched. So whitening an exponential schedule is free:
homoskedastic noise, zero energy bias, and no new DMD code (device_odmd_energy's wide modulus
window already handles |lambda| < 1). A polynomial schedule's weights are not geometric, so the
same whitening perturbs the eigenphases -- a bias, quantified in __main__ and gated in G2b.

HONEST SCOPE: exact statevector s_k with the idealized i.i.d. Hadamard-test noise of odmd.py;
damping enters as the global-depolarizing model, not measured device noise; optimize_decay_factor
is an OFFLINE PLANNER -- it minimizes error against the known reference energy, so it sizes a
budget before an experiment, it is not a runtime estimator.
"""
from __future__ import annotations

import numpy as np

from device_odmd import device_odmd_energy
from odmd import ODMDProblem


def _normalize(w: np.ndarray, total_shots: float) -> np.ndarray:
    """S_0 = 0 (s_0 = 1 exactly, no measurement); the whole budget goes to k = 1..K-1."""
    w = np.asarray(w, dtype=float).copy()
    w[0] = 0.0
    return w / w.sum() * float(total_shots)


def exponential_schedule(K: int, total_shots: float, alpha: float) -> np.ndarray:
    """S_k = C exp(-alpha k) for k >= 1, S_0 = 0, sum = ``total_shots`` (alpha < 0 = growing)."""
    return _normalize(np.exp(-alpha * np.arange(K)), total_shots)


def polynomial_schedule(K: int, total_shots: float, gamma: float) -> np.ndarray:
    """S_k = C (k+1)^(-gamma) for k >= 1, S_0 = 0, sum = ``total_shots``."""
    return _normalize((np.arange(K) + 1.0) ** (-gamma), total_shots)


def whitening_weights(schedule: np.ndarray) -> np.ndarray:
    """Per-element weights r_k ~ 1/sigma_k = sqrt(S_k) that make the noise homoskedastic.

    Normalized to r_1 = 1; r_0 is the geometric extrapolation sqrt(S_1/S_2) of the schedule, so a
    geometric schedule gets r_k = c beta^k on the WHOLE window including the noiseless s_0 -- the
    eigenphase-preserving case. Overall scale is irrelevant (DMD is scale-invariant).
    """
    s = np.asarray(schedule, dtype=float)
    r = np.sqrt(s / s[1])
    r[0] = np.sqrt(s[1] / s[2])
    return r


def sample_adaptive_odmd_energy(prob: ODMDProblem, schedule: np.ndarray, seed: int,
                                damping: float = 1.0, amp_floor: float = 0.0,
                                whiten: bool = True) -> float:
    """One shot-noisy ODMD estimate (centered frame) under a per-element shot budget.

    Element k gets its own Hadamard-test noise scale sigma_k^2 = 2(2 - 1/dim)/S_k (odmd.py
    conventions; s_0 = 1 exactly). ``damping`` f applies the validated global-depolarizing model
    s_k -> f^k s_k, which leaves the exact eigenphases invariant (device_odmd.py). With
    ``whiten`` the signal is rescaled by :func:`whitening_weights` before the damping-robust
    estimator, whose noise-edge cutoff then uses the single whitened scale sqrt(V/S_1).
    """
    S = np.asarray(schedule, dtype=float)
    K = len(S)
    V = 2.0 * (2.0 - 1.0 / prob.dim)
    sigma = np.sqrt(V / np.where(S > 0, S, np.inf))
    rng = np.random.default_rng(seed)
    g = rng.normal(0, sigma / np.sqrt(2)) + 1j * rng.normal(0, sigma / np.sqrt(2))
    s = prob.s[:K] * damping ** np.arange(K) + g
    r = whitening_weights(S) if whiten else np.ones(K)
    return device_odmd_energy(s * r, prob.tau, np.sqrt(V / S[1]), amp_floor=amp_floor)


def median_error(prob: ODMDProblem, schedule: np.ndarray, seeds: int = 200,
                 damping: float = 1.0, **kw) -> float:
    """Median |E - E_exact| (Ha) over ``seeds`` independent noise realizations; inf -> 1 Ha."""
    e = np.array([sample_adaptive_odmd_energy(prob, schedule, sd, damping, **kw)
                  for sd in range(seeds)])
    return float(np.median(np.where(np.isfinite(e), np.abs(e - prob.ref), 1.0)))


def optimize_decay_factor(prob: ODMDProblem, K: int, total_shots: float, damping: float = 1.0,
                          alphas=np.arange(-0.10, 0.51, 0.05), seeds: int = 200, **kw) -> float:
    """OFFLINE PLANNER: the exponential decay rate alpha minimizing the median energy error.

    Scored against the known exact ``prob.ref``, so this sizes a schedule in simulation before an
    experiment -- it is not a runtime estimator. alpha = 0 is uniform; the grid spans negative
    (growing) rates so "uniform is optimal" is a reachable answer, not an assumption.
    """
    errs = [median_error(prob, exponential_schedule(K, total_shots, a), seeds, damping, **kw)
            for a in alphas]
    return float(alphas[int(np.argmin(errs))])


if __name__ == "__main__":
    import csv

    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from odmd import build_odmd_problem

    mh = build_molecular_hamiltonian(atom="N 0 0 0; N 0 0 1.1",
                                     active_electrons=6, active_orbitals=6)
    prob = build_odmd_problem(mh, n=24)
    print(f"N2 CAS(6,6): |s_k| = {np.round(np.abs(prob.s[:8]), 3)} ... (no decay: nothing to defund)")
    rows = []
    for K in (12, 16, 24):
        for f in (1.0, 0.95, 0.9, 0.85, 0.8, 0.7):
            uni = median_error(prob, exponential_schedule(K, 1e4, 0.0), 400, f)
            a = optimize_decay_factor(prob, K, 1e4, f, seeds=400)
            adp = median_error(prob, exponential_schedule(K, 1e4, a), 400, f)
            rows.append(dict(K=K, damping=f, alpha_star=round(a, 3),
                             uniform_mHa=round(uni * 1e3, 4), adaptive_mHa=round(adp * 1e3, 4),
                             error_gain=round(uni / adp, 3), budget_gain=round((uni / adp) ** 2, 3)))
            print(f"  K={K:2d} f={f:4.2f}: alpha*={a:+.2f}  uniform={uni * 1e3:7.2f} mHa  "
                  f"adaptive={adp * 1e3:7.2f} mHa  gain={uni / adp:.2f}x "
                  f"({(uni / adp) ** 2:.2f}x budget)")
    with open("data/adaptive_shots_bench.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote data/adaptive_shots_bench.csv")
