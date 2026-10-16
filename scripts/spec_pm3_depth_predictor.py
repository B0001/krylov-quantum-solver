#!/usr/bin/env python3
"""
chem-j71 driver (specs/SPEC_pm3_depth_predictor.md): a follow-up to chem-pm3
(scripts/spec_pm3_subspace_eta_bound.py), which found that the self-mode floor's overshoot
beta - E_d does NOT order its three escape witnesses by how hard the eta-bound is to satisfy.

Part A: does a "weighted gap" WG = min(gap_lo, gap_hi) / eta -- a band-edge gap divided by the HF
weight it has to exclude -- track the bound's difficulty (b_min, m_pass) better than the overshoot,
across a wider witness set than chem-pm3's three? New escape witnesses are found by the same
guard-passing/invalid-floor search SPEC_subspace_floor_resolvability.md used, over a pre-registered
grid of linear/asym H6 geometries, picking the SMALLEST M at which the guard passes and the floor
overshoots the true E_d (a deterministic rule fixed before any correlation number exists).

Part B: measures (does not just analytically bound) how much the bound moves under realistic shot
noise on the moments S_0k, compared against chem-pm3's analytic worst-case slope and a "typical
case" L2 estimate.

Run: uv run python scripts/spec_pm3_depth_predictor.py
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from pyscf import lib
from scipy.stats import spearmanr

from hf_overlap_subspace import _weinstein_intervals_disjoint
from hybrid_quantum_solver.certified_overlap import residual_norm
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.noise import shot_noise_sigma
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from reachability import _dense_hf_projection
from scripts.spec_pm3_subspace_eta_bound import (
    _asym_chain,
    _build_moment_bound_lp,
    _chain,
    WITNESSES as PM3_WITNESSES,
    run_witness,
)

_REACHABLE_TOL = 1e-8
M_SEARCH = (6, 8, 10, 12, 16, 20, 24)

# Pre-registered candidate geometries for the new-witness search (chosen before any escape/no-escape
# or moment-bound outcome was known): a linear H6 R-sweep avoiding chem-pm3's own R=1.0/1.1, plus two
# asym variants bracketing chem-pm3's own [0.9,0.9,2.2,0.9,0.9] defect size.
CANDIDATES = [
    dict(name="linear H6 R=0.85", atom=_chain(6, 0.85), d=3),
    dict(name="linear H6 R=0.90", atom=_chain(6, 0.90), d=3),
    dict(name="linear H6 R=0.95", atom=_chain(6, 0.95), d=3),
    dict(name="linear H6 R=1.05", atom=_chain(6, 1.05), d=3),
    dict(name="linear H6 R=1.15", atom=_chain(6, 1.15), d=3),
    dict(name="linear H6 R=1.20", atom=_chain(6, 1.20), d=3),
    dict(name="linear H6 R=1.25", atom=_chain(6, 1.25), d=3),
    dict(name="asym H6 [0.9,0.9,2.0,0.9,0.9]", atom=_asym_chain([0.9, 0.9, 2.0, 0.9, 0.9]), d=3),
    dict(name="asym H6 [0.9,0.9,2.5,0.9,0.9]", atom=_asym_chain([0.9, 0.9, 2.5, 0.9, 0.9]), d=3),
]


def find_escape_witnesses(candidates=CANDIDATES, m_search=M_SEARCH):
    """Smallest M in m_search where the disjoint guard passes AND the self-mode floor overshoots
    the true E_d by > 1e-6 Ha (SPEC_subspace_floor_resolvability's escape condition), per candidate
    geometry. Deterministic and pre-registered: a candidate that never escapes contributes nothing.

    pyscf.lib.num_threads(1) is pinned around the SCF/integral build (bd memory chem-7cx-root-cause:
    the free-RHF branch is non-bit-associative across OpenMP thread counts on near-degenerate
    geometries -- pinning makes this reproducible regardless of the ambient thread count).
    """
    found = []
    for c in candidates:
        t0 = time.time()
        ambient = lib.num_threads()
        lib.num_threads(1)
        try:
            mh = build_molecular_hamiltonian(atom=c["atom"])
            wv, _, pops = _dense_hf_projection(mh)
        finally:
            lib.num_threads(ambient)
        solver = QuantumKrylovSolver(mh)
        Hs = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
        offset = mh.energy_offset
        d = c["d"]
        reach = np.where(pops > _REACHABLE_TOL)[0]
        reach_e, reach_w = wv[reach], pops[reach]
        e_d = float(reach_e[d])

        picked = None
        for m in m_search:
            energies, states = solver.eigenstates(m, n_states=d + 1)
            centers = [energies[k] - offset for k in range(d + 1)]
            sigmas = [residual_norm(Hs, states[k], centers[k]) for k in range(d + 1)]
            disjoint = _weinstein_intervals_disjoint(centers, sigmas)
            beta = centers[d] - sigmas[d]
            overshoot = beta - e_d
            if disjoint and overshoot > 1e-6:
                picked = (m, beta, overshoot)
                break
        wall = time.time() - t0
        if picked is None:
            print(f"{c['name']}: no escape in M in {m_search} ({wall:.1f}s)")
            continue
        m, beta, overshoot = picked
        gap_lo = e_d - float(reach_e[d - 1])
        above_idx = np.where(reach_e >= beta)[0][0]
        gap_hi = float(reach_e[above_idx]) - beta
        in_missed = (reach_e >= e_d) & (reach_e < beta)
        eta = float(reach_w[in_missed].sum())
        print(f"{c['name']}: ESCAPE at M={m}, overshoot={overshoot:+.4f}, gap_lo={gap_lo:.4f}, "
              f"gap_hi={gap_hi:.4f}, eta={eta:.4e} ({int(in_missed.sum())} levels)  [{wall:.1f}s]")
        found.append(dict(name=c["name"], atom=c["atom"], d=d, m_floor=m))
    return found


def weighted_gap(res):
    """WG = min(gap_lo, gap_hi) / eta: the tighter band-edge gap, scaled by the HF weight it has
    to exclude. Larger WG -> predicted easier (more room per unit weight); smaller -> predicted
    harder. One cheap, pre-registered candidate -- not a derivation (see spec SS2, SS7)."""
    return min(res["gap_lo"], res["gap_hi"]) / res["eta"]


def rank_correlation_report(results):
    """Spearman rank correlation of WG and of the overshoot against b_min, across all witnesses
    with a valid (non-None) run_witness result."""
    wg = np.array([weighted_gap(r) for r in results])
    overshoot = np.array([r["overshoot"] for r in results])
    b_min = np.array([r["b_min"] for r in results])
    rho_wg, p_wg = spearmanr(wg, b_min)
    rho_os, p_os = spearmanr(overshoot, b_min)
    return dict(n=len(results), rho_wg=float(rho_wg), p_wg=float(p_wg),
                rho_overshoot=float(rho_os), p_overshoot=float(p_os),
                wg=wg.tolist(), overshoot=overshoot.tolist(), b_min=b_min.tolist())


@dataclass
class NoiseTrial:
    shots: float
    sigma: float
    bounds: np.ndarray
    unsound_fraction: float


def noise_sensitivity_sweep(res, moments, dt, e_min, e_max, M=24, cap=1000.0,
                             shots_list=(1e2, 1e3, 1e4, 1e5, 1e6, 1e8), n_trials=40, seed=0):
    """Perturb moments[1:M] (k=0 held exact -- the known self-overlap, mirroring
    QuantumKrylovSolver._hermitian_shot_noise's unit_diagonal=True) with iid complex Gaussian noise
    of std noise.shot_noise_sigma(shots) per real/imaginary part, n_trials times per shots value.
    Recomputes the bound each time at the fixed (M, cap); reports the empirical distribution and the
    fraction of trials where the (no-longer-exact-moment) bound is unsound (< the true eta)."""
    rng = np.random.default_rng(seed)
    band_lo, band_hi = res["e_d"], res["beta"]
    eta_exact = res["eta"]
    trials = []
    for shots in shots_list:
        sigma = shot_noise_sigma(shots)
        bounds = np.empty(n_trials)
        for t in range(n_trials):
            noisy = moments[:M].copy()
            noise = sigma * (rng.normal(size=M - 1) + 1j * rng.normal(size=M - 1))
            noisy[1:] = noisy[1:] + noise
            ub, _ = _build_moment_bound_lp(noisy, e_min, e_max, dt, band_lo, band_hi, coeff_cap=cap)
            bounds[t] = np.inf if ub is None else ub
        trials.append(NoiseTrial(shots=shots, sigma=sigma, bounds=bounds,
                                  unsound_fraction=float(np.mean(bounds < eta_exact - 1e-9))))
    return trials


def _typical_case_estimate(moments, e_min, e_max, dt, band_lo, band_hi, cap, sigma):
    """||x*||_2 * sigma * sqrt(2(M-1)): a CLT-style 'typical' sensitivity estimate at the LP-optimal
    coefficient vector x*, contrasted with chem-pm3's worst-case (2M-1)*cap*sigma slope."""
    from scripts.spec_pm3_subspace_eta_bound import _moment_objective
    from scipy.optimize import linprog

    M = len(moments)
    obj = _moment_objective(moments)
    k = np.arange(1, M)
    n_grid = 2000
    E = np.linspace(e_min, e_max, n_grid)
    kt = np.outer(dt * E, k)
    A = np.empty((n_grid, 2 * M - 1))
    A[:, 0], A[:, 1::2], A[:, 2::2] = 1.0, 2.0 * np.cos(kt), -2.0 * np.sin(kt)
    ind = ((E >= band_lo) & (E <= band_hi)).astype(float)
    res = linprog(obj, A_ub=-A, b_ub=-ind, bounds=(-cap, cap), method="highs")
    x_norm = float(np.linalg.norm(res.x)) if res.success else float("nan")
    return x_norm * sigma * np.sqrt(2.0 * (M - 1)), x_norm


def main():
    print("=" * 100)
    print("PART A: wider witness set + weighted-gap vs overshoot correlation")
    new_defs = find_escape_witnesses()
    all_defs = PM3_WITNESSES + [dict(name=w["name"], atom=w["atom"], d=w["d"], m_floor=w["m_floor"])
                                for w in new_defs]

    results = []
    for w in all_defs:
        r = run_witness(w, m_extended=(), caps=())
        if r is None:
            continue
        results.append(r)

    report = rank_correlation_report(results)
    print("-" * 100)
    print(f"n={report['n']} witnesses")
    for r in results:
        print(f"  {r['name']}: WG={weighted_gap(r):.4f} overshoot={r['overshoot']:+.4f} "
              f"b_min={r['b_min']:.4f} survives={r['survives']}")
    print(f"Spearman rho(WG, b_min)        = {report['rho_wg']:+.3f} (p={report['p_wg']:.3f})")
    print(f"Spearman rho(overshoot, b_min) = {report['rho_overshoot']:+.3f} "
          f"(p={report['p_overshoot']:.3f})")
    supports = abs(report["rho_wg"]) >= 0.7 and abs(report["rho_wg"]) > abs(report["rho_overshoot"])
    print(f"G2 verdict: {'SUPPORTS' if supports else 'DOES NOT SUPPORT'} the weighted-gap hypothesis "
          "under the pre-registered rule (|rho_WG| >= 0.7 and > |rho_overshoot|)")

    print("=" * 100)
    print("PART B: measured noise sensitivity (witness: linear H6 R=1.1, M=24, cap=1e3)")
    target = next(r for r in results if r["name"] == "linear H6 R=1.1, d=3, M=20")
    mh = build_molecular_hamiltonian(atom=_chain(6, 1.1))
    solver = QuantumKrylovSolver(mh)
    solver._ensure_basis(24)
    moments = np.array([np.vdot(solver._basis[0], v) for v in solver._basis[:24]])
    wv, _, pops = _dense_hf_projection(mh)
    e_min, e_max = float(wv[0]), float(wv[-1])
    reach = np.where(pops > _REACHABLE_TOL)[0]
    reach_e = wv[reach]
    e_d = float(reach_e[3])
    beta = e_d + target["overshoot"]
    target = dict(target, e_d=e_d, beta=beta)

    trials = noise_sensitivity_sweep(target, moments, solver.dt, e_min, e_max)
    worst_case_slope = (2 * 24 - 1) * 1000.0
    print(f"{'shots':>10} {'sigma':>10} {'mean':>10} {'std':>10} {'min':>10} {'max':>10} "
          f"{'unsound%':>9} {'worst-case(slope*sigma)':>24} {'typical-case':>14}")
    for tr in trials:
        typical, x_norm = _typical_case_estimate(moments, e_min, e_max, solver.dt, e_d, beta,
                                                  1000.0, tr.sigma)
        print(f"{tr.shots:>10.0e} {tr.sigma:>10.4f} {tr.bounds.mean():>10.4f} "
              f"{tr.bounds.std():>10.4f} {tr.bounds.min():>10.4f} {tr.bounds.max():>10.4f} "
              f"{100 * tr.unsound_fraction:>8.1f}% {worst_case_slope * tr.sigma:>24.2f} "
              f"{typical:>14.4f}")


if __name__ == "__main__":
    main()
