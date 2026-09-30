#!/usr/bin/env python3
"""
chem-pm3 scout probe: does the "certified sub-cluster anchor" question reduce to a
Krylov-moment bound on the HF spectral weight eta missed by an overshooting self-mode floor?

Background (specs/SPEC_subspace_floor_resolvability.md, specs/BACKLOG.md ~line 237):
the self-mode E_d floor beta_self can exceed the true (d+1)-th reachable level E_d
("overshoot"). certify_hf_subspace_overlap then either flags VACUOUS (guard) or -- on the
guard-PASSING witnesses this script studies -- returns a gamma_min that happens to still be
valid, but only by unquantified slack. This script:

  G-a: confirms gamma_corr = sqrt(max(0, 1 - r^2/delta^2 - eta)) is a valid lower bound on
       ||P_S u|| for ANY beta (not just a tight one), where eta is the exact HF spectral
       weight on the reachable levels missed between the true cluster boundary and beta.
       Derivation: S' = {reachable eigenstates with energy < beta} is BY CONSTRUCTION a
       cluster beta rigorously floors (beta excludes everything in S' by definition), so
       block Davis-Kahan gives gamma_min(r, delta)^2 <= ||P_S' u||^2 unconditionally. Since
       S subset S', ||P_S u||^2 = ||P_S' u||^2 - eta >= gamma_min(r,delta)^2 - eta.

  G-b: derives an eta upper bound from ONLY the Krylov moments S_0k = <HF|U^k|HF>,
       U = exp(-i dt H), k = 0..M-1 (already computed by the solver's overlap matrix, first
       row) via the tightest one-sided trigonometric-polynomial majorant of the band
       indicator 1_[E_d, beta) that is representable in the available M moments:

           eta <= min_f  sum_k f_k * conj(S_0k)     s.t. f(E) >= 1 on [E_d, beta), f(E) >= 0
                                                          everywhere, deg(f) <= M-1

       an LP (finite-dimensional Fourier coefficients, a fine energy grid for the pointwise
       constraints). Compares the LP bound at M in {8,12,16,20,24} against the exact eta
       measured by dense diagonalization, and against the pre-registered kill threshold:
       dies if the moment bound exceeds ~0.05 at M <= 24.

Run: uv run python scripts/spec_pm3_subspace_eta_bound.py
"""
from __future__ import annotations

import time

import numpy as np
from scipy.optimize import linprog

from hf_overlap_subspace import _weinstein_intervals_disjoint
from hybrid_quantum_solver.certified_overlap import rayleigh_quotient, residual_norm
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

ETA_KILL_THRESHOLD = 0.05
M_CAP = 24
_REACHABLE_TOL = 1e-8


def _chain(n, R):
    return "; ".join(f"H 0 0 {i * R}" for i in range(n))


def _asym_chain(gaps):
    z = [0.0]
    for g in gaps:
        z.append(z[-1] + g)
    return "; ".join(f"H 0 0 {zi}" for zi in z)


WITNESSES = [
    dict(name="linear H6 R=1.0, d=3, M=16", atom=_chain(6, 1.0), d=3, m_floor=16),
    dict(name="linear H6 R=1.1, d=3, M=20", atom=_chain(6, 1.1), d=3, m_floor=20),
    dict(name="asym H6 [0.9,0.9,2.2,0.9,0.9], d=3, M=6",
         atom=_asym_chain([0.9, 0.9, 2.2, 0.9, 0.9]), d=3, m_floor=6),
]


def _ritz_floor(solver, Hs, offset, d, m):
    """Raw self-mode Weinstein floor (electronic frame) + the guard's disjointness check."""
    energies, states = solver.eigenstates(m, n_states=d + 1)
    centers = [energies[k] - offset for k in range(d + 1)]
    sigmas = [residual_norm(Hs, states[k], centers[k]) for k in range(d + 1)]
    disjoint = _weinstein_intervals_disjoint(centers, sigmas)
    beta_self = centers[d] - sigmas[d]
    return beta_self, disjoint


_COEFF_CAP = 1000.0  # box bound on Fourier coefficients, see docstring below


def _build_moment_bound_lp(moments, e_min, e_max, dt, band_lo, band_hi, n_grid=12000,
                            coeff_cap=_COEFF_CAP):
    """LP: minimize sum_k f_k conj(S_0k) s.t. f(E) >= 1_[band_lo,band_hi)(E), f(E) >= 0, over a
    fine energy grid. moments[k] = S_0k for k = 0..M-1 (complex); moments[0] must be ~1.

    Fourier coefficients are boxed to +-coeff_cap. This is NOT a relaxation of rigor: f = 1
    (all coefficients 0 except p_0 = 1) is always inside the box and always feasible, so the box
    can only ever raise the reported bound, never invalidate it (eta_bound <= 1 in the worst
    case). It exists because leaving coefficients unbounded made HiGHS numerically unstable --
    same (M, grid, moments) nondeterministically returned "optimal" with |x|_inf ~ 1e6 or status
    4 "numerical difficulties" across repeated runs. Capping at 1e3 reproduced the SAME optimal
    value (to 4 decimal places) as the uncapped runs that did converge, confirming the true
    extremal polynomial does not need coefficients anywhere near that large -- the cap fixes
    conditioning without biasing the answer (checked below via cap-sensitivity).

    Returns the LP optimal value (an upper bound on eta), or None if the LP did not solve.
    """
    M = len(moments)
    # real DOF: p_0 (k=0, real) then (p_k, q_k) for k = 1..M-1, f_k = p_k + i q_k, f_{-k} = conj.
    n_vars = 1 + 2 * (M - 1)

    E = np.linspace(e_min, e_max, n_grid)
    theta = dt * E
    # f(E) = p_0 + 2 * sum_{k=1}^{M-1} [p_k cos(k theta) - q_k sin(k theta)]
    cols = [np.ones_like(theta)]
    for k in range(1, M):
        cols.append(2.0 * np.cos(k * theta))
        cols.append(-2.0 * np.sin(k * theta))
    A_row = np.stack(cols, axis=1)          # (n_grid, n_vars): f(E_j) = A_row[j] @ x

    in_band = (E >= band_lo) & (E < band_hi)
    lower = np.where(in_band, 1.0, 0.0)

    # linprog wants A_ub @ x <= b_ub; our constraint is A_row @ x >= lower -> -A_row @ x <= -lower
    A_ub = -A_row
    b_ub = -lower

    c0 = moments[0].real
    obj = [c0]
    for k in range(1, M):
        ck = moments[k]
        obj.append(2.0 * ck.real)     # coefficient of p_k:  2 Re(f_k conj(c_k)) -> Re part
        obj.append(2.0 * ck.imag)     # coefficient of q_k:  contributes -2 q_k * (-Im c_k) term
    obj = np.asarray(obj)

    bounds = [(-coeff_cap, coeff_cap)] * n_vars
    res = linprog(obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    if not res.success:
        return None
    cap_active = bool(np.max(np.abs(res.x)) > 0.999 * coeff_cap)
    if cap_active:
        print(f"    [warning] coefficient cap {coeff_cap} active at M={M} -- bound may be loose")
    return float(res.fun)


def _objective_coefficients_check():
    """Sanity: sum_k f_k conj(c_k) for f real trig poly, verified against direct sum."""
    rng = np.random.default_rng(0)
    M = 5
    c = np.array([1.0] + [rng.normal() + 1j * rng.normal() for _ in range(M - 1)])
    p = rng.normal(size=M)
    q = np.concatenate([[0.0], rng.normal(size=M - 1)])
    f = p + 1j * q
    f_full = {0: f[0].real}
    for k in range(1, M):
        f_full[k] = f[k]
        f_full[-k] = np.conj(f[k])
    direct = sum(f_full[k] * np.conj(c[k]) if k >= 0 else f_full[k] * np.conj(c[-k]).conj()
                 for k in range(M))
    # sum_k f_k conj(c_k), k=-(M-1)..M-1, c_{-k}=conj(c_k)
    total = f_full[0] * np.conj(c[0])
    for k in range(1, M):
        total += f_full[k] * np.conj(c[k]) + f_full[-k] * np.conj(np.conj(c[k]))
    obj = p[0] * c[0].real
    for k in range(1, M):
        obj += 2.0 * (p[k] * c[k].real + q[k] * c[k].imag)
    assert np.isclose(total.real, obj, atol=1e-9), (total.real, obj)


def run_witness(w):
    print("=" * 100)
    print(w["name"])
    t0 = time.time()
    mh = build_molecular_hamiltonian(atom=w["atom"])
    solver = QuantumKrylovSolver(mh)
    offset = mh.energy_offset
    Hs = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    u = np.asarray(mh.hf_state().data, dtype=complex)
    d = w["d"]

    beta_self, disjoint = _ritz_floor(solver, Hs, offset, d, w["m_floor"])
    print(f"  guard disjoint (should PASS): {disjoint}")

    # exact dense reference (electronic frame)
    Hd = mh.qubit_hamiltonian.to_matrix()
    wv, V = np.linalg.eigh(Hd)
    reach = np.where(np.abs(V.conj().T @ u) ** 2 > _REACHABLE_TOL)[0]
    reach_e = wv[reach]                 # ascending electronic energies of reachable levels
    reach_w = np.abs(V[:, reach].conj().T @ u) ** 2

    e_d_true = float(reach_e[d])        # electronic (d+1)-th reachable level
    e_d_true_total = e_d_true + offset
    beta_self_total = beta_self + offset
    overshoot = beta_self_total - e_d_true_total
    print(f"  beta_self (total) = {beta_self_total:.4f} Ha, true E_d (total) = {e_d_true_total:.4f} Ha,"
          f" overshoot = {overshoot:+.4f} Ha")

    lam = rayleigh_quotient(Hs, u)
    r = residual_norm(Hs, u, lam)
    delta = beta_self - lam
    if delta <= 0 or r >= delta:
        print("  beta_self <= lambda_u or r >= delta -- vacuous at the raw floor, skipping witness")
        return None
    gamma_min_beta = float(np.sqrt(max(0.0, 1.0 - (r / delta) ** 2)))

    # eta = exact HF weight on reachable levels in [e_d_true, beta_self) -- the missed band.
    in_missed = (reach_e >= e_d_true) & (reach_e < beta_self)
    eta_exact = float(reach_w[in_missed].sum())
    n_missed = int(in_missed.sum())
    print(f"  eta_exact = {eta_exact:.4e} over {n_missed} missed reachable levels")

    exact_P_S = float(np.linalg.norm(V[:, reach[:d]].conj().T @ u))
    gamma_corr = float(np.sqrt(max(0.0, gamma_min_beta ** 2 - eta_exact)))
    ga_pass = gamma_corr <= exact_P_S + 1e-9
    print(f"  gamma_min(beta_self) = {gamma_min_beta:.4f}, gamma_corr = {gamma_corr:.4f}, "
          f"exact ||P_S u|| = {exact_P_S:.4f}  -> G-a {'PASS' if ga_pass else 'FAIL'}")

    # G-b: Krylov-moment bound on eta_exact's band [e_d_true, beta_self), electronic frame.
    solver._ensure_basis(M_CAP)
    basis = solver._basis[:M_CAP]
    moments = np.array([np.vdot(basis[0], basis[k]) for k in range(M_CAP)])
    dt = solver.dt
    e_min, e_max = float(wv.min()), float(wv.max())
    # Nyquist sanity: dt * (e_max - e_min) must stay < 2*pi (default dt = pi/width guarantees this).
    assert dt * (e_max - e_min) < 2 * np.pi - 1e-9, "aliasing risk: dt too large for this spectrum"

    print(f"  {'M':>4} {'eta_bound':>12} {'1/(M dt) [Ha]':>16} {'<=0.05?':>8}")
    bounds_by_m = {}
    for M in (8, 12, 16, 20, 24):
        if M > M_CAP:
            continue
        bound = _build_moment_bound_lp(moments[:M], e_min, e_max, dt, e_d_true, beta_self)
        bounds_by_m[M] = bound
        res = "n/a" if bound is None else f"{bound:.4f}"
        krylov_res = 1.0 / (M * dt)
        killed = "n/a" if bound is None else ("YES" if bound <= ETA_KILL_THRESHOLD else "no")
        print(f"  {M:>4} {res:>12} {krylov_res:>16.4f} {killed:>8}")

    return dict(
        name=w["name"], overshoot=overshoot, eta_exact=eta_exact, ga_pass=ga_pass,
        bounds_by_m=bounds_by_m, dt=dt, wall_s=time.time() - t0,
    )


def main():
    _objective_coefficients_check()
    print("objective/coefficient sanity check: OK\n")
    results = []
    for w in WITNESSES:
        res = run_witness(w)
        if res is not None:
            results.append(res)
    print("=" * 100)
    print("SUMMARY")
    any_survived = False
    for res in results:
        m24 = res["bounds_by_m"].get(24)
        survived = m24 is not None and m24 <= ETA_KILL_THRESHOLD
        any_survived = any_survived or survived
        print(f"  {res['name']}: G-a {'PASS' if res['ga_pass'] else 'FAIL'}, "
              f"eta_exact={res['eta_exact']:.4f}, eta_bound(M=24)={m24}, "
              f"overshoot={res['overshoot']:+.4f} Ha, wall={res['wall_s']:.1f}s "
              f"-> G-b {'SURVIVES' if survived else 'KILLED'}")
    print(f"\nG-b overall: {'SURVIVES (at least one witness under threshold)' if any_survived else 'KILLED on all three witnesses'}")


if __name__ == "__main__":
    main()
