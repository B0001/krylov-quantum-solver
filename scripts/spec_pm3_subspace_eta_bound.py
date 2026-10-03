#!/usr/bin/env python3
"""
chem-pm3 driver (specs/SPEC_pm3_eta_bound.md): can the Krylov moments bound the HF spectral weight
eta that an overshooting self-mode floor beta misses?

Background (specs/SPEC_subspace_floor_resolvability.md): the self-mode E_d floor beta_self can
exceed the true (d+1)-th reachable level E_d ("overshoot") on guard-PASSING witnesses. On that
spec's three escape witnesses this script checks:

  G-a: gamma_corr = sqrt(max(0, 1 - r^2/delta^2 - eta)) is a valid lower bound on ||P_S u|| for
       ANY beta, eta = the HF weight on reachable levels in [E_d, beta). Derivation: S' =
       {eigenstates with energy < beta} is BY CONSTRUCTION a cluster beta floors rigorously, so the
       sin-theta bound gives 1 - r^2/delta^2 <= ||P_S' u||^2 unconditionally; S is a subset of S',
       so ||P_S u||^2 = ||P_S' u||^2 - eta >= 1 - r^2/delta^2 - eta.

  G-b: bounds eta from ONLY the Krylov moments S_0k = <HF|U^k|HF>, U = exp(-i dt H), k < M (the
       first row of the solver's overlap matrix) with a one-sided trigonometric majorant:

           eta <= min_f sum_k f_k conj(S_0k)  s.t.  f(E) >= 1_[E_d, beta](E) on [e_min, e_max],
                                                     deg f <= M - 1, |coefficients| <= cap

       (a semi-infinite LP, solved soundly by cutting planes in _build_moment_bound_lp), then
       applies the kill rule pre-registered in the spec: a witness is KILLED if the bound exceeds
       ETA_KILL_THRESHOLD at every M in M_SWEEP.

ORACLE INPUTS (spec section 8): the band edge E_d and the support [e_min, e_max] come from dense
diagonalization. A kill is robust to that (any valid oracle-free band/support is wider, which can
only raise the bound); a survival only means "not ruled out".

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
from reachability import _dense_hf_projection

ETA_KILL_THRESHOLD = 0.05
M_SWEEP = (8, 12, 16, 20, 24)          # the pre-registered kill rule uses min over these
M_EXTENDED = tuple(range(28, 65, 4))   # descriptive only: where a killed witness would pass
CAPS = (1e4, 1e5)                      # descriptive box sensitivity of the M = 24 bound
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


_COEFF_CAP = 1000.0  # box bound on Fourier coefficients, see _build_moment_bound_lp


def _moment_objective(moments):
    """c with c @ x = int f dmu for x = [p_0, p_1, q_1, ..., p_{M-1}, q_{M-1}] and
    f(E) = p_0 + 2 sum_k [p_k cos(k dt E) - q_k sin(k dt E)], moments[k] = int e^{-i k dt E} dmu."""
    m = np.asarray(moments)
    return np.concatenate([[m[0].real], 2.0 * np.column_stack([m[1:].real, m[1:].imag]).ravel()])


def _trig(x, theta, order=0):
    """order-th theta-derivative of f = Re[p_0 + 2 sum_k c_k e^{i k theta}], c_k = p_k + i q_k."""
    c = x[1::2] + 1j * x[2::2]
    k = np.arange(1, len(c) + 1)
    out = 2.0 * (np.exp(1j * np.outer(theta, k)) @ ((1j * k) ** order * c)).real
    return out + x[0] if order == 0 else out


def _stationary_thetas(x):
    """Every stationary point of f on the circle: the angles of ALL roots of z^n f'(theta),
    z = e^{i theta} (an off-circle root only adds a harmless extra candidate, so no root-modulus
    cutoff can drop a real one), raw and Newton-polished."""
    c = x[1::2] + 1j * x[2::2]
    n = len(c)
    k = np.arange(n, -n - 1, -1)                              # descending powers z^(k + n)
    ck = np.concatenate([c[::-1], [x[0]], np.conj(c)])        # c_n .. c_0 .. c_{-n}
    raw = np.angle(np.roots(1j * k * ck))
    pol = raw.copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        for _ in range(3):
            pol = pol - _trig(x, pol, 1) / _trig(x, pol, 2)
    return np.concatenate([raw, pol])


def _gap_candidates(x, dt, e_min, e_max, band_lo, band_hi):
    """(E, g): g = f - 1_[lo, hi] at every point where its minimum on [e_min, e_max] can sit --
    the endpoints of the three pieces and every stationary point of f inside the domain."""
    th0 = dt * e_min
    E = (np.mod(_stationary_thetas(x) - th0, 2.0 * np.pi) + th0) / dt
    E = np.concatenate([E[np.isfinite(E) & (E <= e_max)], [e_min, band_lo, band_hi, e_max]])
    return E, _trig(x, dt * E) - ((E >= band_lo) & (E <= band_hi))


def _build_moment_bound_lp(moments, e_min, e_max, dt, band_lo, band_hi, n_grid=2000,
                           coeff_cap=_COEFF_CAP, tol=1e-8, max_iter=20, tie_break=1e-5):
    """Sound upper bound on mu([band_lo, band_hi]) from moments[k] = int e^{-i k dt E} dmu,
    k = 0..M-1 (moments[0] = mu's mass, ~1), for any measure mu supported on [e_min, e_max].

    Any real trig polynomial f of degree M-1 with f >= 1_band on the WHOLE support gives
    mu(band) <= int f dmu. Enforcing that only on a grid is unsound: the LP optimum dips below the
    indicator between grid points, right at the spectral atoms (the chem-pm3 bug -- the grid-only
    LP, i.e. the first iteration below, gives 0.14717 < eta = 0.15 on the synthetic gate at
    n_grid=4000, M=12, and still 0.14998 at n_grid=40000). So, cutting planes:
    solve on the grid, find the CONTINUOUS minimum of g = f - 1_band (domain/band endpoints plus
    every stationary point of f), add the violated minimizers as constraints, re-solve until
    min g >= -tol; then return int f dmu + max(0, -min g) * mass. Soundness never rests on the
    solver: whatever x HiGHS returns (tolerances, early stop, a failed re-solve), f_x plus its own
    continuous gap is >= 1_band everywhere, so the returned value is a valid bound. The solver
    tolerances (1e-9, below the 1e-7 default the loop otherwise plateaus at) only buy tightness.
    The closed band [lo, hi] is the same constraint as [lo, hi) for continuous f.

    The box |x| <= coeff_cap stabilises HiGHS (unboxed, the same LP returned |x| ~ 1e6 or status
    4 "numerical difficulties"). It only shrinks the feasible set, so it can raise but never
    invalidate the bound; f = 1 is always feasible, so the bound is <= mass up to the gap.

    Stall fallback: when mu has fewer atoms than f has coefficients the optimum is hugely
    degenerate and HiGHS returns wild vertices (f ~ 1e4 in atom-free gaps) whose between-grid dips
    wander instead of shrinking. If the plain pass does not converge, a second pass minimises
    int f dmu + tie_break * (grid mean of f), which selects a tame optimal f. Only the LP's choice
    of f changes (never the feasible set or the returned int f dmu + gap), so it costs tightness
    only, at most tie_break times the grid mean of an optimal f; the min over both passes is kept.

    Returns (bound, lower). lower is the plain grid LP's value at its last point set: a relaxation
    of the exact LP, so lower <= exact LP optimum <= bound brackets the remaining solver slack.
    bound is None if no LP solved.
    """
    M = len(moments)
    obj = _moment_objective(moments)
    k = np.arange(1, M)
    best = lower = None
    for eps in (0.0, tie_break):
        E = np.concatenate([np.linspace(e_min, e_max, n_grid), [band_lo, band_hi]])
        for _ in range(max_iter):
            kt = np.outer(dt * E, k)
            A = np.empty((len(E), 2 * M - 1))   # A @ x = f(E)
            A[:, 0], A[:, 1::2], A[:, 2::2] = 1.0, 2.0 * np.cos(kt), -2.0 * np.sin(kt)
            ind = ((E >= band_lo) & (E <= band_hi)).astype(float)
            res = linprog(obj + eps * A[:n_grid].mean(axis=0), A_ub=-A, b_ub=-ind,
                          bounds=(-coeff_cap, coeff_cap), method="highs",
                          options=dict(primal_feasibility_tolerance=1e-9,
                                       dual_feasibility_tolerance=1e-9))
            if not res.success:
                break
            if eps == 0.0:
                lower = float(res.fun)
            cand, g = _gap_candidates(res.x, dt, e_min, e_max, band_lo, band_hi)
            bound = float(obj @ res.x + max(0.0, -g.min()) * moments[0].real)
            best = bound if best is None else min(best, bound)
            if g.min() >= -tol:
                return best, lower
            E = np.concatenate([E, cand[g < -tol]])
    return best, lower


def run_witness(w, m_extended=M_EXTENDED, caps=CAPS):
    print("=" * 100)
    print(w["name"])
    t0 = time.time()
    mh = build_molecular_hamiltonian(atom=w["atom"])
    solver = QuantumKrylovSolver(mh)
    offset = mh.energy_offset
    Hs = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    u = np.asarray(mh.hf_state().data, dtype=complex)
    d = w["d"]

    beta, disjoint = _ritz_floor(solver, Hs, offset, d, w["m_floor"])
    print(f"  guard disjoint (should PASS): {disjoint}")

    # exact dense reference (electronic frame); ZHEEVR fallback for macOS Accelerate (chem-a0y)
    wv, _, pops = _dense_hf_projection(mh)
    reach = np.where(pops > _REACHABLE_TOL)[0]
    reach_e, reach_w = wv[reach], pops[reach]     # ascending reachable levels, HF weights
    e_d = float(reach_e[d])
    overshoot = beta - e_d
    print(f"  beta_self = {beta + offset:.4f} Ha, true E_d = {e_d + offset:.4f} Ha (total), "
          f"overshoot = {overshoot:+.4f} Ha")

    lam = rayleigh_quotient(Hs, u)
    r = residual_norm(Hs, u, lam)
    delta = beta - lam
    if delta <= 0 or r >= delta:
        print("  beta_self <= lambda_u or r >= delta -- vacuous at the raw floor, skipping witness")
        return None
    gamma_beta = float(np.sqrt(1.0 - (r / delta) ** 2))

    in_missed = (reach_e >= e_d) & (reach_e < beta)
    eta = float(reach_w[in_missed].sum())
    exact = float(np.sqrt(reach_w[:d].sum()))     # ||P_S u||, S = lowest d reachable levels
    gamma_corr = float(np.sqrt(max(0.0, gamma_beta ** 2 - eta)))
    ga_pass = gamma_corr <= exact + 1e-9
    print(f"  eta = {eta:.4e} over {int(in_missed.sum())} missed reachable levels")
    print(f"  gamma_min(beta) = {gamma_beta:.4f}, gamma_corr = {gamma_corr:.4f}, exact ||P_S u|| = "
          f"{exact:.4f} -> G-a {'PASS' if ga_pass else 'FAIL'}")

    # G-b on the band [E_d, beta], electronic frame
    above = reach_e[reach_e >= beta]
    gap_lo, gap_hi = e_d - float(reach_e[d - 1]), float(above[0]) - beta
    print(f"  band edges: E_d - E_(d-1) = {gap_lo:.4f} Ha (weight below {reach_w[d - 1]:.3f}), "
          f"next level - beta = {gap_hi:.4f} Ha")
    m_all = M_SWEEP + tuple(m_extended)
    solver._ensure_basis(m_all[-1])
    moments = np.array([np.vdot(solver._basis[0], v) for v in solver._basis[:m_all[-1]]])
    dt = solver.dt
    e_min, e_max = float(wv[0]), float(wv[-1])
    # Nyquist sanity: dt * (e_max - e_min) must stay < 2*pi (default dt = pi/width guarantees this).
    assert dt * (e_max - e_min) < 2 * np.pi - 1e-9, "aliasing risk: dt too large for this spectrum"

    def bound(M, cap=_COEFF_CAP):
        """(sound bound, relaxed-LP lower value); a failed LP counts as no bound (inf)."""
        ub, lb = _build_moment_bound_lp(moments[:M], e_min, e_max, dt, e_d, beta, coeff_cap=cap)
        return (np.inf if ub is None else ub), (np.nan if lb is None else lb)

    def show(M, b):
        print(f"  {M:>4} {b[0]:>10.4f} {b[1]:>10.4f} {1.0 / (M * dt):>14.4f}")

    print(f"  {'M':>4} {'eta_bound':>10} {'lp_lower':>10} {'1/(M dt) [Ha]':>14}")
    bounds = {}
    for M in M_SWEEP:
        bounds[M] = bound(M)
        show(M, bounds[M])
    b_min = min(b[0] for b in bounds.values())
    survives = b_min <= ETA_KILL_THRESHOLD
    cap_bounds = {cap: bound(M_SWEEP[-1], cap) for cap in caps}
    print("  box sensitivity at M=24, [lp_lower, eta_bound]: " + ", ".join(
        f"cap {c:.0e}: [{lb:.4f}, {ub:.4f}]" for c, (ub, lb) in cap_bounds.items()))
    box_dependent = not survives and any(ub <= ETA_KILL_THRESHOLD for ub, _ in cap_bounds.values())
    gamma_corr_bound = float(np.sqrt(max(0.0, gamma_beta ** 2 - b_min)))
    print(f"  gamma_corr with eta -> min eta_bound: {gamma_corr_bound:.4f} (exact {exact:.4f})")
    m_pass = None
    if not survives:                       # descriptive: the M a killed witness would need
        for M in m_extended:
            bounds[M] = bound(M)
            show(M, bounds[M])
            if bounds[M][0] <= ETA_KILL_THRESHOLD:
                m_pass = M
                break
    return dict(
        name=w["name"], overshoot=overshoot, eta=eta, ga_pass=ga_pass, gamma_corr=gamma_corr,
        exact=exact, gamma_beta=gamma_beta, bounds=bounds, b_min=b_min, survives=survives,
        cap_bounds=cap_bounds, box_dependent=box_dependent, gamma_corr_bound=gamma_corr_bound,
        m_pass=m_pass, dt=dt, gap_lo=gap_lo, gap_hi=gap_hi, n_missed=int(in_missed.sum()),
        wall_s=time.time() - t0,
    )


def main():
    results = [res for w in WITNESSES if (res := run_witness(w)) is not None]
    print("=" * 100)
    print("SUMMARY (kill rule: min over M in 8..24 of eta_bound > 0.05)")
    for res in results:
        gb = ("SURVIVES" if res["survives"] else
              "KILLED (box-dependent)" if res["box_dependent"] else "KILLED")
        print(f"  {res['name']}: G-a {'PASS' if res['ga_pass'] else 'FAIL'}, eta={res['eta']:.4f}, "
              f"min eta_bound={res['b_min']:.4f}, overshoot={res['overshoot']:+.4f} Ha, "
              f"wall={res['wall_s']:.1f}s -> G-b {gb}"
              + ("" if res["survives"] else f", passes at M={res['m_pass']}"))
    n = sum(res["survives"] for res in results)
    print(f"\nG-b: survives on {n} of {len(results)} witnesses")


if __name__ == "__main__":
    main()
