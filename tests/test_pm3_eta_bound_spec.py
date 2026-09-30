"""
Regression gates for chem-pm3 (specs/BACKLOG.md ~line 237): does a Krylov-moment bound on the HF
spectral weight eta let gamma_corr = sqrt(max(0, gamma_min(beta)^2 - eta)) certify ||P_S u|| for
an overshooting self-mode floor beta? Full finding: scripts/spec_pm3_subspace_eta_bound.py and
sandbox-handoffs/chem-pm3.md (3-witness H6 sweep, dense eigh, ~4 min total -- too slow for CI).

These gates isolate the two claims algebraically/synthetically so they run in under a second and
fail if either piece of logic regresses, without re-running the expensive chemistry:

  G-a (algebraic identity): for ANY d' >= d with S = levels[:d] subset S' = levels[:d'], and ANY
      gamma_min valid for S' (i.e. gamma_min^2 <= ||P_S'||^2), gamma_corr never exceeds the exact
      ||P_S||. This is pure algebra (Cauchy-Schwarz-free: it only uses eta = ||P_S'||^2 -
      ||P_S||^2 >= 0), independent of chemistry, block2, or any physical Hamiltonian -- so it is
      checked with random synthetic probability vectors.

  G-b (LP soundness): the trig-polynomial-majorant LP bound in _build_moment_bound_lp never
      UNDER-estimates the true band weight eta, for a synthetic discrete spectral measure with
      exactly known moments (no dense eigh needed -- moments are the closed-form finite sum
      S_0k = sum_i w_i exp(-i k dt E_i)). Also regression-guards the coefficient-cap fix (see that
      function's docstring): enlarging the box can only relax the LP, so the bound must be
      monotonically non-increasing in the cap.
"""
import numpy as np
import pytest

from scripts.spec_pm3_subspace_eta_bound import _build_moment_bound_lp


def _random_probability_vector(rng, n):
    w = rng.uniform(0.01, 1.0, size=n)
    return w / w.sum()


@pytest.mark.parametrize("seed", range(20))
def test_Ga_gamma_corr_never_exceeds_exact_overlap(seed):
    """G-a: gamma_corr <= exact ||P_S u|| for ANY valid gamma_min on S' and ANY d' >= d.

    eta = ||P_S'||^2 - ||P_S||^2 by construction (S subset S'), so gamma_min^2 - eta <=
    ||P_S'||^2 - eta = ||P_S||^2 always -- this is the entire chem-pm3 G-a claim, checked here
    without any chemistry: only that eta is genuinely the exact excess weight between S and S'.
    """
    rng = np.random.default_rng(seed)
    n = rng.integers(5, 30)
    w = _random_probability_vector(rng, n)
    d = int(rng.integers(1, n - 2))
    d_prime = int(rng.integers(d + 1, n))

    exact_P_S = np.sqrt(w[:d].sum())
    exact_P_S_prime = np.sqrt(w[:d_prime].sum())
    eta = w[:d_prime].sum() - w[:d].sum()
    assert eta >= -1e-12

    # any gamma_min "valid for S'" satisfies gamma_min <= exact_P_S_prime (block Davis-Kahan's
    # actual guarantee); sweep gamma_min from 0 up to that bound, including the tight case.
    for frac in (0.0, 0.3, 0.7, 1.0):
        gamma_min = frac * exact_P_S_prime
        gamma_corr = np.sqrt(max(0.0, gamma_min**2 - eta))
        assert gamma_corr <= exact_P_S + 1e-9, (
            f"G-a violated: seed={seed} d={d} d'={d_prime} frac={frac} "
            f"gamma_corr={gamma_corr} > exact_P_S={exact_P_S}"
        )


def _synthetic_moments(energies, weights, dt, m):
    theta = dt * np.asarray(energies)
    k = np.arange(m)[:, None]
    return (weights[None, :] * np.exp(-1j * k * theta[None, :])).sum(axis=1)


def test_Gb_lp_bound_never_underestimates_synthetic_eta():
    """G-b soundness: the LP's reported eta_bound must be >= the true band weight, for an exactly
    known synthetic discrete spectral measure (no dense eigh -- moments are closed-form)."""
    energies = np.array([0.5, 1.5, 3.2, 4.1, 6.0, 8.5, 9.2])
    weights = np.array([0.10, 0.20, 0.07, 0.08, 0.30, 0.15, 0.10])
    assert np.isclose(weights.sum(), 1.0)
    e_min, e_max = 0.0, 10.0
    dt = np.pi / (e_max - e_min)
    band_lo, band_hi = 3.0, 5.0
    eta_exact = float(weights[(energies >= band_lo) & (energies < band_hi)].sum())  # 0.07 + 0.08

    for m in (4, 8, 12, 20):
        moments = _synthetic_moments(energies, weights, dt, m)
        bound = _build_moment_bound_lp(moments, e_min, e_max, dt, band_lo, band_hi, n_grid=4000)
        assert bound is not None, f"LP failed to solve at M={m}"
        assert bound >= eta_exact - 1e-6, (
            f"G-b unsound at M={m}: bound={bound} < eta_exact={eta_exact}"
        )
        assert bound <= 1.0 + 1e-9


def test_Gb_bound_tightens_with_more_moments():
    """More moments -> a strictly larger function class -> the LP bound must not get worse."""
    energies = np.array([0.5, 1.5, 3.2, 4.1, 6.0, 8.5, 9.2])
    weights = np.array([0.10, 0.20, 0.07, 0.08, 0.30, 0.15, 0.10])
    e_min, e_max = 0.0, 10.0
    dt = np.pi / (e_max - e_min)
    band_lo, band_hi = 3.0, 5.0

    prev = np.inf
    for m in (4, 8, 12, 20):
        moments = _synthetic_moments(energies, weights, dt, m)
        bound = _build_moment_bound_lp(moments, e_min, e_max, dt, band_lo, band_hi, n_grid=4000)
        assert bound <= prev + 1e-6, f"bound got worse from more moments at M={m}"
        prev = bound


def test_Gb_coefficient_cap_is_monotonic_and_safe():
    """Regression for the coefficient-cap fix: enlarging the box can only relax the LP (bound
    non-increasing in cap), and even the tightest cap (1.0, forcing f close to constant) must
    never report a bound below the true eta -- capping can only loosen, never invalidate."""
    energies = np.array([0.5, 1.5, 3.2, 4.1, 6.0, 8.5, 9.2])
    weights = np.array([0.10, 0.20, 0.07, 0.08, 0.30, 0.15, 0.10])
    e_min, e_max = 0.0, 10.0
    dt = np.pi / (e_max - e_min)
    band_lo, band_hi = 3.0, 5.0
    eta_exact = float(weights[(energies >= band_lo) & (energies < band_hi)].sum())
    moments = _synthetic_moments(energies, weights, dt, 12)

    caps = [1.0, 10.0, 100.0, 1000.0, 10000.0]
    bounds = []
    for cap in caps:
        b = _build_moment_bound_lp(moments, e_min, e_max, dt, band_lo, band_hi, n_grid=4000,
                                    coeff_cap=cap)
        assert b is not None
        assert b >= eta_exact - 1e-6, f"cap={cap} produced an unsound (too-tight) bound {b}"
        bounds.append(b)

    for smaller, larger in zip(bounds, bounds[1:]):
        assert larger <= smaller + 1e-6, f"bound not monotonic non-increasing in cap: {bounds}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
