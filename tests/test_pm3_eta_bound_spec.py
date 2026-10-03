"""
Gates for chem-pm3 (specs/SPEC_pm3_eta_bound.md): does a Krylov-moment bound on the HF spectral
weight eta let gamma_corr = sqrt(max(0, gamma_min(beta)^2 - eta)) certify ||P_S u|| for an
overshooting self-mode floor beta? Driver and full numbers: scripts/spec_pm3_subspace_eta_bound.py.

  G-a (algebraic identity): for ANY d' >= d with S = levels[:d] subset S' = levels[:d'], and ANY
      gamma_min valid for S' (i.e. gamma_min^2 <= ||P_S'||^2), gamma_corr never exceeds the exact
      ||P_S||. Pure algebra (eta = ||P_S'||^2 - ||P_S||^2 >= 0), checked on random synthetic
      probability vectors.

  G-b soundness: the trig-polynomial-majorant LP bound in _build_moment_bound_lp never
      UNDER-estimates the true band weight, for a synthetic discrete spectral measure with exactly
      known moments S_0k = sum_i w_i exp(-i k dt E_i). Its objective is the spectral integral, its
      continuous-minimum finder misses no dip (brute force), and it is monotone in M and in the
      coefficient box (enlarging either only relaxes the LP).

  Witnesses: G-a, bound >= exact eta, and the pre-registered G-b verdict on the three real escape
      witnesses (dense 12-qubit diagonalization; ~3 min for the three).
"""
import functools

import numpy as np
import pytest

from scripts.spec_pm3_subspace_eta_bound import (
    WITNESSES,
    _build_moment_bound_lp,
    _gap_candidates,
    _moment_objective,
    _trig,
    run_witness,
)

ENERGIES = np.array([0.5, 1.5, 3.2, 4.1, 6.0, 8.5, 9.2])
WEIGHTS = np.array([0.10, 0.20, 0.07, 0.08, 0.30, 0.15, 0.10])
E_MIN, E_MAX = 0.0, 10.0
DT = np.pi / (E_MAX - E_MIN)
BAND_LO, BAND_HI = 3.0, 5.0
ETA_EXACT = float(WEIGHTS[(ENERGIES >= BAND_LO) & (ENERGIES < BAND_HI)].sum())  # 0.07 + 0.08


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


def _synthetic_moments(m):
    k = np.arange(m)[:, None]
    return (WEIGHTS[None, :] * np.exp(-1j * k * DT * ENERGIES[None, :])).sum(axis=1)


@functools.cache
def _bound_and_lower(m, coeff_cap):
    return _build_moment_bound_lp(_synthetic_moments(m), E_MIN, E_MAX, DT, BAND_LO, BAND_HI,
                                  n_grid=4000, coeff_cap=coeff_cap)


def _bound(m, coeff_cap=1000.0):
    return _bound_and_lower(m, coeff_cap)[0]


def test_Gb_objective_is_the_spectral_integral():
    """The identity the bound rests on: obj @ x = int f dmu = sum_i w_i f(E_i)."""
    rng = np.random.default_rng(0)
    for m in (2, 5, 12):
        x = 100.0 * rng.normal(size=2 * m - 1)
        assert np.isclose(_moment_objective(_synthetic_moments(m)) @ x,
                          WEIGHTS @ _trig(x, DT * ENERGIES), rtol=1e-12, atol=1e-9)


@pytest.mark.parametrize("seed", range(6))
def test_Gb_continuous_minimum_misses_no_dip(seed):
    """The certified gap is only as good as the minimum finder: the min of f - 1_band over the
    stationary-point candidates must be <= a 2e5-point brute-force min (no dip missed), at the
    coefficient scales the LP actually uses (|x| up to the 1e3 box)."""
    rng = np.random.default_rng(seed)
    x = (1.0 if seed % 2 else 1000.0) * rng.normal(size=2 * 12 - 1)
    _, g = _gap_candidates(x, DT, E_MIN, E_MAX, BAND_LO, BAND_HI)
    E = np.linspace(E_MIN, E_MAX, 200_001)
    brute = (_trig(x, DT * E) - ((E >= BAND_LO) & (E <= BAND_HI))).min()
    assert g.min() <= brute + 1e-12 * np.abs(x).sum(), (g.min(), brute)


def test_Gb_lp_bound_never_underestimates_synthetic_eta():
    """G-b soundness: the LP's reported eta_bound must be >= the true band weight, for an exactly
    known synthetic discrete spectral measure (no dense eigh -- moments are closed-form). The
    grid-only LP failed this at M=12 (0.14717 < 0.15). Also: the plain grid LP is a relaxation of
    the exact LP, so its value can never exceed a sound bound (an under-certified gap would)."""
    for m in (4, 8, 12, 20):
        bound, lower = _bound_and_lower(m, 1000.0)
        assert bound is not None, f"LP failed to solve at M={m}"
        assert bound >= ETA_EXACT - 1e-6, (
            f"G-b unsound at M={m}: bound={bound} < eta_exact={ETA_EXACT}"
        )
        assert bound <= 1.0 + 1e-9
        assert lower <= bound + 1e-6, f"relaxed LP {lower} above the 'sound' bound {bound}"


def test_Gb_bound_tightens_with_more_moments():
    """More moments -> a strictly larger function class -> the LP bound must not get worse."""
    prev = np.inf
    for m in (4, 8, 12, 20):
        bound = _bound(m)
        assert bound <= prev + 1e-6, f"bound got worse from more moments at M={m}"
        prev = bound


def test_Gb_coefficient_cap_is_monotonic_and_safe():
    """Regression for the coefficient-cap fix: enlarging the box can only relax the LP (bound
    non-increasing in cap), and even the tightest cap (1.0, forcing f close to constant) must
    never report a bound below the true eta -- capping can only loosen, never invalidate."""
    caps = [1.0, 10.0, 100.0, 1000.0, 10000.0]
    bounds = []
    for cap in caps:
        b = _bound(12, coeff_cap=cap)
        assert b is not None
        assert b >= ETA_EXACT - 1e-6, f"cap={cap} produced an unsound (too-tight) bound {b}"
        bounds.append(b)

    for smaller, larger in zip(bounds, bounds[1:]):
        assert larger <= smaller + 1e-6, f"bound not monotonic non-increasing in cap: {bounds}"


# name -> (G-b survives, min eta_bound over M in 8..24): the pre-registered verdicts recorded in
# SPEC_pm3_eta_bound.md section 10 (scripts/spec_pm3_subspace_eta_bound.py main()).
RECORDED = {
    "linear H6 R=1.0, d=3, M=16": (True, 0.0396),
    "linear H6 R=1.1, d=3, M=20": (False, 0.0594),
    "asym H6 [0.9,0.9,2.2,0.9,0.9], d=3, M=6": (True, 0.0359),
}


@pytest.fixture(scope="module")
def witnesses():
    return {w["name"]: run_witness(w, m_extended=(), caps=()) for w in WITNESSES}


@pytest.mark.parametrize("name", list(RECORDED))
def test_witness_Ga_soundness_and_Gb_verdict(witnesses, name):
    """On each real escape witness (dense reference, ~1 min each): G-a holds, every moment bound
    is >= the exact eta and >= the relaxed grid LP, and the verdict and min bound are the ones
    recorded in the spec."""
    res = witnesses[name]
    assert res["ga_pass"], f"G-a: gamma_corr {res['gamma_corr']} > exact {res['exact']}"
    for M, (ub, lb) in res["bounds"].items():
        assert ub >= res["eta"], f"M={M}: bound {ub} below the exact eta {res['eta']}"
        assert lb <= ub + 1e-6, f"M={M}: relaxed LP {lb} above the sound bound {ub}"
    survives, b_min = RECORDED[name]
    assert res["survives"] == survives
    assert res["b_min"] == pytest.approx(b_min, abs=5e-4)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
