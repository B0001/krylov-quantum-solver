"""
Acceptance gates G1-G4 for specs/SPEC_odmd_uq.md (coverage-gated single-signal error bars).

Test-first: ``odmd_uq`` does not exist yet, so this file is RED until the spec is implemented.
Claim: from ONE noisy survival signal (no ground truth), the union of a parametric bootstrap
(refit modes -> rebuild -> re-noise at the KNOWN sigma) and a BOP-DMD-style bagging ensemble
(random Hankel-column subsets) gives a 90% CI whose empirical coverage meets nominal. Error bars
are falsifiable in exactly one way -- coverage over independent realizations -- and that is the
gate. Recorded findings: each arm alone is broken in a complementary regime (parametric is
anti-conservative up to 18x -- the fit absorbs realized noise and never sees threshold
rank-switching; bagging under-spreads few-mode signals), and NO resampling of one signal can see
model-misspecification bias (the K=8 boundary).

All RNGs seeded -> deterministic gates. PySCF/qiskit, no block2; `make gates` isolates it.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from odmd import build_odmd_problem, odmd_energy
from odmd_uq import odmd_confidence_interval

SYSTEMS = {
    "h2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "h4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "n2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
CONFIGS = [("h2", 1e4), ("h2", 1e5), ("h4", 1e4), ("h4", 1e5), ("n2", 1e4), ("n2", 1e5)]
_PROBS, _COV = {}, {}


def _prob(key):
    if key not in _PROBS:
        _PROBS[key] = build_odmd_problem(build_molecular_hamiltonian(**SYSTEMS[key]), n=24)
    return _PROBS[key]


def _coverage(key, K, shots, trials=200):
    """(union, parametric, bagging coverage; median half-width; median |err|) -- cached."""
    ck = (key, K, shots)
    if ck not in _COV:
        prob = _prob(key)
        sigma = np.sqrt(2.0 * (2.0 - 1.0 / prob.dim) / shots)
        cu = cp = cb = 0
        widths, errs = [], []
        for sd in range(trials):
            rng = np.random.default_rng(10_000 + sd)
            g = rng.normal(0, sigma / np.sqrt(2), K) + 1j * rng.normal(0, sigma / np.sqrt(2), K)
            g[0] = 0.0
            s = prob.s[:K] + g
            ci = odmd_confidence_interval(s, prob.tau, sigma, seed=sd)
            cu += int(ci.lower <= prob.ref <= ci.upper)
            cp += int(ci.parametric[0] <= prob.ref <= ci.parametric[1])
            cb += int(ci.bagging[0] <= prob.ref <= ci.bagging[1])
            widths.append(ci.half_width)
            errs.append(abs(odmd_energy(s, prob.tau, svd_threshold=5 * sigma)[0] - prob.ref))
        _COV[ck] = (cu / trials, cp / trials, cb / trials,
                    float(np.median(widths)), float(np.median(errs)))
    return _COV[ck]


def test_G1_union_coverage_meets_nominal_everywhere():
    """DEFINITION OF DONE: nominal-90% union CI covers the true energy >= 85% of the time on
    every system x budget (measured 0.895-1.000)."""
    for key, shots in CONFIGS:
        cu, _, _, _, _ = _coverage(key, 24, shots)
        assert cu >= 0.85, (key, shots, cu)


def test_G2_width_is_informative():
    """Median half-width within [1x, 4x] the median true |error| -- conservative, not vacuous."""
    for key, shots in CONFIGS:
        _, _, _, w, e = _coverage(key, 24, shots)
        assert 1.0 <= w / e <= 4.0, (key, shots, w / e)


def test_G3_union_is_load_bearing():
    """Each arm alone fails in a complementary regime: parametric anti-conservative on N2 1e5
    (< 0.5 at nominal 0.9 -- the fit absorbs realized noise); bagging under-spreads the 2-mode
    H2 signal (< 0.85). Neither component suffices."""
    _, cp, _, _, _ = _coverage("n2", 24, 1e5)
    assert cp < 0.5, cp
    _, _, cb, _, _ = _coverage("h2", 24, 1e5)
    assert cb < 0.85, cb


def test_G4_bias_is_invisible_the_boundary():
    """N2 at K=8, 1e6 shots: truncation bias (~8.7 mHa) >> noise -- NO resampling of the single
    signal can see it, and coverage collapses (< 0.1). Pair intervals with a depth-convergence
    check before trusting them."""
    cu, _, _, _, _ = _coverage("n2", 8, 1e6)
    assert cu < 0.1, cu
