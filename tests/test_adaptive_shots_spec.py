"""
Acceptance gates G1-G4 for specs/SPEC_adaptive_shots_planning.md (non-uniform ODMD shot budgets).

The spec's headline hypothesis -- a decaying schedule S_k ~ e^{-alpha k} resolves the ground state
at a 5-10x lower budget than uniform -- was FALSIFIED during implementation, and the gates below
encode the falsification rather than the claim (specs/README.md step 5). What survives:

  G1  budget preservation (unchanged).
  G2  THE KILL: on an undamped survival amplitude there is no noise-dominated tail to defund;
      |s_k| ~ 0.9 for every k, total noise power sum V/S_k is convex, so uniform is the
      constrained optimum and NO decaying schedule beats it. Gated as a ceiling, not a floor.
  G2b the re-weighting theorem that makes the surviving case work: whitening a GEOMETRIC schedule
      is a geometric rescale, which moves DMD eigenvalue moduli and leaves phases exact (the
      depolarizing-immunity mechanism of device_odmd.py) -- bit-identical energies across alpha.
      A polynomial schedule's weights are not geometric and shift the energy by tens of mHa.
  G3  THE SURVIVOR: under the validated global-depolarizing damping s_k -> f^k s_k a dead tail
      does exist, and the optimal decay rate alpha* rises with the damping rate -log f.
  G4  the size of that win, honestly: ~1.3x in median error (~1.8x in budget), NOT the spec's 5x
      variance reduction -- the error variance is set by a heavy tail of mode misidentifications
      that shot allocation does not control.

All RNG seeded -> deterministic. PySCF/qiskit, no block2; `make gates` isolates it.
"""
import numpy as np

from adaptive_shots import (
    exponential_schedule,
    median_error,
    optimize_decay_factor,
    polynomial_schedule,
    sample_adaptive_odmd_energy,
    whitening_weights,
)
from device_odmd import device_odmd_energy
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from odmd import build_odmd_problem

SHOTS = 1e4
_CACHE = {}


def _n2():
    """N2 CAS(6,6) at 1.1 A -- the spec's system, built once (n=24 covers every gated depth)."""
    if "n2" not in _CACHE:
        mh = build_molecular_hamiltonian(atom="N 0 0 0; N 0 0 1.1",
                                         active_electrons=6, active_orbitals=6)
        _CACHE["n2"] = build_odmd_problem(mh, n=24)
    return _CACHE["n2"]


def test_G1_strict_budget_preservation():
    """Every schedule spends exactly the target budget, with S_0 = 0 (s_0 = 1 needs no shots).

    Exactness is what makes the uniform-vs-adaptive comparison in G2/G4 a matched-budget one.
    """
    for K in (8, 12, 16, 24):
        for sched in ([exponential_schedule(K, SHOTS, a) for a in (-0.2, 0.0, 0.15, 0.5, 1.0)]
                      + [polynomial_schedule(K, SHOTS, g) for g in (0.25, 0.5, 1.0, 2.0)]):
            assert abs(sched.sum() - SHOTS) <= 1.0, (K, sched.sum())
            assert sched[0] == 0.0
            assert np.all(sched[1:] > 0.0)


def test_G2_uniform_is_optimal_on_an_undamped_signal():
    """FALSIFIES the spec's headline: no decaying schedule beats uniform when nothing decays.

    The spec predicted adaptive < 1.0 mHa against uniform > 8.0 mHa at K=12, S=1e4. Measured:
    the survival amplitude does not decay at all (|s_k| >= 0.85 out to k=23), uniform sits at a
    few mHa, and every exponential/polynomial schedule tried is WORSE -- exactly as the convexity
    of sum_k V/S_k under a fixed budget requires. Gated as a ceiling on the achievable gain.
    """
    prob = _n2()
    assert np.abs(prob.s).min() > 0.85, np.abs(prob.s).min()      # no tail to defund

    uniform = median_error(prob, exponential_schedule(12, SHOTS, 0.0), seeds=400)
    assert 2e-3 < uniform < 8e-3, uniform * 1e3                   # spec said > 8 mHa; it is not

    decaying = ([exponential_schedule(12, SHOTS, a) for a in (0.05, 0.1, 0.2, 0.3, 0.5, 0.8)]
                + [polynomial_schedule(12, SHOTS, g) for g in (0.25, 0.5, 1.0, 2.0)])
    errs = np.array([median_error(prob, s, seeds=400) for s in decaying])
    assert errs.min() > 1e-3, errs.min() * 1e3                    # the 1.0 mHa target is unreached
    assert uniform / errs.min() < 1.05, uniform / errs.min()      # and nothing beats uniform


def test_G2b_geometric_whitening_preserves_eigenphases():
    """The re-weighting theorem: geometric weights move |lambda|, never arg(lambda).

    Whitening an exponential schedule rescales s_k by c*beta^k, a similarity that multiplies every
    DMD eigenvalue by beta -- so the noiseless energy is INDEPENDENT of alpha to machine
    precision (residual < 1e-7 Ha is SVD conditioning of the ill-conditioned noiseless Hankel --
    it does not grow with alpha). Polynomial weights are not geometric and bias the energy by tens
    of mHa, five orders larger, which is why the adaptive path is specced on the exponential
    family.
    """
    prob = _n2()
    for K in (12, 16, 20):
        ref = device_odmd_energy(prob.s[:K], prob.tau, 1e-12, amp_floor=0.0)
        for a in (0.1, 0.15, 0.3, 0.5):
            r = whitening_weights(exponential_schedule(K, SHOTS, a))
            e = device_odmd_energy(prob.s[:K] * r, prob.tau, 1e-12, amp_floor=0.0)
            assert abs(e - ref) < 1e-7, (K, a, e - ref)      # SVD conditioning, not bias
        biases = [abs(device_odmd_energy(prob.s[:K] * whitening_weights(
            polynomial_schedule(K, SHOTS, g)), prob.tau, 1e-12, amp_floor=0.0) - ref)
            for g in (0.5, 1.0)]
        assert max(biases) > 1e-2, (K, biases)                    # > 10 mHa: not phase-preserving


def test_G3_optimal_decay_tracks_the_damping_rate():
    """THE SURVIVING CLAIM: alpha* is ~0 with no damping and rises with the damping rate -log f.

    Damping is the repo's validated global-depolarizing model s_k -> f^k s_k (device_odmd.py),
    which leaves the exact eigenphases invariant, so alpha* responds to the SIGNAL-TO-NOISE
    profile alone and not to a moving target. Gated as monotone non-decreasing plus a strict
    separation between the undamped and damped ends -- the optimum is broad, so the grid
    position, not a fitted vertex, is the robust statistic (stable across 150-600 seeds).
    """
    prob = _n2()
    alphas = [optimize_decay_factor(prob, 24, SHOTS, damping=f, seeds=300)
              for f in (1.0, 0.9, 0.8)]
    assert alphas == sorted(alphas), alphas
    assert alphas[0] <= 0.05, alphas
    assert min(alphas[1], alphas[2]) >= 0.10, alphas


def test_G4_the_win_is_1p3x_in_error_and_zero_in_variance():
    """Sizes the surviving win honestly against the spec's 5x variance-reduction claim.

    1000 independent realizations, K=24, matched budget S=1e4. Under damping the adaptive
    schedule cuts the MEDIAN error by >= 1.25x (~1.6x in budget at the 1/sqrt(S) scaling) --
    real, but far from 5-10x. The error VARIANCE barely moves (< 1.5x, not the specced >= 5x):
    it is dominated by rare mode misidentifications, which no shot allocation fixes. Undamped,
    the same schedule is strictly worse -- the G2 kill, re-measured on the tail statistics.
    """
    prob = _n2()

    def errs(alpha, f):
        e = np.array([sample_adaptive_odmd_energy(prob, exponential_schedule(24, SHOTS, alpha),
                                                  sd, damping=f) for sd in range(1000)])
        return np.abs(np.where(np.isfinite(e), e, prob.ref + 1.0) - prob.ref)

    damped_u, damped_a = errs(0.0, 0.8), errs(0.15, 0.8)
    assert np.median(damped_u) / np.median(damped_a) > 1.25, \
        np.median(damped_u) / np.median(damped_a)
    assert np.var(damped_u) / np.var(damped_a) < 1.5, np.var(damped_u) / np.var(damped_a)

    undamped_u, undamped_a = errs(0.0, 1.0), errs(0.15, 1.0)
    assert np.median(undamped_u) / np.median(undamped_a) < 1.0, \
        np.median(undamped_u) / np.median(undamped_a)
