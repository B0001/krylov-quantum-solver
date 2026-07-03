#!/usr/bin/env python3
"""
Coverage-gated error bars for ODMD -- a union bootstrap from a single signal.

Every ODMD number elsewhere in this repo is validated as a median over noise seeds against a
known FCI answer. A real experiment gets ONE noisy signal and no ground truth. This module puts
an error bar on that single signal, and the error bar itself is validated the only way error
bars can be: empirical COVERAGE over independent realizations (tests/test_odmd_uq_spec.py).

The interval is the UNION of two resampling arms, because each alone is broken in a
complementary regime (the recorded finding, specs/SPEC_odmd_uq.md):

  * parametric bootstrap -- fit modes, rebuild the clean signal, re-noise at the KNOWN
    per-element sigma, re-estimate. Anti-conservative by up to 18x on many-mode signals: the DMD
    fit absorbs the realized noise into its modes, so the resampled ensemble is artificially
    stable and never sees the threshold rank-switching instability (mean-shift bias correction
    does not fix it -- probed).
  * BOP-DMD-style bagging (cf. arXiv:2107.10878) -- random subsets of the Hankel columns (the
    DMD least-squares is over columns, so X'_sub = A X_sub stays exact). Sees the fit
    instability, but under-spreads few-mode signals (H2).

HONEST SCOPE: the union is conservative (~2-3x the true spread) -- an anti-conservative error
bar is worse than none, so that trade is deliberate. sigma must be KNOWN (shot budget). And no
resampling of one signal can see model-misspecification bias: a too-shallow K leaves a
truncation bias the interval will confidently miss (coverage 0 -- the gated boundary). Pair any
interval with a depth-convergence check. Validated at K=24, three systems, Gaussian noise,
alpha = 0.1.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy.linalg import eig, hankel, svd

from odmd import _dmd_modes, odmd_energy


@dataclass
class ODMDInterval:
    """A single-signal confidence interval (centered frame, same units as the signal's H)."""
    estimate: float                     # the pinned odmd_energy point estimate
    lower: float                        # union interval
    upper: float
    parametric: Tuple[float, float]     # component intervals (for diagnostics/gates)
    bagging: Tuple[float, float]
    sigma: float
    alpha: float

    @property
    def half_width(self) -> float:
        return 0.5 * (self.upper - self.lower)


def _bag_energy(s, tau, sigma, cols):
    """Truncated-SVD DMD ground eigenphase using only the Hankel columns in ``cols``."""
    K = len(s)
    d = K // 2
    X = hankel(s[:d], s[d - 1:K - 1])[:, cols]
    Xp = hankel(s[1:d + 1], s[d:K])[:, cols]
    U, sig, Vh = svd(X, full_matrices=False)
    rank = max(int(np.sum(sig > 5.0 * sigma * sig[0])), 1)
    A = U[:, :rank].conj().T @ Xp @ Vh[:rank, :].conj().T / sig[:rank]
    lam = eig(A, right=False)
    keep = np.abs(np.abs(lam) - 1.0) < 0.2
    if not np.any(keep):                                  # degenerate guard (as odmd_energy)
        keep = np.ones_like(lam, dtype=bool)
    return float(np.min(-np.angle(lam[keep]) / tau))


def odmd_confidence_interval(s, tau: float, sigma: float, n_resamples: int = 200,
                             alpha: float = 0.1, seed: int = 0) -> ODMDInterval:
    """Union-bootstrap CI for the ODMD ground energy from ONE noisy signal.

    ``sigma`` is the known per-element noise scale (~sqrt(2(2-1/d)/shots), the repo's
    Hadamard-test convention). The point estimate and both arms use the pinned noise-aware
    threshold 5*sigma. Deterministic given ``seed``.
    """
    s = np.asarray(s, dtype=complex)
    K = len(s)
    m = K - K // 2
    e_hat = odmd_energy(s, tau, svd_threshold=5.0 * sigma)[0]
    rng = np.random.default_rng(seed)

    # parametric arm: refit -> rebuild -> re-noise at the known sigma
    lam, _ = _dmd_modes(s, tau, mod_window=0.2, cutoff_rel=5.0 * sigma)
    V = np.vander(lam, N=K, increasing=True).T
    a, *_ = np.linalg.lstsq(V, s, rcond=None)
    s_hat = V @ a
    s_hat[0] = 1.0
    boots = np.empty(n_resamples)
    for b in range(n_resamples):
        g = rng.normal(0, sigma / np.sqrt(2), K) + 1j * rng.normal(0, sigma / np.sqrt(2), K)
        g[0] = 0.0
        boots[b] = odmd_energy(s_hat + g, tau, svd_threshold=5.0 * sigma)[0]
    p_lo, p_hi = (float(q) for q in np.quantile(boots, [alpha / 2, 1 - alpha / 2]))

    # bagging arm: random 60% column subsets, recentered at the full-data estimate
    bags = np.empty(n_resamples)
    for b in range(n_resamples):
        cols = np.sort(rng.choice(m, size=max(3, int(0.6 * m)), replace=False))
        bags[b] = _bag_energy(s, tau, sigma, cols)
    b_lo, b_hi = np.quantile(bags, [alpha / 2, 1 - alpha / 2])
    med = float(np.median(bags))
    b_lo, b_hi = e_hat + float(b_lo) - med, e_hat + float(b_hi) - med

    return ODMDInterval(estimate=e_hat, lower=min(p_lo, b_lo), upper=max(p_hi, b_hi),
                        parametric=(p_lo, p_hi), bagging=(b_lo, b_hi),
                        sigma=float(sigma), alpha=float(alpha))


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from odmd import build_odmd_problem

    mh = build_molecular_hamiltonian(atom="N 0 0 0; N 0 0 1.1",
                                     active_electrons=6, active_orbitals=6)
    prob = build_odmd_problem(mh, n=24)
    shots = 1e5
    sigma = np.sqrt(2.0 * (2.0 - 1.0 / prob.dim) / shots)
    print("N2 CAS(6,6), K=24, 1e5 shots -- one signal, one interval (nominal 90%):")
    rng = np.random.default_rng(42)
    g = rng.normal(0, sigma / np.sqrt(2), 24) + 1j * rng.normal(0, sigma / np.sqrt(2), 24)
    g[0] = 0.0
    ci = odmd_confidence_interval(prob.s + g, prob.tau, sigma, seed=42)
    print(f"  E = {ci.estimate + prob.offset:.6f}  "
          f"CI = [{ci.lower + prob.offset:.6f}, {ci.upper + prob.offset:.6f}] Ha "
          f"(half-width {ci.half_width * 1e3:.3f} mHa)")
    truth = prob.ref + prob.offset
    print(f"  exact (would be unknown in an experiment): {truth:.6f} Ha -- "
          f"{'inside' if ci.lower + prob.offset <= truth <= ci.upper + prob.offset else 'OUTSIDE'}")
