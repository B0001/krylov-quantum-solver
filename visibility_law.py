#!/usr/bin/env python3
"""
The visibility law, made predictive -- a calibrated shot-cost law for spectral lines.

Three specs recorded the same qualitative rule (SPEC_odmd_excited: a mode needs p*sqrt(dm) above
the noise edge; SPEC_device_odmd: immunity ends at the shot floor; SPEC_odmd_spectral: weak
satellites need deeper K). This module turns it into an experiment-planning tool. For the
physical (unnormalized) correlation signal C_k = <psi0|O e^{-i k tau H} O|psi0>, a line of
weight w is detectable iff its Hankel singular value ~ w*sqrt(dm) clears the noise edge
c*sigma*(sqrt d + sqrt m):

    sigma* = w sqrt(dm) / (c (sqrt d + sqrt m))    =>    shots* = 2(2 - 1/dim)/sigma*^2
                                                          ~  1 / (w^2 K).

GATED (tests/test_visibility_law_spec.py): the 50%-detection crossover follows the -2 log-log
slope over FOUR orders of magnitude in w (eight orders in shots: the Nb3I8 optical line costs
~2 shots/element, the near-dark Nb3F8 line ~2e8); a single calibration of the prefactor on one
line predicts every other -- including all three components of a multi-line signal -- to ~10%;
depth buys shots linearly (shots* ~ 1/K). BOUNDARY: line ATTRIBUTION needs a tolerance below
the line spacing, or the strong neighbor masquerades as the weak line (found in probing: 45x
too-early apparent detection at a sloppy tolerance).

HONEST SCOPE (specs/SPEC_visibility_law.md): Gaussian per-element noise with known sigma (the
repo's conventions; fold measured device damping into w before budgeting hardware); the law
prices DETECTION, not precision; c = 1.2 is calibrated, not derived (Hankel noise is
correlated); uniform shot allocation (adaptive schemes could beat it -- a hypothesis, not a bug).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eig, hankel, svd


def detect_line(C, tau: float, sigma: float, e_line: float, seed: int, c: float = 1.2,
                tol_frac: float = 0.03) -> bool:
    """One detection trial on the noisy correlator: the top Hankel singular value must clear the
    noise edge AND a retained DMD pole must sit within ``tol_frac * pi/tau`` of the line
    (centered frame). ``tol_frac`` must keep the tolerance below the line spacing -- see the
    gated attribution boundary."""
    C = np.asarray(C, dtype=complex)
    K = len(C)
    d = K // 2
    rng = np.random.default_rng(seed)
    g = rng.normal(0, sigma / np.sqrt(2), K) + 1j * rng.normal(0, sigma / np.sqrt(2), K)
    s = C + g
    X = hankel(s[:d], s[d - 1:K - 1])
    Xp = hankel(s[1:d + 1], s[d:K])
    U, sig, Vh = svd(X, full_matrices=False)
    edge = c * sigma * (np.sqrt(d) + np.sqrt(K - d))
    if sig[0] <= edge:
        return False
    rank = max(int(np.sum(sig > edge)), 1)
    A = U[:, :rank].conj().T @ Xp @ Vh[:rank, :].conj().T / sig[:rank]
    lam = eig(A, right=False)
    E = -np.angle(lam) / tau
    return bool(np.any(np.abs(E - e_line) < tol_frac * np.pi / tau))


def detection_rate(C, tau: float, shots: float, e_line: float, dim: int, seeds: int = 60,
                   **kw) -> float:
    """Fraction of seeded trials that detect the line at this shot budget."""
    sigma = np.sqrt(2.0 * (2.0 - 1.0 / dim) / shots)
    return float(np.mean([detect_line(C, tau, sigma, e_line, sd, **kw) for sd in range(seeds)]))


def crossover_shots(C, tau: float, e_line: float, dim: int, seeds: int = 60, lo: float = 1e0,
                    hi: float = 1e14, **kw) -> float:
    """The 50%-detection shot budget, by bisection in log-shots (deterministic given seeds)."""
    llo, lhi = np.log10(lo), np.log10(hi)
    if detection_rate(C, tau, 10 ** llo, e_line, dim, seeds, **kw) >= 0.5:
        return float(10 ** llo)
    for _ in range(30):
        mid = 0.5 * (llo + lhi)
        if detection_rate(C, tau, 10 ** mid, e_line, dim, seeds, **kw) >= 0.5:
            lhi = mid
        else:
            llo = mid
    return float(10 ** lhi)


def predicted_shots(w: float, K: int, dim: int, c: float = 1.2,
                    calibration: float = 1.0) -> float:
    """The law: shot budget to detect a line of weight ``w`` in a depth-``K`` correlator.

    ``calibration`` is the measured-over-law ratio from ONE reference line (the transfer gate
    shows a single calibration predicts every other line to ~10%)."""
    d = K // 2
    sigma_star = w * np.sqrt(d * (K - d)) / (c * (np.sqrt(d) + np.sqrt(K - d)))
    return calibration * 2.0 * (2.0 - 1.0 / dim) / sigma_star ** 2


if __name__ == "__main__":
    from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
    from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals
    from odmd_optical import dimer_polarization
    from odmd_spectral import reference_signal

    print("Shot budget to detect each Nb3X8 optical line (K=16, per element, nominal c=1.2):")
    rows = []
    for name, p in NB3X8_LT_BULK.items():
        base = dimer_cluster_integrals(**p)
        mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
        weig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
        pops = np.abs(V.conj().T @ psi_hf) ** 2
        psi0 = V[:, int(np.flatnonzero(pops > 1e-8)[0])]
        s, tau, mu, nrm2 = reference_signal(mh, dimer_polarization() @ psi0, 16)
        star = crossover_shots(nrm2 * s, tau, 0.0, 16)
        rows.append((name, nrm2, star, predicted_shots(nrm2, 16, 16)))
        print(f"  {name}: w={nrm2:9.3e}  measured shots*={star:9.3e}  law={rows[-1][3]:9.3e}"
              f"  ratio={star / rows[-1][3]:.2f}")
    slope = np.polyfit(np.log10([r[1] for r in rows]), np.log10([r[2] for r in rows]), 1)[0]
    print(f"  log-log slope = {slope:.3f}  (law: -2)")
