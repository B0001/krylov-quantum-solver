#!/usr/bin/env python3
"""
The certified energy bracket under shot noise -- when the guarantee becomes probabilistic.

The certified arc (`temple_bounds`, `certified_gaps`, `certified_dipole`, `certified_thermochem`) is
exact-statevector: it answers "can you certify without FCI?" but not "does the certificate survive
finite sampling?" -- the question that decides whether any of it is usable on hardware. This module
answers it, and the answer is surprising.

Both the variational upper bound rho_0 = <H> and the Temple lower bound tau_0 need <H> and <H^2>.
Under N shots these are estimates with i.i.d. standard errors set by the Hamiltonian 1-norms:

    se(<H>)   ~ lambda_H  / sqrt(N)          se(<H^2>) ~ lambda_{H^2} / sqrt(N)

(lambda = sum |Pauli coefficients|). A Monte-Carlo over noise realizations measures coverage of the
exact reachable E_0.

THE FINDINGS (specs/SPEC_certified_noise.md):
  * SAMPLING BREAKS THE CERTIFICATE. At converged depth the raw bracket covers E_0 only ~0.40 of the
    time, and the *variational upper bound holds ~0.50* -- a coin flip -- because rho_0 -> E_0 makes
    symmetric noise push the estimate below E_0 half the time. The tighter the Ritz state, the more
    fragile: a converged variational "bound" is NOT a bound under sampling.
  * SHOTS DO NOT BUY COVERAGE. Raw coverage is ~N-INDEPENDENT (unchanged from N=1e4 to 1e8): the
    variational knife-edge is structural, not a finite-sample effect you can shoot your way out of.
  * INFLATION BUYS COVERAGE, SHOTS BUY TIGHTNESS. Widening the bracket by z * standard-error restores
    coverage >= 0.9 (conservative by design, cf. `odmd_uq`), and the inflated half-width scales as
    z * lambda_H / sqrt(N) -- so more shots shrink the certified interval as 1/sqrt(N) (the shot-cost
    law, cf. the visibility law) without ever restoring the raw guarantee.

HONEST SCOPE: i.i.d. Gaussian shot noise with lambda-1-norm standard errors (a standard idealization;
real grouped-Pauli measurement with covariances differs by O(1) factors); an ORACLE gap E_1 feeds
Temple (making eps_1 noisy only worsens the break, so the finding is conservative); the model cannot
see systematic (Trotter/basis) bias. <H^2> carries the larger 1-norm (lambda_{H^2} >> lambda_H), so
the Temple lower bound is the noise-expensive side.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import _mean_and_variance


def hamiltonian_one_norms(mh: MolecularHamiltonian) -> Tuple[float, float]:
    """(lambda_H, lambda_{H^2}) -- the Pauli 1-norms of H and H^2 (electronic frame)."""
    Hop = mh.qubit_hamiltonian
    lam_h = float(np.abs(Hop.coeffs).sum())
    lam_h2 = float(np.abs((Hop @ Hop).simplify().coeffs).sum())
    return lam_h, lam_h2


def reachable_E0_E1(mh: MolecularHamiltonian) -> Tuple[float, float]:
    """REFERENCE (dense): the two lowest HF-reachable eigenvalues (electronic frame) -- the coverage
    target E_0 and the oracle gap input E_1."""
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    reach = np.sort(w[np.abs(V.conj().T @ hf) ** 2 > 1e-10])
    return float(reach[0]), float(reach[1])


def certified_half_width(lam_h: float, shots: float, z: float = 2.0) -> float:
    """The inflated upper half-width z * lambda_H / sqrt(N) -- the shot-cost of a probabilistic
    certificate (Ha). Shrinks as 1/sqrt(N)."""
    return z * lam_h / np.sqrt(shots)


def shot_noise_coverage(mh: MolecularHamiltonian, m: int, shots: float, trials: int = 4000,
                        z: float = 2.0, seed: int = 0,
                        solver: Optional[QuantumKrylovSolver] = None) -> dict:
    """Monte-Carlo coverage of the exact E_0 by the noisy certified bracket at N=``shots``.

    Returns raw two-sided coverage, the variational-upper-bound hit rate, the z*se inflated
    coverage, and the inflated half-width (Ha). Uses the exact Ritz state and an oracle gap; only
    the <H>, <H^2> measurements are sampled."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    psi0 = solver.eigenstates(m, n_states=1)[1][0]
    th0, var0 = _mean_and_variance(H, psi0)
    h2_exact = var0 + th0 * th0
    E0, E1 = reachable_E0_E1(mh)
    lam_h, lam_h2 = hamiltonian_one_norms(mh)
    se_h, se_h2 = lam_h / np.sqrt(shots), lam_h2 / np.sqrt(shots)

    rng = np.random.default_rng(seed)
    th = th0 + rng.normal(0.0, se_h, trials)
    h2 = h2_exact + rng.normal(0.0, se_h2, trials)
    var = np.clip(h2 - th * th, 0.0, None)
    tau = np.where(E1 > th, th - var / (E1 - th), -np.inf)          # noisy Temple lower
    cov_raw = float(np.mean((tau <= E0) & (E0 <= th)))
    cov_upper = float(np.mean(E0 <= th))                            # variational bound hit rate
    th_lo, th_hi = th - z * se_h, th + z * se_h                     # inflated
    var_hi = np.clip(h2 + z * se_h2 - th_lo * th_lo, 0.0, None)
    tau_i = np.where(E1 > th_lo, th_lo - var_hi / (E1 - th_lo), -np.inf)
    cov_inflated = float(np.mean((tau_i <= E0) & (E0 <= th_hi)))
    return dict(cov_raw=cov_raw, cov_upper=cov_upper, cov_inflated=cov_inflated,
                half_width=certified_half_width(lam_h, shots, z), lam_h=lam_h, lam_h2=lam_h2)


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H2": ("H 0 0 0; H 0 0 0.74", 6),
        "H4": ("H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7", 12),
    }
    for name, (atom, m) in cases.items():
        mh = build_molecular_hamiltonian(atom=atom)
        solver = QuantumKrylovSolver(mh)
        lam_h, lam_h2 = hamiltonian_one_norms(mh)
        print("=" * 78)
        print(f"{name}: lambda_H={lam_h:.2f}, lambda_H2={lam_h2:.2f} (H^2 noisier)")
        print("   N shots | raw cov | var-upper holds | inflated cov (z=2) | inflated +/- (mHa)")
        for shots in (1e4, 1e6, 1e8):
            r = shot_noise_coverage(mh, m, shots, solver=solver)
            print(f"   {shots:.0e} |  {r['cov_raw']:.3f}  |      {r['cov_upper']:.3f}      "
                  f"|       {r['cov_inflated']:.3f}        | {r['half_width'] * 1e3:8.2f}")
    print("=" * 78)
    print("Raw coverage ~0.4 and N-independent: shots do not buy back the variational guarantee;")
    print("inflation buys coverage (>=0.9), shots buy tightness (half-width ~ lambda/sqrt(N)).")
