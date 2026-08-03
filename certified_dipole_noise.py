#!/usr/bin/env python3
"""
Certified dipole under shot noise -- inflation cannibalizes the gap margin it needs to stay finite.

`certified_dipole`'s property bracket rests on Davis-Kahan: sin(theta) <= sigma_0/Delta_lo, so its
half-width is finite only while Delta_lo (the certified gap lower bound, `certified_gaps`) exceeds
sigma_0. `certified_noise`/`certified_thermochem_noise`/`gap_selfcheck_noise` already showed z*se
inflation is the standard fix for shot-noise coverage collapse, monotonically in both prior
compositions (difference: less z needed; intersection: more z needed). This module asks the same
question of a THIRD composition -- one that spends the SAME resource (Delta_lo's margin) twice: once
to keep the certificate finite, once as the thing inflation must conservatively shrink to stay
honest. The result is not monotonic in z.

One noisy realization per trial of (theta0, var0, theta1, var1) gives raw sigma_0, Delta_lo (the
self-mode `certified_gaps.gap_bracket` formula, unmodified); one noisy realization of (mu, <A^2>)
(the dipole operator's own Pauli 1-norms, `operator_one_norms`, mirroring
`certified_noise.hamiltonian_one_norms`) gives raw mu, sigma_A. Padding is applied post-hoc on the
already-computed raw quantities -- the same safe, unambiguous recipe `gap_selfcheck_noise` settled
on after finding internal directional perturbation genuinely ambiguous for composed quantities (here
even more so: three composed noisy terms, not two). A trial is vacuous (s = sigma_0/Delta_lo >= 1
after padding) is scored as NOT covered and tracked separately via ``finite_frac`` -- the same
honest convention `gap_selfcheck_noise.frac_empty` uses.

THE FINDING (specs/SPEC_certified_dipole_noise.md): raw coverage is broken on both HeH+ (healthy
Delta_lo margin) and LiH (fragile margin at M=16), and moderate inflation (z=1) restores >=0.9
coverage on HeH+ -- but at a tight shot budget (1e4), MORE inflation (z=3) leaves FEWER trials with
a finite bracket than z=1 (0.554 vs 0.728), an inflation CEILING neither prior noise spec showed
(both monotonic in z up to z=5-6). On LiH the margin is too thin at this depth for any tested z to
help at all (finite_frac stays < 0.3 throughout).

HONEST SCOPE: the coverage target is this module's own noiseless M=16 Krylov point estimate, not
the dense-diagonalization FCI dipole `certified_dipole` ultimately targets -- isolates the noise
question from the depth-convergence question those gates already own; i.i.d. Gaussian shot noise;
self-mode (oracle-free) eps1; two systems, full orbital space (`build_dipole_operators`'s scope).
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from certified_dipole import spectral_width
from certified_noise import certified_half_width, hamiltonian_one_norms
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import mean_and_variance


def operator_one_norms(op: SparsePauliOp) -> Tuple[float, float]:
    """(lambda_A, lambda_{A^2}) -- Pauli 1-norms of a Hermitian operator, mirroring
    `certified_noise.hamiltonian_one_norms` (same construction, generic to any SparsePauliOp)."""
    lam = float(np.abs(op.coeffs).sum())
    lam2 = float(np.abs((op @ op).simplify().coeffs).sum())
    return lam, lam2


def dipole_noise_coverage(mh: MolecularHamiltonian, a_op: SparsePauliOp, m: int, shots: float,
                          z: float = 2.0, trials: int = 6000, seed: int = 0,
                          solver: Optional[QuantumKrylovSolver] = None) -> dict:
    """Monte-Carlo coverage of the noiseless M-dim Krylov dipole estimate by the certified
    property bracket under i.i.d. shot noise, at inflation ``z`` (z=0 is raw/unpadded). A padded
    trial with s=sigma_0/Delta_lo >= 1 is vacuous: scored as not-covered, tracked in
    ``finite_frac``.
    """
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    a_sparse = a_op.to_matrix(sparse=True)
    W = spectral_width(a_sparse)

    _, states = solver.eigenstates(m, n_states=2)
    psi0, psi1 = states[0], states[1]
    th0, var0 = mean_and_variance(H, psi0)
    th1, var1 = mean_and_variance(H, psi1)
    ax = a_sparse @ psi0
    mu_exact = float((psi0.conj() @ ax).real)
    h2a_exact = float((ax.conj() @ ax).real)

    lam_h, lam_h2 = hamiltonian_one_norms(mh)
    lam_a, lam_a2 = operator_one_norms(a_op)
    se_h, se_h2 = lam_h / np.sqrt(shots), lam_h2 / np.sqrt(shots)
    se_a, se_a2 = lam_a / np.sqrt(shots), lam_a2 / np.sqrt(shots)
    hw_h = certified_half_width(lam_h, shots, z)
    hw_a = certified_half_width(lam_a, shots, z)

    rng = np.random.default_rng(seed)
    covered = 0
    finite = 0
    for _ in range(trials):
        n_th0 = th0 + rng.normal(0.0, se_h)
        n_h20 = (var0 + th0 * th0) + rng.normal(0.0, se_h2)
        n_var0 = max(n_h20 - n_th0 * n_th0, 0.0)
        n_th1 = th1 + rng.normal(0.0, se_h)
        n_h21 = (var1 + th1 * th1) + rng.normal(0.0, se_h2)
        n_var1 = max(n_h21 - n_th1 * n_th1, 0.0)
        sig1 = float(np.sqrt(n_var1))
        eps1 = n_th1 - sig1
        dlo_raw = eps1 - n_th0
        sigma0_raw = float(np.sqrt(n_var0))
        dlo_padded = dlo_raw - hw_h                 # conservative: shrink the margin
        sigma0_padded = sigma0_raw + hw_h            # conservative: grow the residual

        n_mu = mu_exact + rng.normal(0.0, se_a)
        n_h2a = h2a_exact + rng.normal(0.0, se_a2)
        n_sigA = float(np.sqrt(max(n_h2a - n_mu * n_mu, 0.0)))
        sigA_padded = n_sigA + hw_a

        if dlo_padded <= 0:
            continue
        s = sigma0_padded / dlo_padded
        if s >= 1.0:
            continue
        half = 2.0 * sigA_padded * s + W * s * s + hw_a
        finite += 1
        if (n_mu - half) <= mu_exact <= (n_mu + half):
            covered += 1

    return {
        "coverage": covered / trials, "finite_frac": finite / trials,
        "mu_exact": mu_exact, "lam_h": lam_h, "lam_a": lam_a, "spectral_width": W,
    }


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_dipole_operators, build_molecular_hamiltonian

    cases = {
        "HeH+": dict(atom="He 0 0 0; H 0 0 0.772", charge=1),
        "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    }
    print("=" * 78)
    print("Certified dipole under shot noise: coverage / finite-bracket rate vs inflation z")
    print("  system | shots   | z   | coverage | finite_frac")
    for name, spec in cases.items():
        mh = build_molecular_hamiltonian(**spec)
        az = build_dipole_operators(**spec)[2]
        for shots in (1e4, 1e5, 1e6):
            for z in (0.0, 1.0, 2.0, 3.0):
                r = dipole_noise_coverage(mh, az, 16, shots, z=z)
                print(f"  {name:5s} | {shots:.0e} | {z:.1f} |  {r['coverage']:.3f}   |  "
                      f"{r['finite_frac']:.3f}")
    print("=" * 78)
    print("HeH+: z=1 restores coverage but z=3 shrinks the finite-bracket rate below z=1's --")
    print("the inflation ceiling. LiH: finite_frac stays low regardless of z (thin margin at M=16).")
