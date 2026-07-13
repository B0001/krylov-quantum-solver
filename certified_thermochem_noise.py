#!/usr/bin/env python3
"""
The certified relative-energy bracket under shot noise -- composition beats either endpoint alone.

`certified_thermochem` composes two exact-statevector Temple/Ritz brackets (one per geometry) into
a rigorous interval on Delta = E(B) - E(A). `certified_noise` showed a SINGLE such bracket's raw
coverage of the true energy collapses to ~0.40 under i.i.d. shot noise -- a coin flip, not a bound,
because the Ritz value sits so close to E0 that symmetric noise pushes it below E0 half the time.
This module asks the composed question: does the collapse survive composition, and does fixing it
still cost the same z=2 inflation `certified_noise` needed for one bracket?

Each geometry is noised independently exactly as `certified_noise.shot_noise_coverage` does (own
lambda_H/lambda_H2 1-norms, own oracle E1), then composed the same way `certified_thermochem` does
noiselessly: Delta_lower = tau_B - rho_A, Delta_upper = rho_B - tau_A. Each geometry's own
`energy_offset` is added before composing -- geometries have different nuclear-repulsion offsets,
which do NOT cancel in a cross-geometry difference (unlike a single-geometry gap, where they do).

THE FINDING (specs/SPEC_certified_thermochem_noise.md): raw coverage of the composed bracket is
broken (~0.70, vs ~0.40 for either endpoint alone) and shot-count-independent, same mechanism as
`certified_noise` -- BUT the composed bracket needs markedly LESS inflation to fix: z* ~ 0.5-0.6
restores 90% coverage here, vs z=2 for a single bracket. Two independent one-sided coin-flips,
composed into a two-sided difference, are not as broken as either coin-flip alone.

HONEST SCOPE: oracle E1 per geometry (isolates the noise-composition question from the self-mode
Temple-premise question `certified_gaps`/`gap_selfcheck` own); i.i.d. Gaussian shot noise with
lambda-1-norm standard errors (same idealization as `certified_noise`); two composed geometries
only; the reduced z* is a measurement on this system (R1 in the spec), not a derived constant.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from certified_noise import hamiltonian_one_norms, reachable_E0_E1
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import _mean_and_variance


def _endpoint_stats(mh: MolecularHamiltonian, m: int):
    """(th0, h2_exact) -- noiseless Ritz mean and <H^2> of the ground Krylov state (electronic
    frame, matching `certified_noise`'s convention: offset added later, at composition time)."""
    solver = QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    psi0 = solver.eigenstates(m, n_states=1)[1][0]
    th0, var0 = _mean_and_variance(H, psi0)
    return th0, var0 + th0 * th0


def _noisy_rho_tau(th0: float, h2_exact: float, eps1_exact: float, se_h: float, se_h2: float,
                   z: float, n: int, rng: np.random.Generator):
    """One geometry's noisy (rho, tau) pair at inflation z, matching `certified_noise`'s per-sample
    model: rho = th (inflated by z*se_h), tau = th - var/(eps1-th) at th=th (uninflated point,
    inflated separately by z*se_h on the low side), var from the inflated h2 -- same construction as
    `certified_noise.shot_noise_coverage`'s `tau_i`/`th_lo`/`th_hi`."""
    th = th0 + rng.normal(0.0, se_h, n)
    h2 = h2_exact + rng.normal(0.0, se_h2, n)
    th_lo, th_hi = th - z * se_h, th + z * se_h
    var_hi = np.clip(h2 + z * se_h2 - th_lo * th_lo, 0.0, None)
    denom = eps1_exact - th_lo
    tau = np.where(denom > 0, th_lo - var_hi / np.where(denom > 0, denom, 1.0), -np.inf)
    return th_hi, tau  # rho (upper), tau (lower)


def thermochem_noise_coverage(mh_a: MolecularHamiltonian, mh_b: MolecularHamiltonian, m: int,
                              shots: float, e1_a: Optional[float] = None,
                              e1_b: Optional[float] = None, trials: int = 20000,
                              z: float = 2.0, seed: int = 0) -> dict:
    """Monte-Carlo coverage of the exact relative energy Delta = E(B)-E(A) by the composed
    certified bracket under i.i.d. shot noise, at inflation ``z`` (z=0 is the raw/uninflated
    bracket). ``e1_a``/``e1_b``: oracle E1 per geometry (electronic frame from
    `certified_noise.reachable_E0_E1`), defaulting to the exact value if not supplied.
    """
    offset_a, offset_b = mh_a.energy_offset, mh_b.energy_offset
    delta_exact = mh_b.ground_state_energy() - mh_a.ground_state_energy()

    _, e1_a_exact = reachable_E0_E1(mh_a)
    _, e1_b_exact = reachable_E0_E1(mh_b)
    e1_a = e1_a_exact if e1_a is None else e1_a
    e1_b = e1_b_exact if e1_b is None else e1_b

    th0_a, h2_a = _endpoint_stats(mh_a, m)
    th0_b, h2_b = _endpoint_stats(mh_b, m)
    lam_h_a, lam_h2_a = hamiltonian_one_norms(mh_a)
    lam_h_b, lam_h2_b = hamiltonian_one_norms(mh_b)
    se_h_a, se_h2_a = lam_h_a / np.sqrt(shots), lam_h2_a / np.sqrt(shots)
    se_h_b, se_h2_b = lam_h_b / np.sqrt(shots), lam_h2_b / np.sqrt(shots)

    rng = np.random.default_rng(seed)
    rho_a, tau_a = _noisy_rho_tau(th0_a, h2_a, e1_a, se_h_a, se_h2_a, z, trials, rng)
    rho_b, tau_b = _noisy_rho_tau(th0_b, h2_b, e1_b, se_h_b, se_h2_b, z, trials, rng)
    rho_a, tau_a = rho_a + offset_a, tau_a + offset_a
    rho_b, tau_b = rho_b + offset_b, tau_b + offset_b

    lo, hi = tau_b - rho_a, rho_b - tau_a
    cov = float(np.mean((lo <= delta_exact) & (delta_exact <= hi)))
    key = "cov_raw" if z == 0.0 else "cov_inflated"
    return {
        key: cov, "delta_exact": delta_exact,
        "lam_h_a": lam_h_a, "lam_h_b": lam_h_b,
        "lam_h2_a": lam_h2_a, "lam_h2_b": lam_h2_b,
    }


def minimal_z_for_coverage(mh_a: MolecularHamiltonian, mh_b: MolecularHamiltonian, m: int,
                           shots: float, target: float = 0.9, resolution: float = 0.05,
                           z_max: float = 3.0, **kwargs) -> float:
    """Smallest ``z`` (grid search at ``resolution``, up to ``z_max``) whose composed bracket
    reaches ``target`` coverage. Returns ``z_max`` if the target is never reached (a boundary, not
    a silent pass -- callers should treat that as the composition-helps claim failing there)."""
    z = 0.0
    while z <= z_max:
        cov = thermochem_noise_coverage(mh_a, mh_b, m, shots, z=z, **kwargs)
        cov_val = cov["cov_raw"] if z == 0.0 else cov["cov_inflated"]
        if cov_val >= target:
            return z
        z += resolution
    return z_max


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    def chain(r):
        return f"H 0 0 0; H 0 0 {r}; H 0 0 {2 * r}; H 0 0 {3 * r}"

    mh_eq = build_molecular_hamiltonian(atom=chain(0.9))
    mh_st = build_molecular_hamiltonian(atom=chain(2.3))
    print("=" * 78)
    print("H4 symmetric stretch: composed relative-energy bracket under shot noise (M=16)")
    print("  shots   | raw cov (z=0) | z=2 cov | minimal z for 90% coverage")
    for shots in (1e4, 1e5, 1e6):
        raw = thermochem_noise_coverage(mh_eq, mh_st, 16, shots, z=0.0)["cov_raw"]
        z2 = thermochem_noise_coverage(mh_eq, mh_st, 16, shots, z=2.0)["cov_inflated"]
        zstar = minimal_z_for_coverage(mh_eq, mh_st, 16, shots)
        print(f"  {shots:.0e} |     {raw:.3f}     |  {z2:.3f}  |  {zstar:.2f}")
    print("=" * 78)
    print("Composition partially self-corrects the certified_noise coin-flip collapse:")
    print("less inflation needed than the single-bracket z=2 rule.")
