#!/usr/bin/env python3
"""
Gap self-check under shot noise -- intersection concentrates noise where composition diluted it.

`gap_selfcheck` builds an oracle-free trustworthiness certificate by INTERSECTING several
independently-noisy-at-depth brackets (corroborated iff a bracket overlaps the deep anchor).
`certified_thermochem_noise` asked the noise question of the OTHER composition operator this repo
uses -- DIFFERENCE (Delta = E_B - E_A) -- and found composition partially self-corrects the
certified_noise coin-flip collapse (needs LESS inflation than a single bracket). This module asks
the same question of INTERSECTION, and finds the opposite: intersecting brackets needs MORE
inflation than a single bracket, not less -- composing certified intervals by AND concentrates
noise-induced miscalibration exactly where composing by difference diluted it.

Each depth in the ladder gets ONE noisy realization of (theta0, var0, theta1, sigma1), built with
the exact `certified_gaps.gap_bracket` formula (self-mode: eps1 = theta1 - sigma1) -- no internal
directional "worst-case" substitution, which the spec's probe found is genuinely ambiguous for a
two-eigenstate gap bracket (theta0 plays opposing conservative roles in gap_lower vs. the Temple
term feeding gap_upper). Inflation is applied the simplest, least ambiguous way instead: pad each
depth's already-computed `[gap_lower, gap_upper]` post-hoc by
`+/- certified_noise.certified_half_width(lambda_H, shots, z)`, then hand the padded ladder to
`gap_selfcheck.self_checked_gap` completely unmodified.

THE FINDING (specs/SPEC_gap_selfcheck_noise.md): raw (z=0) coverage of the self-checked interval is
broken (~0.15-0.20 H4, ~0.06-0.08 LiH) and shot-count-independent, same mechanism as
`certified_noise`. But the z=2 single-bracket rule -- already MORE than sufficient for
`certified_thermochem_noise`'s composed bracket -- does NOT restore 90% coverage here (0.70 H4,
0.53 LiH at shots=1e5); the minimal z that does is 3.25 (H4) / 4.00 (LiH), roughly 1.6-2x the
single-bracket rule. Padding does fix "inconclusive" (the self-checked interval coming back empty
because no bracket corroborates) well before it fixes "correct": frac_empty drops from ~0.1-0.13 at
z=0 to < 0.005 at z=2.0, even though coverage is still broken there.

HONEST SCOPE: self-mode (oracle-free) eps1 only; i.i.d. Gaussian shot noise with lambda-1-norm
standard errors (same idealization as `certified_noise`); two systems (H4, LiH) -- the specific z*
numbers are a measurement on these systems (R1 in the spec), the DIRECTION (more, not less,
inflation needed) is the falsifiable claim; no repair for the miscalibration is attempted (a
follow-up spec).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from certified_gaps import GapBracket
from certified_noise import certified_half_width, hamiltonian_one_norms
from gap_selfcheck import self_checked_gap
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import mean_and_variance


def _ladder_stats(mh: MolecularHamiltonian, dims: Sequence[int]) -> dict:
    """Noiseless (theta0, var0, theta1, var1) at each Krylov dimension (one solver, reused)."""
    solver = QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    stats = {}
    for m in dims:
        _, states = solver.eigenstates(m, n_states=2)
        th0, var0 = mean_and_variance(H, states[0])
        th1, var1 = mean_and_variance(H, states[1])
        stats[m] = (th0, var0, th1, var1)
    return stats


def _raw_noisy_ladder(stats: dict, dims: Sequence[int], se_h: float, se_h2: float,
                      rng: np.random.Generator) -> list[tuple[int, float, float]]:
    """ONE noisy realization per depth of (gap_lower, gap_upper), exactly the self-mode
    `certified_gaps.gap_bracket` formula with a single noisy (th0, var0, th1, var1) -- no
    directional inflation surgery (see module docstring: ambiguous for this two-eigenstate case)."""
    out = []
    for m in dims:
        th0, var0, th1, var1 = stats[m]
        h20, h21 = var0 + th0 * th0, var1 + th1 * th1
        n_th0 = th0 + rng.normal(0.0, se_h)
        n_h20 = h20 + rng.normal(0.0, se_h2)
        n_th1 = th1 + rng.normal(0.0, se_h)
        n_h21 = h21 + rng.normal(0.0, se_h2)
        n_var0 = max(n_h20 - n_th0 * n_th0, 0.0)
        n_var1 = max(n_h21 - n_th1 * n_th1, 0.0)
        sig1 = float(np.sqrt(n_var1))
        eps1 = n_th1 - sig1
        denom = eps1 - n_th0
        tau0 = n_th0 - n_var0 / denom if denom > 0 else -np.inf
        gap_upper = n_th1 - tau0 if np.isfinite(tau0) else np.inf
        gap_lower = eps1 - n_th0
        out.append((m, gap_lower, gap_upper))
    return out


def _pad(raw: list[tuple[int, float, float]], half_width: float) -> list[GapBracket]:
    """Post-hoc symmetric pad on each depth's bracket -- the simplest, least ambiguous inflation."""
    return [
        GapBracket(m=m, gap_lower=lo - half_width, gap_upper=hi + half_width,
                  width=(hi - lo) + 2 * half_width, theta0=0.0, theta1=0.0, sigma1=0.0,
                  eps1=0.0, eps1_source="self")
        for m, lo, hi in raw
    ]


def self_check_noise_coverage(mh: MolecularHamiltonian, dims: Sequence[int], shots: float,
                              z: float = 2.0, trials: int = 4000, seed: int = 0,
                              k: int = 2) -> dict:
    """Monte-Carlo coverage of the exact reachable gap by the oracle-free self-checked interval
    under i.i.d. shot noise, at inflation ``z`` (z=0 is the raw/unpadded ladder). Independent noise
    per depth (different measurements at different Krylov dimensions).
    """
    from certified_gaps import reachable_gap

    gap_exact = reachable_gap(mh)
    lam_h, lam_h2 = hamiltonian_one_norms(mh)
    se_h, se_h2 = lam_h / np.sqrt(shots), lam_h2 / np.sqrt(shots)
    half_width = certified_half_width(lam_h, shots, z)
    stats = _ladder_stats(mh, dims)

    rng = np.random.default_rng(seed)
    covered = 0
    empty = 0
    for _ in range(trials):
        raw = _raw_noisy_ladder(stats, dims, se_h, se_h2, rng)
        brackets = _pad(raw, half_width)
        try:
            lo, hi = self_checked_gap(brackets, k=k)
        except ValueError:
            empty += 1
            continue
        if lo <= gap_exact <= hi:
            covered += 1
    return {
        "coverage": covered / trials, "frac_empty": empty / trials,
        "lam_h": lam_h, "lam_h2": lam_h2, "gap_exact": gap_exact,
    }


def minimal_z_for_selfcheck_coverage(mh: MolecularHamiltonian, dims: Sequence[int], shots: float,
                                     target: float = 0.9, resolution: float = 0.25,
                                     z_max: float = 6.0, **kwargs) -> float:
    """Smallest ``z`` (grid search at ``resolution``, up to ``z_max``) whose self-checked interval
    reaches ``target`` coverage. Returns ``z_max`` if never reached (a recorded boundary, not a
    silent pass)."""
    z = 0.0
    while z <= z_max:
        cov = self_check_noise_coverage(mh, dims, shots, z=z, **kwargs)["coverage"]
        if cov >= target:
            return z
        z += resolution
    return z_max


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
        "LiH": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
    }
    dims = (6, 8, 12, 16, 20, 24)
    print("=" * 78)
    print("Self-checked gap under shot noise: raw coverage, z=2 (single-bracket rule), minimal z")
    print("  system | shots   | raw cov | z=2 cov | minimal z for 90%")
    for name, spec in cases.items():
        mh = build_molecular_hamiltonian(**spec)
        for shots in (1e4, 1e5, 1e6):
            raw = self_check_noise_coverage(mh, dims, shots, z=0.0)["coverage"]
            z2 = self_check_noise_coverage(mh, dims, shots, z=2.0)["coverage"]
            print(f"  {name:5s} | {shots:.0e} |  {raw:.3f}  |  {z2:.3f}  |", end=" ")
            if shots == 1e5:
                zstar = minimal_z_for_selfcheck_coverage(mh, dims, shots)
                print(f"{zstar:.2f}")
            else:
                print("(only computed at 1e5)")
    print("=" * 78)
    print("Intersection composition needs MORE inflation than a single bracket -- the opposite")
    print("of certified_thermochem_noise's difference composition, which needed less.")
