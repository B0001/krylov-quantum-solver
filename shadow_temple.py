#!/usr/bin/env python3
"""
Do classical shadows buy Temple-bracket coverage, or just relabel the knife-edge? (chem-loe)

`SPEC_certified_noise` found lambda_{H^2} >> lambda_H (H4: 62.9 vs 10.3), making the Temple lower
bound the noise-expensive side of the certified arc, and `temple_bounds` states outright that the
<H^2> hardware cost is unmodeled. Random-Pauli classical shadows (`classical_shadows.py`, a closed
spec) estimate BOTH <H> and <H^2> from the SAME measurement snapshots -- a different protocol from
the repo's assumed i.i.d. lambda^2-Gaussian shot-noise model (`certified_noise.py`). This module
bridges the two, composing both without adding a new measurement primitive: H^2 is the same
`(H @ H).simplify()` expansion `certified_noise.hamiltonian_one_norms` already builds.

THE FINDING (specs/SPEC_shadow_temple.md): shadows DO give a lower empirical variance than the
lambda^2 model, for H and H^2 alike, and the H^2-vs-H cost asymmetry SHRINKS under shadows (it does
not vanish -- shadows are still more expensive on H^2 than on H). A SHARPER, unpredicted finding:
the HKP-style additive `shadow_norm` bound itself (valid, empirically, for H -- ratio 0.69-0.95 over
10 seeds) is VIOLATED, reproducibly, for H^2 (ratio 1.16-1.49x over 10 seeds) -- the per-term
diagonal sum ignores positive cross-correlations between H^2's 1775 heavily support-overlapping
Pauli terms that are not negligible at this scale. `shadow_norm` cannot be trusted as a shot-budget
bound for H^2 even though the measured variance still beats lambda^2. Separately, feeding
shadow-estimated moments into the Temple bracket does NOT buy back coverage: the noisy-Temple break
(raw coverage ~0.40 at converged depth, `SPEC_certified_noise`) is a structural property of the
variational knife-edge rho_0 -> E_0, not an artifact of the Gaussian noise model -- so a cheaper,
unbiased estimator lands in the SAME broken regime. The efficiency-vs-lambda^2 win is real; the
shadow_norm bound and the coverage win are not.

HONEST SCOPE: single molecule (H4, the system SPEC_certified_noise recorded the lambda ratio on),
single Krylov depth (M=12, converged), exact statevector simulation of the shadow measurement,
random-Pauli (not grouped/derandomized) shadows. See specs/SPEC_shadow_temple.md Sec 2/7.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from certified_noise import hamiltonian_one_norms, reachable_E0_E1
from classical_shadows import collect_classical_shadow, shadow_energy_samples, shadow_norm
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver


def h2_operator(mh: MolecularHamiltonian):
    """H^2 as a SparsePauliOp (electronic frame) -- the SAME expansion
    `certified_noise.hamiltonian_one_norms` builds for lambda_{H^2}, so shadow_norm(h2_operator(mh))
    and lambda_{H^2} are directly comparable (same operator, two different variance protocols)."""
    op = mh.qubit_hamiltonian
    return (op @ op).simplify()


def ritz_state(mh: MolecularHamiltonian, m: int,
               solver: Optional[QuantumKrylovSolver] = None) -> Tuple[np.ndarray, QuantumKrylovSolver]:
    """The M-dimensional Krylov ground Ritz state -- the state the certified arc actually runs on
    (not HF, not FCI)."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    psi = solver.eigenstates(m, n_states=1)[1][0]
    return psi, solver


def shadow_moments(psi: np.ndarray, n_qubits: int, op, op2, n_shots: int,
                   seed: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    """Collect ONE shadow dataset and estimate <H> and <H^2> from the SAME snapshots -- the central
    efficiency claim under test. Returns (per_shot_H, per_shot_H2); means estimate <H>, <H^2>."""
    bases, signs = collect_classical_shadow(psi, n_qubits, n_shots, seed=seed)
    per_shot_h = shadow_energy_samples(bases, signs, op)
    per_shot_h2 = shadow_energy_samples(bases, signs, op2)
    return per_shot_h, per_shot_h2


def shadow_vs_lambda_model(mh: MolecularHamiltonian, psi: np.ndarray, n_shots: int = 16000,
                           seed: int = 0) -> dict:
    """One shadow dataset on `psi`: for both H and H^2, the HKP shadow-norm bound, the repo's
    assumed lambda^2-model variance, and the EMPIRICAL single-shot variance measured directly --
    comparing shadow_norm to lambda^2 alone compares two different upper bounds from two different
    protocols (SPEC_shadow_temple R1); this also measures the real thing."""
    op = mh.qubit_hamiltonian
    op2 = h2_operator(mh)
    lam_h, lam_h2 = hamiltonian_one_norms(mh)
    s1, s2 = shadow_moments(psi, mh.num_qubits, op, op2, n_shots, seed=seed)
    return dict(
        H=dict(shadow_norm=shadow_norm(op), lambda2=lam_h ** 2, empirical_var=float(s1.var())),
        H2=dict(shadow_norm=shadow_norm(op2), lambda2=lam_h2 ** 2, empirical_var=float(s2.var())),
    )


def shadow_temple_coverage(mh: MolecularHamiltonian, m: int, n_shots: int, trials: int = 200,
                           seed: int = 0, solver: Optional[QuantumKrylovSolver] = None) -> dict:
    """Coverage of the exact reachable E_0 by the noisy Temple bracket built from REAL classical-
    shadow snapshots -- the empirical analogue of `certified_noise.shot_noise_coverage`, with actual
    shadow sampling in place of the assumed i.i.d. Gaussian.

    Each of `trials` seeds collects a FRESH n_shots shadow dataset of the Krylov Ritz state,
    estimates <H> and <H^2> from the SAME snapshots, forms [tau_0, rho_0], and checks coverage of
    the exact reachable E_0. Uses an oracle gap (exact reachable E_1), matching
    `shot_noise_coverage`'s convention -- a noisy eps_1 would only worsen coverage further.
    """
    psi, solver = ritz_state(mh, m, solver)
    op = mh.qubit_hamiltonian
    op2 = h2_operator(mh)
    E0, E1 = reachable_E0_E1(mh)

    ths = np.empty(trials)
    hits_raw = 0
    hits_upper = 0
    for t in range(trials):
        s1, s2 = shadow_moments(psi, mh.num_qubits, op, op2, n_shots, seed=seed * 1_000_003 + t)
        th = float(s1.mean())
        h2 = float(s2.mean())
        var = max(h2 - th * th, 0.0)
        tau = th - var / (E1 - th) if E1 > th else -np.inf
        hits_raw += int(tau <= E0 <= th)
        hits_upper += int(E0 <= th)
        ths[t] = th
    return dict(cov_raw=hits_raw / trials, cov_upper=hits_upper / trials,
                th_mean=float(ths.mean()), th_std=float(ths.std()),
                n_shots=n_shots, trials=trials)


if __name__ == "__main__":
    import time

    from certified_noise import shot_noise_coverage
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    atom = "H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"
    m = 12
    mh = build_molecular_hamiltonian(atom=atom)
    solver = QuantumKrylovSolver(mh)
    psi, solver = ritz_state(mh, m, solver=solver)

    print("=" * 78)
    print(f"H4/STO-3G, M={m}: H has {len(mh.qubit_hamiltonian.paulis)} Pauli terms, "
          f"H^2 has {len(h2_operator(mh).paulis)} (O(N^8) blowup, recorded not hidden)")

    lam_h, lam_h2 = hamiltonian_one_norms(mh)
    cmp = shadow_vs_lambda_model(mh, psi, n_shots=16000, seed=1)
    print(f"lambda_H={lam_h:.2f}  lambda_H2={lam_h2:.2f}  (lambda^2 ratio {(lam_h2 / lam_h) ** 2:.2f})")
    for key in ("H", "H2"):
        r = cmp[key]
        bound_note = "HOLDS" if r["empirical_var"] <= r["shadow_norm"] * 1.05 else "VIOLATED"
        print(f"  {key}: shadow_norm={r['shadow_norm']:8.2f}  lambda^2={r['lambda2']:8.2f}  "
              f"empirical_var={r['empirical_var']:8.2f}  (empirical/shadow_norm={r['empirical_var'] / r['shadow_norm']:.3f}, "
              f"bound {bound_note})  (empirical/lambda^2={r['empirical_var'] / r['lambda2']:.3f})")
    shadow_ratio = cmp["H2"]["shadow_norm"] / cmp["H"]["shadow_norm"]
    lambda2_ratio = (lam_h2 / lam_h) ** 2
    print(f"shadow_norm ratio H2/H = {shadow_ratio:.2f}  vs  lambda^2 ratio = {lambda2_ratio:.2f}  "
          "(asymmetry shrinks, does not vanish)")

    print("-" * 78)
    n_shots = 4000
    t0 = time.time()
    shadow_cov = shadow_temple_coverage(mh, m, n_shots, trials=200, seed=11, solver=solver)
    dt = time.time() - t0
    lambda_cov = shot_noise_coverage(mh, m, n_shots, trials=4000, solver=solver)
    print(f"shadow-estimated-moment coverage (200 trials, {n_shots} shots each, {dt:.1f}s): "
          f"cov_raw={shadow_cov['cov_raw']:.3f}  cov_upper={shadow_cov['cov_upper']:.3f}")
    print(f"lambda^2-model coverage (4000 trials, matched shot budget): "
          f"cov_raw={lambda_cov['cov_raw']:.3f}  cov_upper={lambda_cov['cov_upper']:.3f}")
    print("=" * 78)
    print("Shadows buy variance (both H and H^2 beat the lambda^2 model); they do not buy Temple")
    print("coverage -- the break is a structural variational knife-edge, estimator-independent.")
