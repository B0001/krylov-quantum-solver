"""
Acceptance gates G1-G4 for specs/SPEC_shadow_temple.md (chem-loe).

Claim: random-Pauli classical shadows estimate <H> and <H^2> from the SAME measurement snapshots,
with a variance (bound AND measured empirically) below the repo's assumed lambda^2 shot-noise
model -- but feeding shadow-estimated moments into the Temple bracket does NOT buy back coverage:
the noisy-Temple break SPEC_certified_noise found (raw coverage ~0.40) is a structural variational
knife-edge (rho_0 -> E_0), not an artifact of the lambda^2-Gaussian noise model, so shadows just
relabel it rather than removing it.

H4/STO-3G (the SPEC_certified_noise case with the recorded lambda_H2 >> lambda_H boundary; 8
qubits, H: 185 / H^2: 1775 Pauli terms). Seeded (deterministic). PySCF/qiskit only, no block2;
`make gates` runs it in its own process.
"""
import numpy as np

from certified_noise import hamiltonian_one_norms, reachable_E0_E1, shot_noise_coverage
from classical_shadows import shadow_norm
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from shadow_temple import (
    h2_operator,
    ritz_state,
    shadow_moments,
    shadow_temple_coverage,
    shadow_vs_lambda_model,
)
from temple_bounds import mean_and_variance

_H4 = "H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"
_M = 12


def _mh_solver():
    mh = build_molecular_hamiltonian(atom=_H4)
    return mh, QuantumKrylovSolver(mh)


def test_G1_unbiased_on_krylov_ritz_state():
    """The shadow mean of BOTH <H> and <H^2> is within 4*stderr of the exact value on the M=12
    Krylov Ritz state (not HF, not FCI -- the state the certified arc actually runs on)."""
    mh, solver = _mh_solver()
    psi, solver = ritz_state(mh, _M, solver=solver)
    op = mh.qubit_hamiltonian
    op2 = h2_operator(mh)
    H = op.to_matrix(sparse=True).tocsc()
    th_exact, var_exact = mean_and_variance(H, psi)
    h2_exact = var_exact + th_exact ** 2

    s1, s2 = shadow_moments(psi, mh.num_qubits, op, op2, 16000, seed=3)
    se1 = s1.std() / np.sqrt(s1.size)
    se2 = s2.std() / np.sqrt(s2.size)
    assert abs(s1.mean() - th_exact) < 4.0 * se1, (s1.mean(), th_exact, se1)
    assert abs(s2.mean() - h2_exact) < 4.0 * se2, (s2.mean(), h2_exact, se2)


def test_G2_shadow_bound_beats_lambda_model_but_keeps_the_asymmetry():
    """shadow_norm < lambda^2 for BOTH H and H^2 (an efficiency win); the shadow-norm H^2/H ratio is
    smaller than the lambda-based ratio (the asymmetry SHRINKS) but stays > 1 (it does not vanish)."""
    mh, _ = _mh_solver()
    op = mh.qubit_hamiltonian
    op2 = h2_operator(mh)
    lam_h, lam_h2 = hamiltonian_one_norms(mh)
    norm_h, norm_h2 = shadow_norm(op), shadow_norm(op2)

    assert norm_h < lam_h ** 2, (norm_h, lam_h ** 2)
    assert norm_h2 < lam_h2 ** 2, (norm_h2, lam_h2 ** 2)

    shadow_ratio = norm_h2 / norm_h
    lambda_ratio = lam_h2 / lam_h
    assert shadow_ratio < lambda_ratio, (shadow_ratio, lambda_ratio)   # asymmetry SHRINKS
    assert shadow_ratio > 1.0, shadow_ratio                            # but does NOT vanish


def test_G3_empirical_variance_measured_directly():
    """The empirical single-shot variance (measured, not inferred from bounds alone) is <=
    shadow_norm (HKP bound holds) AND < lambda^2 (beats the model it is compared against) for both
    H and H^2 -- the apples-to-apples check the bound-only comparison in G2 cannot provide alone."""
    mh, solver = _mh_solver()
    psi, solver = ritz_state(mh, _M, solver=solver)
    result = shadow_vs_lambda_model(mh, psi, n_shots=16000, seed=5)

    for key in ("H", "H2"):
        r = result[key]
        assert 0 < r["empirical_var"] <= r["shadow_norm"] * 1.05, (key, r)   # HKP bound (5% seed slack)
        assert r["empirical_var"] < r["lambda2"], (key, r)                   # beats the lambda^2 model


def test_G4_shadow_moments_do_not_buy_temple_coverage():
    """THE REAL TEST: 200-trial Temple coverage of exact E_0 using shadow-estimated moments lands in
    the SAME broken regime (~0.40) as the lambda^2-Gaussian model -- shadows relabel the knife-edge,
    they do not remove it. Falsified if shadow coverage is materially (>0.15) higher than the
    lambda^2-model coverage at a matched shot budget (a real free lunch)."""
    mh, solver = _mh_solver()
    n_shots = 4000
    shadow_result = shadow_temple_coverage(mh, _M, n_shots, trials=200, seed=11, solver=solver)
    lambda_result = shot_noise_coverage(mh, _M, n_shots, trials=4000, solver=solver)

    assert 0.0 <= shadow_result["cov_raw"] <= 1.0
    assert shadow_result["trials"] == 200
    # The structural break: shadow-based coverage stays low (same knife-edge regime as the
    # lambda^2 model), not near-perfect.
    assert shadow_result["cov_raw"] < 0.65, shadow_result
    # THE FALSIFIABLE CLAIM: shadows do not buy back coverage relative to the lambda^2 model.
    assert shadow_result["cov_raw"] - lambda_result["cov_raw"] < 0.15, (shadow_result, lambda_result)


def test_G4_reachable_E0_used_is_finite():
    """Sanity: the coverage target E_0/E_1 used by shadow_temple_coverage matches the same reachable
    sector certified_noise uses (shared reference, not a re-derivation)."""
    mh, _ = _mh_solver()
    E0, E1 = reachable_E0_E1(mh)
    assert np.isfinite(E0) and np.isfinite(E1) and E1 > E0
