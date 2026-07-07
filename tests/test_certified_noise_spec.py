"""
Acceptance gates G1-G4 for specs/SPEC_certified_noise.md.

Claim: under i.i.d. shot noise (standard errors set by the Hamiltonian 1-norms lambda_H, lambda_H2),
the certified energy bracket's guarantee BREAKS -- raw coverage of E_0 falls to ~0.4 and the
variational upper bound to a ~0.5 coin flip at converged depth -- and this is N-INDEPENDENT (shots do
not restore it). Inflating by z*standard-error restores coverage >= 0.9 (conservative), with the
inflated half-width scaling as z*lambda_H/sqrt(N): shots buy tightness, not coverage.

Idealized i.i.d. Gaussian shot noise, oracle gap, exact statevector Ritz state. PySCF/qiskit, no
block2; `make gates` runs it in its own process.
"""
from certified_noise import (
    certified_half_width,
    hamiltonian_one_norms,
    shot_noise_coverage,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

_CASES = {
    "H2": ("H 0 0 0; H 0 0 0.74", 6),
    "H4": ("H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7", 12),
}


def _mh_solver(atom):
    mh = build_molecular_hamiltonian(atom=atom)
    return mh, QuantumKrylovSolver(mh)


def test_G1_sampling_breaks_the_certificate():
    """At converged depth the noisy certificate is broken: raw two-sided coverage < 0.65 and the
    variational UPPER bound holds < 0.65 (a coin flip, not a bound) -- symmetric noise on rho_0 -> E_0
    lands below E_0 about half the time."""
    for atom, m in _CASES.values():
        mh, solver = _mh_solver(atom)
        r = shot_noise_coverage(mh, m, 1e6, solver=solver)
        assert r["cov_raw"] < 0.65, (atom, r["cov_raw"])
        assert r["cov_upper"] < 0.65, (atom, r["cov_upper"])


def test_G2_shots_do_not_buy_coverage():
    """DEFINITION OF DONE (the surprise): raw coverage is ~N-INDEPENDENT -- spanning N from 1e4 to
    1e8 moves it by < 0.12. The broken guarantee is structural; you cannot shoot your way out."""
    mh, solver = _mh_solver(_CASES["H4"][0])
    m = _CASES["H4"][1]
    covs = [shot_noise_coverage(mh, m, N, solver=solver)["cov_raw"] for N in (1e4, 1e6, 1e8)]
    assert max(covs) - min(covs) < 0.12, covs
    assert all(c < 0.65 for c in covs), covs                        # and all broken


def test_G3_inflation_restores_coverage_and_width_scales():
    """Inflation by z*se restores coverage >= 0.9 (conservative) at every budget; and the inflated
    half-width scales as lambda_H/sqrt(N) -- a 100x shot increase shrinks it ~10x (1/sqrt(N))."""
    for atom, m in _CASES.values():
        mh, solver = _mh_solver(atom)
        for N in (1e4, 1e6, 1e8):
            assert shot_noise_coverage(mh, m, N, solver=solver)["cov_inflated"] >= 0.9, (atom, N)
    lam_h, _ = hamiltonian_one_norms(build_molecular_hamiltonian(atom=_CASES["H4"][0]))
    w4 = certified_half_width(lam_h, 1e4)
    w6 = certified_half_width(lam_h, 1e6)
    assert abs(w4 / w6 - 10.0) < 0.2, (w4, w6)                       # 1/sqrt(N): sqrt(100)=10


def test_G4_H2_side_is_the_noise_expensive_boundary():
    """The honest boundary: <H^2> carries the larger 1-norm (lambda_H2 > lambda_H), so the Temple
    lower bound is the noise-expensive side; and coverage numbers are deterministic (fixed seed)."""
    for atom, _ in _CASES.values():
        lam_h, lam_h2 = hamiltonian_one_norms(build_molecular_hamiltonian(atom=atom))
        assert lam_h2 > lam_h, (atom, lam_h, lam_h2)
    mh, solver = _mh_solver(_CASES["H2"][0])
    a = shot_noise_coverage(mh, 6, 1e6, solver=solver, seed=0)["cov_raw"]
    b = shot_noise_coverage(mh, 6, 1e6, solver=solver, seed=0)["cov_raw"]
    assert a == b                                                    # deterministic under fixed seed
