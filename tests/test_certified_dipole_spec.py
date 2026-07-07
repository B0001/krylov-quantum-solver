"""
Acceptance gates G1-G4 for specs/SPEC_certified_dipole.md.

Claim: the ground-state dipole moment gets a certified interval [mu +/- half_width] from Krylov data
alone (no FCI), via the Davis-Kahan eigenvector bound sin theta <= sigma_0 / Delta_lo (Delta_lo the
certified GAP lower bound of certified_gaps) and the SHARP fluctuation bound half_width =
2 sigma_A s + W_A s^2. The exact FCI dipole lies inside at every depth (zero escapes); the interval
closes with depth; and it is finite iff s < 1 -- so the property certificate INHERITS the gap
certificate (vacuous where Delta_lo is weak).

Full orbital space, exact statevector, sector-restricted. PySCF/qiskit, no block2; `make gates` runs
it in its own process.
"""
import numpy as np
import pytest

from certified_dipole import certified_dipole_ladder, spectral_width
from hybrid_quantum_solver.molecular_hamiltonian import (
    build_dipole_operators,
    build_molecular_hamiltonian,
)
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

CASES = {
    "HeH+": dict(atom="He 0 0 0; H 0 0 0.772", charge=1),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
}
_DIMS = (6, 8, 12, 16, 20, 24)


@pytest.fixture(scope="module")
def ladders():
    out = {}
    for name, spec in CASES.items():
        mh = build_molecular_hamiltonian(**spec)
        Az = build_dipole_operators(**spec)[2].to_matrix(sparse=True)
        # reference = the HF-REACHABLE ground state (correct particle-number sector), NOT the global
        # lowest eigenvector (which for a charged species lives in a different sector).
        w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        hf = np.asarray(mh.hf_state().data, dtype=complex)
        reach = np.where(np.abs(V.conj().T @ hf) ** 2 > 1e-10)[0]
        psi_ex = V[:, reach[np.argmin(w[reach])]]
        mu_exact = float((psi_ex.conj() @ (Az @ psi_ex)).real)
        ladder = certified_dipole_ladder(mh, Az, _DIMS, solver=QuantumKrylovSolver(mh))
        out[name] = (mu_exact, Az, ladder)
    return out


def test_G1_zero_escapes(ladders):
    """The exact FCI dipole lies inside [mu +/- half_width] at EVERY depth (certification), and at
    the deepest depth the interval is FINITE and still contains it -- so it is not passing only via
    vacuous (infinite) intervals. One escape kills the claim."""
    for name, (mu_exact, _, ladder) in ladders.items():
        for cd in ladder:
            assert cd.mu - cd.half_width - 1e-9 <= mu_exact <= cd.mu + cd.half_width + 1e-9, (
                name, cd.m, cd.mu, cd.half_width, mu_exact)
        deepest = ladder[-1]
        assert deepest.finite, (name, "deepest interval must be non-vacuous")
        assert (deepest.mu - deepest.half_width - 1e-9
                <= mu_exact <= deepest.mu + deepest.half_width + 1e-9), (name, deepest)


def test_G2_interval_closes_and_is_useful(ladders):
    """DEFINITION OF DONE: at M=24 the certified dipole is a genuinely useful error bar -- LiH lands
    -1.818 +/- < 0.15 a.u. containing the exact -1.817; HeH+ is certified to < 1e-2 a.u. (its
    reachable subspace saturates, sigma_0 -> 0)."""
    lih = ladders["LiH"][2][-1]
    assert lih.finite and lih.half_width < 0.15, lih.half_width
    assert abs(ladders["LiH"][0] - (-1.8174)) < 1e-2                  # exact dipole sanity
    heh = ladders["HeH+"][2][-1]
    assert heh.finite and heh.half_width < 1e-2, heh.half_width


def test_G3_fluctuation_bound_beats_operator_norm(ladders):
    """The SHARP bound uses the dipole fluctuation sigma_A, not the operator norm ||A||: for LiH
    sigma_A << ||mu_z||, so the fluctuation half-width is many times tighter than the naive
    2 ||A|| sin theta bound would give (> 3x)."""
    mu_exact, Az, ladder = ladders["LiH"]
    cd = ladder[-1]                                                   # M=24, finite
    assert cd.finite
    a_dense = Az.toarray()
    op_norm = float(np.max(np.abs(np.linalg.eigvalsh(a_dense))))      # ||A||
    W = spectral_width(Az)
    naive_half = 2.0 * op_norm * cd.sin_theta_bound + W * cd.sin_theta_bound ** 2
    assert cd.half_width < 0.3 * naive_half, (cd.half_width, naive_half)
    assert cd.sigma_A < 0.25 * op_norm, (cd.sigma_A, op_norm)         # fluctuation << norm


def test_G4_property_inherits_gap_certificate(ladders):
    """THE BOUNDARY: half_width is finite iff s < 1 iff sigma_0 < Delta_lo, so the property
    certificate inherits the certified-gap lower bound. LiH exhibits BOTH regimes -- vacuous at the
    depths where Delta_lo is weak and sharp where it is healthy -- proving the dependence."""
    for name, (_, _, ladder) in ladders.items():
        for cd in ladder:
            assert cd.finite == (cd.sin_theta_bound < 1.0), (name, cd.m)
    lih = ladders["LiH"][2]
    assert any(not cd.finite for cd in lih), "expected some vacuous (weak Delta_lo) depths for LiH"
    assert any(cd.finite for cd in lih), "expected some sharp (healthy Delta_lo) depths for LiH"
    # where finite, Delta_lo is positive and exceeds the residual bound implied by s<1
    for cd in lih:
        if cd.finite:
            assert cd.gap_lower > 0.0, cd.m
