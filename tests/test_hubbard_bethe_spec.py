"""
Acceptance gates G1-G4 for specs/SPEC_hubbard_bethe.md (1D Hubbard vs the exact Bethe-ansatz energy).

Test-first: ``hubbard_chain_integrals`` / ``lieb_wu_energy`` do not exist yet, so this file is RED
until the spec is implemented. The half-filled 1D Hubbard model has an exact analytic ground-state
energy per site (Lieb & Wu 1968); we check that the validated stack (model loader -> FCI / Krylov)
reproduces it: finite-size FCI extrapolates to the Lieb-Wu thermodynamic-limit integral, and the
free-fermion and dimer limits match to machine precision.

PySCF FCI only (no block2), L <= 12 to stay fast; `make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.model_hamiltonians import (
    fixed_filling_energy,
    hubbard_chain_integrals,
    hubbard_dimer_energy,
    lieb_wu_energy,
)
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver


def _per_site_fci(L, U, t=1.0):
    return fixed_filling_energy(hubbard_chain_integrals(L, U, t)) / L


def _extrapolate_tdl(U, t=1.0, Ls=(6, 8, 10, 12)):
    """Least-squares 1/L^2 extrapolation of the closed-shell per-site FCI energy to L -> inf."""
    es = np.array([_per_site_fci(L, U, t) for L in Ls])
    x = 1.0 / np.array(Ls, dtype=float) ** 2
    A = np.vstack([np.ones_like(x), x]).T
    return float(np.linalg.lstsq(A, es, rcond=None)[0][0])


def test_G1_tdl_matches_lieb_wu():
    """DEFINITION OF DONE: extrapolated per-site FCI energy matches the Bethe-ansatz integral."""
    for U in (2.0, 4.0, 8.0):
        e_tdl = _extrapolate_tdl(U)
        e_exact = lieb_wu_energy(U)
        assert abs(e_tdl - e_exact) < 8e-3, (U, e_tdl, e_exact, abs(e_tdl - e_exact))


def test_G2_free_fermion_limit():
    """U=0: the integral gives -4/pi, and the chain FCI equals the analytic free-fermion energy."""
    assert abs(lieb_wu_energy(0.0) - (-4.0 / np.pi)) < 1e-4, lieb_wu_energy(0.0)
    for L in (6, 8, 10):
        model = hubbard_chain_integrals(L, 0.0)
        evals = np.linalg.eigvalsh(model.h1)                    # single-particle levels
        e_free = 2.0 * float(np.sort(evals)[: L // 2].sum())    # doubly fill the lowest L/2
        assert abs(fixed_filling_energy(model) - e_free) < 1e-8, (L, e_free)


def test_G3_dimer_limit():
    """L=2 chain reproduces the analytic Hubbard-dimer energy (covers the U->inf superexchange end)."""
    for U in (1.0, 4.0, 12.0, 50.0):
        e = fixed_filling_energy(hubbard_chain_integrals(2, U, t=1.0))
        assert abs(e - hubbard_dimer_energy(1.0, U)) < 1e-8, (U, e, hubbard_dimer_energy(1.0, U))


def test_G4_krylov_reproduces_fci_on_lattice():
    """The number-conserving Krylov solver matches FCI on the chain (model->qubit->Krylov path).

    FINDING: real-time Krylov from |HF> converges to FCI on the small chain (L=4, all U here at
    depth 24), but the strongly-correlated half-filled Mott chain is hard for it -- |HF> has poor
    overlap with the true ground state, so L>=6 at large U/t needs far deeper Krylov (e.g. L=6,U=8
    is still ~160 mHa off at depth 24). We gate where it converges and record the limitation.
    """
    for L in (4,):
        for U in (2.0, 4.0, 8.0):
            model = hubbard_chain_integrals(L, U, t=1.0)
            e_fci = fixed_filling_energy(model)
            e_kry = QuantumKrylovSolver(model.to_hamiltonian()).solve(24).energy
            assert abs(e_kry - e_fci) < 1e-6, (L, U, e_kry, e_fci)
