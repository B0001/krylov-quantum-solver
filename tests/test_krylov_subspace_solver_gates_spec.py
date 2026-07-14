"""
Acceptance gates G1-G4 for specs/SPEC_krylov_subspace_solver_gates.md (krylov_subspace_solver --
the documented bug fix gets a regression test, and the cross-check finally happens).

Deliberately no new library code: the checks live entirely here, reusing
`krylov_subspace_solver.py`'s existing `krylov_ground_state`/`krylov_convergence_sweep`/
`_dense_active_H` unmodified. G2's absolute-cutoff comparison is a deliberate REIMPLEMENTATION of
the module's former, buggy behavior (never imported from the module) -- the point is to prove the
fix matters, not to exercise the current code a second time.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from krylov_subspace_solver import _dense_active_H, krylov_convergence_sweep, krylov_ground_state

SYSTEMS = {
    "H2": ("H 0 0 0; H 0 0 0.74", 2, 2, 0),
    "H4": ("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4, 0),
    "O2_triplet": ("O 0 0 0; O 0 0 1.21", 4, 4, 2),
}


def _reference(label):
    atom, norb, ne, spin = SYSTEMS[label]
    mol = gto.M(atom=atom, basis="sto-3g", spin=spin)
    mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
    mf.verbose = 0
    mf.kernel()
    na, nb = (ne + spin) // 2, (ne - spin) // 2
    cas = mcscf.CASCI(mf, norb, (na, nb))
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return h1, eri, float(e_core), (na, nb), norb, float(cas.e_tot)


@pytest.mark.parametrize("label", SYSTEMS)
def test_G1_variational_and_convergent(label):
    """Every step of the convergence sweep stays >= CASCI (variational), and the deepest tested
    dimension is within 5 mHa -- pins the informal __main__ assertion into a real gate, including
    the open-shell O2 triplet case."""
    h1, eri, e_core, nelec, norb, casci = _reference(label)
    rows = krylov_convergence_sweep(h1, eri, e_core, nelec, norb, dims=(2, 4, 8, 12),
                                    casci_energy=casci)
    assert all(r["variational_ok"] for r in rows), (label, rows)
    assert rows[-1]["delta_mHa"] < 5.0, (label, rows[-1])


def _independent_S(h1, eri, e_core, nelec, norb, dt, krylov_dim, phi0_scale=1.0):
    """Reconstructs S exactly as `krylov_ground_state` does internally, but with the reference
    state scaled by `phi0_scale` -- a benign renormalization that must not change the physics."""
    Hd, dim = _dense_active_H(h1, eri, norb, nelec)
    w, V = np.linalg.eigh(Hd)
    phi0 = np.zeros(dim)
    phi0[0] = phi0_scale
    c0 = V.T @ phi0
    times = dt * np.arange(krylov_dim)
    phase = np.exp(-1j * np.outer(times, w))
    B = (V @ (phase * c0).T).T
    S = B.conj() @ B.T
    return 0.5 * (S + S.conj().T)


def test_G2_the_collapse_fix_is_proven_not_asserted():
    """THE FINDING / definition of done: under a benign rescaling of the reference state, the
    RELATIVE cutoff (the actual fix) keeps an IDENTICAL number of basis vectors at every scale --
    proving scale-invariance, required for a well-posed generalized eigenproblem. An independently
    REIMPLEMENTED absolute cutoff at the identical numeric threshold -- the module's former, buggy
    behavior, documented in its docstring but never regression-tested -- collapses to ZERO kept
    vectors at small scale, reproducing "dropped almost every vector and nulled the eigenproblem"."""
    h1, eri, e_core, nelec, norb, _casci = _reference("H4")
    cutoff = 1e-8
    kept_relative, kept_absolute = [], []
    for scale in (1.0, 1e-3, 1e-6, 1e-9):
        S = _independent_S(h1, eri, e_core, nelec, norb, dt=0.5, krylov_dim=12, phi0_scale=scale)
        s_eig = np.linalg.eigvalsh(S).real
        kept_relative.append(int(np.sum(s_eig > cutoff * s_eig.max())))
        kept_absolute.append(int(np.sum(s_eig > cutoff)))  # the OLD, buggy logic

    assert len(set(kept_relative)) == 1, kept_relative
    assert kept_relative[0] > 0, kept_relative
    assert kept_absolute[-1] == 0, kept_absolute       # scale=1e-9: fully collapsed
    assert kept_absolute[-2] == 0, kept_absolute       # scale=1e-6: fully collapsed
    assert kept_absolute[0] > kept_absolute[-1], kept_absolute  # scale=1: not collapsed


def test_G3_cross_check_against_the_independent_qubit_mapped_implementation():
    """The cross-check the module's own docstring says is its purpose, run for the first time:
    this module's converged energy (FCI-direct) agrees with
    hybrid_quantum_solver.QuantumKrylovSolver's (qubit-mapped statevector, an entirely independent
    code path) to sub-mHa precision on the same active space."""
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

    h1, eri, e_core, nelec, norb, _casci = _reference("H4")
    e_fci_direct, _info = krylov_ground_state(h1, eri, e_core, nelec, norb, krylov_dim=16, dt=0.5)

    mh = build_molecular_hamiltonian(atom=SYSTEMS["H4"][0], basis="sto-3g",
                                     active_electrons=4, active_orbitals=4)
    step = QuantumKrylovSolver(mh).solve(16)

    assert abs(e_fci_direct - step.energy) * 1e3 < 1.0, (e_fci_direct, step.energy)


def test_G4_condition_number_claim_is_measured_not_asserted():
    """The docstring's "condition numbers of 1e6+ are normal and harmless" claim, checked: sweeping
    dt reaches condition numbers >= 1e6 while energy error stays bounded at every tested dt."""
    h1, eri, e_core, nelec, norb, casci = _reference("H4")
    saw_large_condition = False
    for dt in (0.1, 0.3, 0.5, 1.0, 2.0):
        e, info = krylov_ground_state(h1, eri, e_core, nelec, norb, krylov_dim=12, dt=dt)
        err_mha = abs(e - casci) * 1e3
        assert err_mha < 3.0, (dt, err_mha, info)
        if info["overlap_condition"] >= 1e6:
            saw_large_condition = True
    assert saw_large_condition
