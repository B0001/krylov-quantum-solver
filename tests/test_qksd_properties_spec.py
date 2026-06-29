"""
Acceptance gates G1-G4 for specs/SPEC_qksd_properties.md (molecular properties from QKSD).

Test-first: ``QuantumKrylovSolver.eigenstates`` and ``hybrid_quantum_solver.qksd_properties`` do
not exist yet, so this file is RED until the spec is implemented. The same Krylov subspace that
gives the excited-state energies also gives their eigenvectors, hence dipoles, transition dipoles,
and oscillator strengths -- validated against the same matrix elements between the dense-diagonalized
exact eigenstates (the FCI reference for properties).

HeH+ (polar, 4 qubits) carries bright transitions and a nonzero permanent dipole; H2 (centrosymmetric,
4 qubits) is the convention-free anchor: zero permanent dipole, dark g->g transition.

Uses only pyscf/qiskit (no block2), so it runs in the non-DMRG process group -- but `make gates`
runs every test_*_spec.py in its own process anyway.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import (
    build_dipole_operators,
    build_molecular_hamiltonian,
)
from hybrid_quantum_solver.qksd_properties import (
    oscillator_strengths,
    property_matrix,
    transition_dipoles,
)
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

DIP_TOL = 1e-3  # a.u.
DEPTH = 12      # reachable subspace (<=3 dim here) is saturated well before this


def _heh_plus():
    return dict(atom="He 0 0 0; H 0 0 0.772", charge=1)


def _h2():
    return dict(atom="H 0 0 0; H 0 0 0.74")


def _reachable_states(mh, overlap_tol=1e-8):
    """Dense-diagonalization reference: exact eigen-energies + eigenstates reachable from |HF>."""
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    overlaps = np.abs(V.conj().T @ hf) ** 2
    idx = np.where(overlaps > overlap_tol)[0]
    idx = idx[np.argsort(w[idx].real)]
    energies = w[idx].real + mh.energy_offset
    states = V[:, idx].T                              # rows are |v_m>
    return energies, states


def _dip_matrices(spec):
    return [op.to_matrix(sparse=True) for op in build_dipole_operators(**spec)]


def test_G1_permanent_dipole_vs_fci():
    """Ground-state permanent dipole: nonzero & FCI-accurate for HeH+, exactly zero for H2."""
    heh = _heh_plus()
    mh = build_molecular_hamiltonian(**heh)
    dmats = _dip_matrices(heh)
    _, states = QuantumKrylovSolver(mh).eigenstates(DEPTH, n_states=1)
    _, ref_states = _reachable_states(mh)
    qksd = transition_dipoles(states, dmats)[:, 0, 0].real
    exact = transition_dipoles(ref_states[:1], dmats)[:, 0, 0].real
    assert np.linalg.norm(qksd - exact) < DIP_TOL, (qksd, exact)
    assert abs(qksd[2]) > 0.5, qksd                  # polar molecule: nonzero z-dipole

    mh2 = build_molecular_hamiltonian(**_h2())
    _, s2 = QuantumKrylovSolver(mh2).eigenstates(DEPTH, n_states=1)
    perm_h2 = transition_dipoles(s2, _dip_matrices(_h2()))[:, 0, 0]
    assert np.max(np.abs(perm_h2)) < 1e-6, perm_h2   # centrosymmetric: zero by symmetry


def test_G2_transition_dipoles_vs_fci():
    """Ground->excited transition-dipole magnitudes match FCI; first HeH+ transition is bright."""
    heh = _heh_plus()
    mh = build_molecular_hamiltonian(**heh)
    dmats = _dip_matrices(heh)
    ref_e, ref_states = _reachable_states(mh)
    k = len(ref_e)
    _, states = QuantumKrylovSolver(mh).eigenstates(DEPTH, n_states=k)
    qksd = np.linalg.norm(transition_dipoles(states, dmats)[:, 0, :], axis=0)     # |mu_0n|
    exact = np.linalg.norm(transition_dipoles(ref_states, dmats)[:, 0, :], axis=0)
    assert np.allclose(qksd, exact, atol=DIP_TOL), (qksd, exact)
    assert qksd[1] > 0.5, qksd                        # the bright transition (measured ~0.85)


def test_G3_oscillator_strengths_and_dark_recovery():
    """DEFINITION OF DONE: oscillator strengths match FCI; H2's g->g transition stays dark."""
    heh = _heh_plus()
    mh = build_molecular_hamiltonian(**heh)
    dmats = _dip_matrices(heh)
    ref_e, ref_states = _reachable_states(mh)
    k = len(ref_e)
    energies, states = QuantumKrylovSolver(mh).eigenstates(DEPTH, n_states=k)
    f_qksd = oscillator_strengths(energies, states, dmats)
    f_exact = oscillator_strengths(ref_e, ref_states, dmats)
    assert np.allclose(f_qksd, f_exact, atol=1e-3), (f_qksd, f_exact)
    assert f_qksd[1] > 1e-3, f_qksd                   # bright state carries oscillator strength

    # H2: the reachable excitation is dipole-forbidden -> dark
    mh2 = build_molecular_hamiltonian(**_h2())
    e2, s2 = QuantumKrylovSolver(mh2).eigenstates(DEPTH, n_states=2)
    f_h2 = oscillator_strengths(e2, s2, _dip_matrices(_h2()))
    assert f_h2[1] < 1e-6, f_h2


def test_G4_hermiticity_and_normalization():
    """Can't-be-faked invariants: property matrices are Hermitian, eigenstates are normalized."""
    for spec in (_heh_plus(), _h2()):
        mh = build_molecular_hamiltonian(**spec)
        dmats = _dip_matrices(spec)
        _, states = QuantumKrylovSolver(mh).eigenstates(DEPTH, n_states=2)
        norms = np.einsum("mi,mi->m", states.conj(), states).real
        assert np.allclose(norms, 1.0, atol=1e-9), norms
        for op in dmats:
            mat = property_matrix(states, op)
            assert np.linalg.norm(mat - mat.conj().T) < 1e-9, np.linalg.norm(mat - mat.conj().T)
