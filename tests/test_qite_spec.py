"""
Acceptance gates G1-G4 for specs/SPEC_qite.md (quantum imaginary-time evolution).

Test-first: ``qite`` does not exist yet, so this file is RED until the spec is implemented. QITE
replaces the non-unitary imaginary-time step e^{-dtau H} with a unitary e^{-i dtau A}, A found from
the McLachlan system S a = b. We validate: exact imaginary-time evolution is variational and reaches
FCI; the full-domain QITE update reproduces it (its equations are correct); the step error vanishes
with dtau; and a truncated (low-weight) operator domain stalls at Hartree-Fock -- QITE's accuracy is
set by the domain. Reference: exact ITE and FCI.

PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from qite import exact_imaginary_time, pauli_operators, qite_evolve


def _h2():
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")


def test_G1_exact_ite_variational_and_converges():
    """Exact imaginary-time energy is monotone non-increasing and variational (>= FCI); it converges
    to FCI at a beta set by the gap (H2 fast, LiH's small gap needs larger beta)."""
    lih = build_molecular_hamiltonian(atom="Li 0 0 0; H 0 0 1.6")
    for mh in (_h2(), lih):
        fci = mh.ground_state_energy()
        es = exact_imaginary_time(mh, [0.5 * k for k in range(1, 9)])   # beta 0.5..4
        assert all(es[i + 1] <= es[i] + 1e-12 for i in range(len(es) - 1)), es   # monotone
        assert all(e >= fci - 1e-9 for e in es), es                              # variational
    assert abs(exact_imaginary_time(_h2(), [4.0])[0] - _h2().ground_state_energy()) < 1e-4
    assert abs(exact_imaginary_time(lih, [15.0])[0] - lih.ground_state_energy()) < 1e-3


def test_G2_full_domain_qite_reaches_fci():
    """DEFINITION OF DONE: full-domain QITE reproduces exact ITE and lands at FCI (update is correct)."""
    mh = _h2()
    fci = mh.ground_state_energy()
    full = pauli_operators(mh.num_qubits)
    e_qite = qite_evolve(mh, 0.1, 40, full)
    e_exact = exact_imaginary_time(mh, [0.1 * k for k in range(1, 41)])
    assert abs(e_qite[-1] - fci) < 1e-4, (e_qite[-1], fci)
    assert abs(e_qite[-1] - e_exact[-1]) < 1e-4, (e_qite[-1], e_exact[-1])


def test_G3_step_error_vanishes_with_dtau():
    """Smaller dtau tracks exact ITE better (O(dtau^2) single-step error) at matched beta=2."""
    mh = _h2()
    full = pauli_operators(mh.num_qubits)
    e_exact = exact_imaginary_time(mh, [2.0])[0]
    err_coarse = abs(qite_evolve(mh, 0.2, 10, full)[-1] - e_exact)   # beta = 2.0
    err_fine = abs(qite_evolve(mh, 0.05, 40, full)[-1] - e_exact)    # beta = 2.0
    assert err_fine < err_coarse, (err_fine, err_coarse)


def test_G4_truncated_domain_stalls_at_hf():
    """THE FINDING: a weight-<=2 domain stalls at Hartree-Fock; only the full domain reaches FCI."""
    mh = _h2()
    fci = mh.ground_state_energy()
    hf_err = mh.hf_energy - fci                                       # ~ +20.5 mHa
    e_trunc = qite_evolve(mh, 0.05, 80, pauli_operators(mh.num_qubits, max_weight=2))
    e_full = qite_evolve(mh, 0.05, 80, pauli_operators(mh.num_qubits))
    assert abs((e_trunc[-1] - fci) - hf_err) < 1e-3, (e_trunc[-1] - fci, hf_err)   # stuck at HF
    assert abs(e_full[-1] - fci) < 1e-4, (e_full[-1], fci)                          # full reaches FCI
