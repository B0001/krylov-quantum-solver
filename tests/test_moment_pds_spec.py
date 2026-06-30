"""
Acceptance gates G1-G4 for specs/SPEC_moment_pds.md (Hamiltonian-moment energies: PDS and CMX).

Test-first: ``moment_expansion`` does not exist yet, so this file is RED until the spec is
implemented. PDS(K) estimates the ground-state energy from the Hamiltonian moments
<phi|H^n|phi> of the HF reference and is a *variational upper bound* converging to FCI; CMX(2) is
the same moment data resummed non-variationally (it can dip below FCI). Reference is FCI (dense
diagonalization of the same qubit Hamiltonian).

PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from moment_expansion import cmx2_energy, hamiltonian_moments, pds_energy

CHEM_ACC = 1.6e-3


def _cases():
    return {
        "H2": build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74"),
        "H4": build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
        "LiH": build_molecular_hamiltonian(atom="Li 0 0 0; H 0 0 1.6"),
    }


def test_G1_pds1_is_hf_expectation():
    """PDS(1) equals <H> = the Rayleigh quotient of the HF reference (RHF energy)."""
    for mh in _cases().values():
        mu, off = hamiltonian_moments(mh, 2)
        assert abs(pds_energy(mu, 1, off) - (mu[1] + off)) < 1e-9


def test_G2_pds_variational_and_converges():
    """DEFINITION OF DONE: PDS(K) >= FCI at every K, and PDS(4) reaches chemical accuracy."""
    for mh in _cases().values():
        fci = mh.ground_state_energy()
        mu, off = hamiltonian_moments(mh, 7)
        for K in (1, 2, 3, 4):
            e = pds_energy(mu, K, off)
            assert e >= fci - 1e-9, (mh.num_qubits, K, e, fci)            # variational
        assert abs(pds_energy(mu, 4, off) - fci) < CHEM_ACC, (mh.num_qubits, pds_energy(mu, 4, off), fci)


def test_G3_pds_monotone_tightening():
    """Higher PDS order gives a tighter (non-looser) upper bound."""
    for mh in _cases().values():
        mu, off = hamiltonian_moments(mh, 7)
        es = [pds_energy(mu, K, off) for K in (1, 2, 3, 4)]
        assert all(es[k + 1] <= es[k] + 1e-9 for k in range(3)), es


def test_G4_cmx_is_not_variational():
    """THE BOUNDARY: on H2 the connected-moment expansion dips BELOW FCI, unlike PDS."""
    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    fci = mh.ground_state_energy()
    mu, off = hamiltonian_moments(mh, 4)
    assert cmx2_energy(mu, off) < fci - 1e-4, (cmx2_energy(mu, off), fci)   # non-variational
    assert pds_energy(mu, 2, off) >= fci - 1e-9, (pds_energy(mu, 2, off), fci)
