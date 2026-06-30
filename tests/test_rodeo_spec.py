"""
Acceptance gates G1-G4 for specs/SPEC_rodeo.md (the rodeo algorithm spectral filter).

Test-first: ``rodeo`` does not exist yet, so this file is RED until the spec is implemented. The
rodeo algorithm filters the spectrum with K cycles of random-time evolution; the expected survival
probability P_bar(E) peaks at the eigenvalues with height = reference overlap, the ground energy is
the dominant low-energy peak, and the filter sharpens / suppresses off-resonance as K grows.
Reference: FCI (dense diagonalization of the same qubit Hamiltonian).

PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from rodeo import reference_spectrum, rodeo_ground_energy, rodeo_survival


def _h2():
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")


def _h4():
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")


def test_G1_ground_recovery():
    """DEFINITION OF DONE: the dominant rodeo peak recovers the FCI ground-state energy."""
    for mh in (_h2(), _h4()):
        e = rodeo_ground_energy(mh, sigma=2.0, n_cycles=12)
        assert abs(e - mh.ground_state_energy()) < 3e-3, (mh.num_qubits, e, mh.ground_state_energy())


def test_G2_filter_sharpens_with_cycles():
    """More cycles narrow the peak: H4 ground-energy error at K=12 beats K=3."""
    mh = _h4()
    fci = mh.ground_state_energy()
    err_lo = abs(rodeo_ground_energy(mh, n_cycles=3) - fci)
    err_hi = abs(rodeo_ground_energy(mh, n_cycles=12) - fci)
    assert err_hi < err_lo, (err_hi, err_lo)


def test_G3_off_resonance_suppression():
    """P_bar at a non-eigenvalue decreases monotonically with K (the (<1)^K filter)."""
    mh = _h4()
    w, ov, _ = reference_spectrum(mh)
    e0 = w[ov > 1e-8].min()
    e_off = e0 + 0.7                                  # between eigenvalues
    bg = [rodeo_survival(w, ov, e_off, sigma=2.0, n_cycles=K) for K in (3, 6, 12)]
    assert bg[2] < bg[1] < bg[0], bg


def test_G4_peak_height_is_reference_overlap():
    """THE FINDING: the ground-peak height equals the reference overlap |<HF|E0>|^2."""
    for mh in (_h2(), _h4()):
        w, ov, _ = reference_spectrum(mh)
        i0 = w.argmin()
        e0, overlap0 = w[i0], ov[i0]
        peak = rodeo_survival(w, ov, e0, sigma=2.0, n_cycles=12)
        assert abs(peak - overlap0) < 0.02, (mh.num_qubits, peak, overlap0)
