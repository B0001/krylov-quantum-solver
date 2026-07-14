"""
Acceptance gates G1-G4 for specs/SPEC_taper_spectrum.md (Z2 qubit tapering preserves the FULL
sector spectrum, independently verified -- not just the ground energy `taper_qubits.py`'s own
`__main__` checks).

Deliberately no new library code: the independent verification lives entirely here, reusing
`taper_qubits.py`'s existing `find_symmetries`/`_gf2_independent`/`taper_hamiltonian` unmodified,
so it is a genuine external check rather than a round-trip through the same machinery.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from adapt_vqe import hf_state
from qubitization_blueprint import PAULI, _kron, build_qubit_hamiltonian, pauli_decompose
from taper_qubits import _gf2_independent, find_symmetries, taper_hamiltonian

SYSTEMS = ("H2", "LiH", "H3_radical")


def _reference(label):
    if label == "H2":
        mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g")
        mf = scf.RHF(mol)
        norb, ne, nelec = 2, 2, (1, 1)
    elif label == "LiH":
        mol = gto.M(atom="Li 0 0 0; H 0 0 1.6", basis="sto-3g")
        mf = scf.RHF(mol)
        norb, ne, nelec = 2, 2, (1, 1)
    else:  # open-shell radical -- not exercised by taper_qubits.py's own __main__
        mol = gto.M(atom="H 0 0 0; H 0 0 1; H 0 0 2", basis="sto-3g", spin=1)
        mf = scf.ROHF(mol)
        norb, ne, nelec = 3, 3, (2, 1)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, norb, ne)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return h1, eri, float(e_core), nelec, norb, float(cas.e_tot)


def _z_indep(h1, eri, e_core, nelec, norb):
    """The qubit Hamiltonian, HF state, and the independent Z-type symmetry generators --
    everything `taper_hamiltonian` computes internally, recomputed here for an external check."""
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    hf = hf_state(nelec[0], nelec[1], n)
    terms = pauli_decompose(H, n)
    syms = find_symmetries(terms, n)
    z_syms = [z for (x, z) in syms if not x.any() and z.any()]
    z_indep = _gf2_independent(np.array(z_syms)) if z_syms else []
    return H, n, hf, z_indep


def _independent_sector_spectrum(H, n, ref_state, z_indep):
    """Sector-projected spectrum built directly from computational-basis parity of each Z-type
    generator on ``ref_state`` -- no Clifford rotation, so this is NOT `taper_hamiltonian`'s own
    code path."""
    P = np.eye(2 ** n, dtype=complex)
    for g in z_indep:
        Zg = _kron([PAULI["Z"] if g[q] else PAULI["I"] for q in range(n)])
        s = float(np.real(np.vdot(ref_state, Zg @ ref_state)))
        assert abs(abs(s) - 1.0) < 1e-9, "reference state is not a Z_g eigenstate"
        P = P @ (np.eye(2 ** n) + s * Zg) / 2.0
    w, V = np.linalg.eigh((P + P.conj().T) / 2)
    Vk = V[:, w > 0.5]
    Hp = Vk.conj().T @ H @ Vk
    return np.sort(np.linalg.eigvalsh(Hp).real)


def _flipped_basis_state(hf, n):
    """A different computational basis string (flip the top qubit) -- a wrong reference for the
    independent projector, used only to prove G1 can fail (G3)."""
    idx = int(np.argmax(np.abs(hf)))
    flipped = np.zeros_like(hf)
    flipped[idx ^ (1 << (n - 1))] = 1.0
    return flipped


@pytest.mark.parametrize("label", SYSTEMS)
def test_G1_full_spectrum_matches_independent_projection(label):
    """The tapered Hamiltonian's FULL spectrum (not just the ground state) matches an
    independently-constructed sector projection to high precision."""
    h1, eri, e_core, nelec, norb, _casci = _reference(label)
    H, n, hf, z_indep = _z_indep(h1, eri, e_core, nelec, norb)
    r = taper_hamiltonian(h1, eri, e_core, nelec, norb)
    tapered_spec = np.sort(np.linalg.eigvalsh(r["H_tapered"]).real)
    indep_spec = _independent_sector_spectrum(H, n, hf, z_indep)
    assert tapered_spec.shape == indep_spec.shape, (label, tapered_spec.shape, indep_spec.shape)
    assert np.allclose(tapered_spec, indep_spec, atol=1e-7), (label, tapered_spec, indep_spec)


@pytest.mark.parametrize("label", SYSTEMS)
def test_G2_reduction_matches_independent_symmetry_count(label):
    """Qubits removed == the independent Z-symmetry count exactly -- the reduction mechanism is
    internally consistent, not silently over/under-reducing."""
    h1, eri, e_core, nelec, norb, _casci = _reference(label)
    _H, n, _hf, z_indep = _z_indep(h1, eri, e_core, nelec, norb)
    r = taper_hamiltonian(h1, eri, e_core, nelec, norb)
    assert n - r["n_qubits_tapered"] == len(z_indep), (label, n, r["n_qubits_tapered"], z_indep)


@pytest.mark.parametrize("label", ("H2", "LiH"))
def test_G3_wrong_sector_gives_a_different_spectrum(label):
    """THE FINDING / definition of done: the same independent construction with a WRONG reference
    state gives a spectrum that does NOT match H_tapered's -- proof G1 can fail, not a vacuous
    check that trivially matches any sector."""
    h1, eri, e_core, nelec, norb, _casci = _reference(label)
    H, n, hf, z_indep = _z_indep(h1, eri, e_core, nelec, norb)
    r = taper_hamiltonian(h1, eri, e_core, nelec, norb)
    tapered_spec = np.sort(np.linalg.eigvalsh(r["H_tapered"]).real)

    wrong_ref = _flipped_basis_state(hf, n)
    wrong_spec = _independent_sector_spectrum(H, n, wrong_ref, z_indep)
    same_shape = wrong_spec.shape == tapered_spec.shape
    assert not (same_shape and np.allclose(wrong_spec, tapered_spec, atol=1e-6)), label


def test_G4_open_shell_reference_exercised():
    """The H3 radical case (spin=1, na != nb) is genuinely open-shell -- taper_qubits.py's own
    __main__ only ever runs closed-shell RHF systems; this is the first check on that scope."""
    _h1, _eri, _e_core, nelec, _norb, _casci = _reference("H3_radical")
    assert nelec[0] != nelec[1], "H3 radical fixture is not actually open-shell"
