"""Gates for specs/SPEC_reachability_tolerance.md.

The certified arc thresholds the HF-reachable sector at |<HF|psi_k>|^2 > tol, with `tol` written as
a magic number independently in each file -- 1e-10 in certified_gaps / hf_overlap_certificate /
certified_dipole / certified_noise, 1e-8 in hf_overlap_subspace. At square H4 a = 1.1 A the two
pick DIFFERENT ground states, so the two specs whose head-to-head is SPEC_hf_overlap_subspace's
headline are certifying different targets.

Found while verifying the krylov_refine chained-overlap bound: the bound appeared to violate its
reference, and the cause was the reference, not the bound.
"""
import numpy as np
import pytest

from hf_overlap_certificate import REACHABLE_TOL_CERTIFIED, exact_reachable_overlap
from hf_overlap_subspace import _REACHABLE_TOL, exact_hf_subspace_overlap
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

# Verbatim from hf_overlap_subspace.py's own __main__ sweep -- atom ORDER matters (R2).
def _square_h4(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"

WITNESS_A = 1.1                    # the geometry where the thresholds diverge
AGREE_A = (1.0, 1.2, 1.3, 1.4)     # where they do not -- the boundary


def _amp2_and_spectrum(a):
    mh = build_molecular_hamiltonian(atom=_square_h4(a))
    H = mh.qubit_hamiltonian.to_matrix()
    w, V = np.linalg.eigh(H)
    u = np.asarray(mh.hf_state().data, dtype=complex)
    return w, np.abs(V.conj().T @ u) ** 2, mh


def _lowest_reachable_index(w, amp2, tol):
    reach = np.where(amp2 > tol)[0]        # eigh returns ascending eigenvalues
    return int(reach[0])


# --- G1: the witness (DEFINITION OF DONE) ---------------------------------------------------------

def test_G1_two_thresholds_select_different_ground_states():
    w, amp2, _ = _amp2_and_spectrum(WITNESS_A)
    i10 = _lowest_reachable_index(w, amp2, 1e-10)
    i08 = _lowest_reachable_index(w, amp2, 1e-8)
    assert i10 != i08, (i10, i08)
    ov10, ov08 = np.sqrt(amp2[i10]), np.sqrt(amp2[i08])
    assert ov08 / ov10 > 1000.0, (ov10, ov08)
    # ov08 is physics -- the Ag ground state's HF weight -- so it is safe to pin.
    assert ov08 == pytest.approx(0.667, rel=0.05), ov08
    # ov10 is NOT pinned. The original version of this gate asserted ov10 ~ 2.25e-5, which pinned
    # an SCF CONVERGENCE RESIDUE (see G6): it moves 19 orders of magnitude with PySCFDriver's
    # conv_tol and would silently break on any pyscf/qiskit-nature bump. Assert only that it is
    # far below the physical scale, which is the durable statement.
    assert ov10 < 1e-3, ov10


# --- G2: the two SHIPPED modules disagree, measured through their public references ---------------

def test_G2_the_two_specs_d1_references_disagree():
    """`exact_reachable_overlap` (SPEC_hf_overlap_certificate, 1e-10) vs
    `exact_hf_subspace_overlap(..., 1)` (SPEC_hf_overlap_subspace, 1e-8) -- the d=1 references of
    the two specs whose comparison is the latter's headline finding.
    """
    mh = build_molecular_hamiltonian(atom=_square_h4(WITNESS_A))
    a = exact_reachable_overlap(mh)                 # tol 1e-10
    b = exact_hf_subspace_overlap(mh, 1)            # tol _REACHABLE_TOL = 1e-8
    assert b / max(a, 1e-300) > 1000.0, (a, b)


# --- G3: the offending level really is between the thresholds -------------------------------------

def test_G3_a_level_sits_between_the_two_thresholds():
    _, amp2, _ = _amp2_and_spectrum(WITNESS_A)
    between = np.where((amp2 > 1e-10) & (amp2 < 1e-8))[0]
    assert between.size >= 1, amp2[amp2 > 1e-12]


# --- G4: the boundary -- near-threshold, not pervasive --------------------------------------------

@pytest.mark.parametrize("a", AGREE_A)
def test_G4_thresholds_agree_away_from_the_witness(a):
    """Bounds the claim to existence, not prevalence. Killed if they disagree everywhere, which
    would be a far larger claim than the evidence supports."""
    w, amp2, _ = _amp2_and_spectrum(a)
    assert _lowest_reachable_index(w, amp2, 1e-10) == _lowest_reachable_index(w, amp2, 1e-8)


# --- G5: the constant is pinned to the literal it documents ---------------------------------------

def test_G5_certified_tolerance_constant_is_pinned():
    assert REACHABLE_TOL_CERTIFIED == 1e-10
    assert _REACHABLE_TOL == 1e-8
    assert REACHABLE_TOL_CERTIFIED != _REACHABLE_TOL     # the finding, asserted not assumed


def test_G5_constant_governs_the_certified_reference():
    """`exact_reachable_overlap` must actually USE the constant -- otherwise G5 pins a decoration."""
    mh = build_molecular_hamiltonian(atom=_square_h4(WITNESS_A))
    w, amp2, _ = _amp2_and_spectrum(WITNESS_A)
    expected = float(np.sqrt(amp2[_lowest_reachable_index(w, amp2, REACHABLE_TOL_CERTIFIED)]))
    assert exact_reachable_overlap(mh) == pytest.approx(expected, rel=1e-12)


# =================================================================================================
# FALSIFICATION (2026-07-31). The gates above stand as OBSERVATION -- the two thresholds really do
# select different states. But this spec's original physics framing ("neither value is wrong, they
# answer different questions") is FALSE, and the gates below are why. See §2b of the spec.
# =================================================================================================

def _p0_at_conv_tol(conv_tol):
    """Lowest-level HF population, as a function of how tightly the SCF was converged."""
    from qiskit.quantum_info import Statevector
    from qiskit_nature.second_q.circuit.library import HartreeFock
    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import JordanWignerMapper

    prob = PySCFDriver(atom=_square_h4(WITNESS_A), basis="sto-3g", conv_tol=conv_tol).run()
    mapper = JordanWignerMapper()
    H = mapper.map(prob.hamiltonian.second_q_op()).to_matrix()
    w, V = np.linalg.eigh(H)
    hf = np.asarray(Statevector(HartreeFock(prob.num_spatial_orbitals, prob.num_particles,
                                            mapper)).data, dtype=complex)
    return float(w[0]), float((np.abs(V.conj().T @ hf) ** 2)[0])


def test_G6_the_disputed_amplitude_is_an_SCF_convergence_residue():
    """THE FALSIFIER. A physical overlap does not depend on how tightly the SCF was converged.
    This one moves ~19 orders of magnitude while the eigenvalue is unchanged to 10 digits.
    """
    e_loose, p_loose = _p0_at_conv_tol(1e-6)
    e_tight, p_tight = _p0_at_conv_tol(1e-13)
    assert e_loose == pytest.approx(e_tight, abs=1e-9), (e_loose, e_tight)   # same state
    assert p_loose > 1e-9, p_loose                                            # admitted at 1e-10
    assert p_tight < 1e-20, p_tight                                           # collapses to zero
    assert p_loose / p_tight > 1e10, (p_loose, p_tight)


def test_G7_the_state_is_symmetry_forbidden_to_the_HF_determinant():
    """THE MECHANISM. Square H4 is D2h; the HF determinant is Ag, the disputed level is B1g, and
    its HF-determinant coefficient is EXACTLY zero. So the true overlap is 0, not 2.25e-5.
    """
    from pyscf import fci, gto, scf

    mol = gto.M(atom=_square_h4(WITNESS_A), basis="sto-3g", symmetry=True, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.kernel()
    coeffs = {}
    for wfnsym in ("Ag", "B1g"):
        solver = fci.FCI(mf)
        solver.wfnsym = wfnsym
        solver.nroots = 1
        e, civec = solver.kernel()
        c0 = float(np.atleast_2d(np.asarray(civec))[0, 0])
        coeffs[wfnsym] = (float(np.atleast_1d(e)[0]), c0 ** 2)
    # The B1g state is LOWER in energy but carries exactly zero HF weight.
    assert coeffs["B1g"][0] < coeffs["Ag"][0], coeffs
    assert coeffs["B1g"][1] < 1e-20, coeffs["B1g"]
    # ...while the Ag ground state -- the one tol=1e-8 selects -- carries the real weight.
    assert coeffs["Ag"][1] == pytest.approx(0.4451, rel=0.05), coeffs["Ag"]


def test_G8_neither_constant_is_safe():
    """The killer for 'just pick the right constant'. There is a geometry in the same family where
    the artifact exceeds even the LOOSER 1e-8 threshold, so no fixed tolerance separates physics
    from SCF residue. The fix has to be symmetry/sector-aware, not a better number.
    """
    offenders = []
    for a in (1.185, 1.190, 1.195):
        _, amp2, _ = _amp2_and_spectrum(a)
        if amp2[0] > 1e-8:
            offenders.append((a, float(amp2[0])))
    assert offenders, "no geometry breached 1e-8 -- the 'no safe constant' claim would die here"
