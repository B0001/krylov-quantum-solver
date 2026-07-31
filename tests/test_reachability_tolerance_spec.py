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
    # Pin the measured witness so a silent drift in either value is loud.
    assert ov10 == pytest.approx(2.25e-5, rel=0.05), ov10
    assert ov08 == pytest.approx(0.667, rel=0.05), ov08


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
