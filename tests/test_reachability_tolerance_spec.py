"""Gates for specs/SPEC_reachability_tolerance.md.

The certified arc thresholds the HF-reachable sector at |<HF|psi_k>|^2 > tol, with `tol` written as
a magic number independently in each file -- 1e-10 in certified_gaps / hf_overlap_certificate /
certified_dipole / certified_noise, 1e-8 in hf_overlap_subspace. At square H4 a = 1.1 A the two
pick DIFFERENT ground states, so the two specs whose head-to-head is SPEC_hf_overlap_subspace's
headline are certifying different targets.

Found while verifying the krylov_refine chained-overlap bound: the bound appeared to violate its
reference, and the cause was the reference, not the bound.

chem-8ls (2026-11-05): that residue level existed only while the free SCF stopped on square H4's
D2h saddle point, and which sides it stopped on depended on the platform's rounding. The builder
now follows the instability to the stable broken-symmetry minimum
(specs/SPEC_scf_stability_rule.md). The lowest level is then physical at every side, so G1, G2, G6
and G8 assert that the two cuts agree, and G3 is retired (SPEC section 11).
"""
import numpy as np
import pytest

from hf_overlap_certificate import REACHABLE_TOL_CERTIFIED, exact_reachable_overlap
from hf_overlap_subspace import _REACHABLE_TOL, exact_hf_subspace_overlap
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from reachability import scf_symmetry_status

# Verbatim from hf_overlap_subspace.py's own __main__ sweep -- atom ORDER matters (R2).
def _square_h4(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"

WITNESS_A = 1.1                    # the geometry where the thresholds diverged until chem-8ls
AGREE_A = (1.0, 1.2, 1.3, 1.4)     # where they did not -- the boundary

# G1-G3 used to pin a platform-specific SCF STOPPING POINT and were Linux-reference-only (macOS put
# no residue on the a=1.1 lowest level, section 10). The stopping point was the D2h saddle; the
# stability rule removed it on every platform, so G1-G2 now run everywhere (section 11).


def _amp2_and_spectrum(a):
    mh = build_molecular_hamiltonian(atom=_square_h4(a))
    H = mh.qubit_hamiltonian.to_matrix()
    w, V = np.linalg.eigh(H)
    u = np.asarray(mh.hf_state().data, dtype=complex)
    return w, np.abs(V.conj().T @ u) ** 2, mh


def _lowest_reachable_index(w, amp2, tol):
    reach = np.where(amp2 > tol)[0]        # eigh returns ascending eigenvalues
    return int(reach[0])


# --- G1: the former witness -- both cuts now select the same, physical ground state (chem-8ls) ----

def test_G1_two_thresholds_select_the_same_ground_state():
    """Was the witness (DEFINITION OF DONE): at a = 1.1 the cuts selected different ground states,
    HF overlap 2.25e-5 vs 0.667 (Linux freeze, section 2). The 1e-10 cut's level carried SCF residue
    from a free SCF stopped on the D2h saddle. With the stability rule the reference is the
    broken-symmetry minimum on every platform, and both cuts select its lowest level, which HF
    genuinely overlaps (0.685; section 11)."""
    w, amp2, _ = _amp2_and_spectrum(WITNESS_A)
    assert _lowest_reachable_index(w, amp2, 1e-10) == _lowest_reachable_index(w, amp2, 1e-8) == 0
    assert np.sqrt(amp2[0]) == pytest.approx(0.685, abs=0.005), amp2[0]


# --- G2: the two SHIPPED modules now agree, measured through their public references --------------

def test_G2_the_two_specs_d1_references_agree():
    """`exact_reachable_overlap` (SPEC_hf_overlap_certificate, 1e-10) vs
    `exact_hf_subspace_overlap(..., 1)` (SPEC_hf_overlap_subspace, 1e-8) -- the d=1 references of
    the two specs whose comparison is the latter's headline finding. They differed by > 1000x here
    until chem-8ls; the d=1 vs d=2 head-to-head now certifies one target at this geometry.
    """
    mh = build_molecular_hamiltonian(atom=_square_h4(WITNESS_A))
    a = exact_reachable_overlap(mh)                 # tol 1e-10
    b = exact_hf_subspace_overlap(mh, 1)            # tol _REACHABLE_TOL = 1e-8
    assert a == pytest.approx(b, rel=1e-8), (a, b)


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
# FALSIFICATION (2026-07-31). The gates above stood as OBSERVATION -- the two thresholds really did
# select different states. But this spec's original physics framing ("neither value is wrong, they
# answer different questions") is FALSE, and the gates below are why. See §2b of the spec. Since
# chem-8ls (section 11) G6 and G8 gate the remedy at the source: the builder's reference no longer
# carries the residue.
# =================================================================================================

def _p0_at_conv_tol(conv_tol):
    """(lowest eigenvalue, its HF population) of the builder's reference at this conv_tol."""
    mh = build_molecular_hamiltonian(atom=_square_h4(WITNESS_A), conv_tol=conv_tol)
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    return float(w[0]), float((np.abs(V.conj().T @ hf) ** 2)[0])


def test_G6_the_lowest_population_does_not_depend_on_conv_tol():
    """Was THE FALSIFIER: through a raw PySCFDriver this population moved ~19 orders of magnitude
    with conv_tol while the eigenvalue stayed put (1.49e-9 at 1e-6, 1.53e-28 at 1e-13, Linux
    freeze), so it was an SCF residue, not an overlap (section 2b(i)). That driver had neither the
    thread pin nor the stability rule, so which side of the D2h saddle it stopped on varied by
    platform and by run (chem-8ls). The builder's reference now has a physical overlap here, and a
    physical overlap does not depend on how tightly the SCF was converged (section 11)."""
    e_loose, p_loose = _p0_at_conv_tol(1e-6)
    e_tight, p_tight = _p0_at_conv_tol(1e-13)
    assert e_loose == pytest.approx(e_tight, abs=1e-9), (e_loose, e_tight)   # same state
    assert p_loose == pytest.approx(p_tight, abs=1e-4), (p_loose, p_tight)    # does not move
    assert p_tight > 0.4, p_tight                                             # and is physical


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
        if not solver.converged:
            raise RuntimeError(f"FCI Davidson did not converge for wfnsym={wfnsym!r} "
                                f"(chem-60i: an unchecked Davidson call can silently report a wrong energy)")
        c0 = float(np.atleast_2d(np.asarray(civec))[0, 0])
        coeffs[wfnsym] = (float(np.atleast_1d(e)[0]), c0 ** 2)
    # The B1g state is LOWER in energy but carries exactly zero HF weight.
    assert coeffs["B1g"][0] < coeffs["Ag"][0], coeffs
    assert coeffs["B1g"][1] < 1e-20, coeffs["B1g"]
    # ...while the Ag ground state -- the one tol=1e-8 selects -- carries the real weight.
    assert coeffs["Ag"][1] == pytest.approx(0.4451, rel=0.05), coeffs["Ag"]


def test_G8_no_symmetric_offender_remains():
    """Was 'neither constant is safe': some side in the family had a SYMMETRIC-SCF lowest level whose
    residue cleared even the looser 1e-8 cut (a = 1.22: 1.26e-8 on Linux; chem-qhi), so no fixed
    constant separated physics from SCF residue. Every such offender was the free SCF stopped on the
    D2h saddle point. With the stability rule each candidate is the broken-symmetry minimum, whose
    lowest level is physical, so there is no lowest-level residue left for a constant to separate.
    Excited levels can still carry a conv_tol-dependent residue between the cuts (section 11); this
    gate does not cover them."""
    for a in (1.185, 1.190, 1.195, 1.22):
        _, amp2, _ = _amp2_and_spectrum(a)
        assert amp2[0] > 0.4, (a, float(amp2[0]))
        assert scf_symmetry_status(_square_h4(a))[0], a
