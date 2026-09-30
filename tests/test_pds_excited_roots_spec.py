"""
Acceptance gates G1 (+G1b), G2-G4 for specs/SPEC_pds_excited_roots.md (higher roots of PDS's
P_K(E) as excited-state estimates).

Test-first: ``pds_roots`` did not exist before this bead (chem-pc1). ``pds_energy`` used to discard
every root but the smallest; the j-th smallest root of P_K(E) is a variational upper bound on the
j-th smallest *reachable* eigenvalue (Gauss-quadrature/Lanczos-Ritz equivalence) -- but this file's
own gates show that bound converges far too slowly to be *useful* at attainable K on H4, and that
naive index-matching on near-degenerate levels skips levels. Reference: dense diagonalization of the
same qubit Hamiltonian, filtered to nonzero-HF-overlap ("reachable") eigenvalues -- identical
construction to tests/test_qksd_excited_spec.py.

Thread pinning below (same pattern as tests/test_shift_both_sides_spec.py): at K=7-8 the moment
matrix M has cond(M) ~ 1e17-1e20, past double precision, and BLAS thread count selects between
discrete floating-point code paths there. Pinned to 1 thread, the equilibrium-ish H4 case is stable
over 8+ repeated fresh-process trials; the near-degenerate stretched H4 case is NOT -- even pinned,
K=7 there lands on one of two reproducible-per-path outcomes, one of which VIOLATES the variational
bound by ~7.5 mHa (not float noise -- see G1's docstring and specs/SPEC_pds_excited_roots.md R1/R3).
G1 is scoped to exclude that combination rather than assert something empirically false.

PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from moment_expansion import hamiltonian_moments, pds_roots

CHEM_ACC = 1.6e-3  # Ha (1 kcal/mol)


def _h4_equilibrium():
    """Same geometry as SPEC_moment_pds's H4 case -- 8 qubits."""
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")


def _h4_stretched():
    """Stretched H4 -- several reachable levels crowd within ~20 mHa of each other."""
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0")


def _reachable_spectrum(mh, overlap_tol=1e-8):
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    overlaps = np.abs(V.conj().T @ hf) ** 2
    reachable = w[overlaps > overlap_tol].real + mh.energy_offset
    return np.sort(reachable)


def test_G1_variational_bound_holds_for_every_root():
    """DEFINITION OF DONE (the can't-be-faked invariant): root_j >= reachable_j - 1e-9 Ha.

    On the equilibrium-ish H4, checked for K in 3..8 (verified reproducible over 8+ pinned-thread
    fresh-process trials, including K=7,8 where cond(M) ~ 1e17-1e20).

    On the near-degenerate stretched H4, checked only for K in 3..6 (cond(M) <~ 1.5e14, safely
    inside double precision). K=7-8 are DELIBERATELY EXCLUDED here, not because the bound is assumed
    to hold -- it measurably does NOT always hold there. Even under thread-pinning, K=7 on this
    geometry reproducibly lands on one of two floating-point code paths depending on process state:
    one gives a small positive margin (+196 uHa), the other VIOLATES the bound by -7.49 mHa (~4700x
    the 1e-9 Ha tolerance -- not roundoff dust). This is the honest finding: near-degenerate orbitals
    plus cond(M) > 1e17 can make PDS's variational guarantee fail outright in floating point, not
    just degrade to uselessness. See specs/SPEC_pds_excited_roots.md R3.
    """
    mh_eq = _h4_equilibrium()
    reachable_eq = _reachable_spectrum(mh_eq)
    mu_eq, off_eq = hamiltonian_moments(mh_eq, 16)
    for K in range(3, 9):
        roots = pds_roots(mu_eq, K, off_eq)
        n = min(len(roots), len(reachable_eq))
        for j in range(n):
            assert roots[j] >= reachable_eq[j] - 1e-9, (K, j, roots[j], reachable_eq[j])

    mh_st = _h4_stretched()
    reachable_st = _reachable_spectrum(mh_st)
    mu_st, off_st = hamiltonian_moments(mh_st, 16)
    for K in range(3, 7):
        roots = pds_roots(mu_st, K, off_st)
        n = min(len(roots), len(reachable_st))
        for j in range(n):
            assert roots[j] >= reachable_st[j] - 1e-9, (K, j, roots[j], reachable_st[j])


def test_G1b_bound_can_break_at_K7_on_near_degenerate_stretched_h4():
    """The excluded case from G1, made explicit and falsifiable rather than swept under the rug:
    cond(M) at K=7 on the stretched H4 is far past double precision (>1e14, the point past which
    G1 stops trusting the strict bound). This does not assert the bound fails on every run (it does
    not -- it is genuinely bimodal, see G1's docstring); it asserts the matrix is provably too
    ill-conditioned for the strict bound to be a reliable claim at this K on this geometry, which is
    the actual, falsifiable reason G1 excludes it."""
    mh = _h4_stretched()
    mu, off = hamiltonian_moments(mh, 16)
    K = 7
    M = np.array([[mu[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])
    cond = np.linalg.cond(M)
    assert cond > 1e14, cond


def test_G2_K8_first_excited_root_misses_chemical_accuracy():
    """KILLED: SPEC_moment_pds section 7 / the backlog scout probe predicted PDS(8)'s first excited
    root would land near reachable E_1. It does not -- the gap is 8-9x the 1.6 mHa chemical-accuracy
    bar (measured 12-14 mHa depending on process; see SPEC_pds_excited_roots.md R1). Asserted with a
    5 mHa margin well clear of that run-to-run float noise, so a future accidental "fix" that hides
    the conditioning problem instead of solving it gets caught by this gate flipping green, not by
    silence."""
    mh = _h4_equilibrium()
    reachable = _reachable_spectrum(mh)
    mu, off = hamiltonian_moments(mh, 16)
    roots8 = pds_roots(mu, 8, off)
    err_root1 = abs(roots8[1] - reachable[1])
    assert err_root1 > 5e-3, (err_root1, "usefulness claim unexpectedly holds -- update the spec")
    assert err_root1 > CHEM_ACC  # explicit: this is a documented miss of the accuracy bar


def test_G3_ground_root_converges_orders_of_magnitude_faster_than_first_excited():
    """Quantifies the backlog's '~1000x slower' claim: at fixed K=7, the ground root is >= 100x
    more converged than the first excited root, from the same linear solve."""
    mh = _h4_equilibrium()
    reachable = _reachable_spectrum(mh)
    mu, off = hamiltonian_moments(mh, 16)
    roots = pds_roots(mu, 7, off)
    err0 = abs(roots[0] - reachable[0])
    err1 = abs(roots[1] - reachable[1])
    assert err0 > 0
    assert err1 / err0 >= 100, (err0, err1, err1 / err0)


def test_G4_index_matching_skips_a_level_on_near_degenerate_stretched_h4():
    """Explicitly tests and reports the predicted failure mode: on near-degenerate stretched H4 at
    K=6 (chosen so cond(M) ~ 1.5e14, inside double precision, unlike K=7/8's cond(M) > 1e17),
    nearest-neighbor matching each of the first 6 roots to its closest reachable level is NOT
    injective -- a level gets skipped because a PDS root lands between two true levels. This is the
    informative finding (PDS roots track spectral density, not level index), not a bug to hide:
    the gate asserts the skip is reproducible, not that matching works."""
    mh = _h4_stretched()
    reachable = _reachable_spectrum(mh)
    mu, off = hamiltonian_moments(mh, 16)
    K = 6
    roots = pds_roots(mu, K, off)
    n_show = min(len(roots), 6)
    candidates = reachable[:8]

    # G1 still holds here -- the variational bound survives even where index-matching does not.
    for j in range(min(n_show, len(reachable))):
        assert roots[j] >= reachable[j] - 1e-9, (K, j, roots[j], reachable[j])

    assigned = [int(np.argmin(np.abs(candidates - roots[j]))) for j in range(n_show)]
    skipped = sorted(set(range(n_show)) - set(assigned))
    assert len(set(assigned)) < n_show, (
        "expected a non-injective (skipping) nearest-neighbor assignment on this near-degenerate "
        "geometry -- if this now passes, the index-skip finding needs re-verifying, not deleting",
        assigned,
    )
    assert skipped, (assigned, skipped)
