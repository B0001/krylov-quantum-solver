"""
Acceptance gates G1-G4 for specs/SPEC_odmd.md (observable dynamic mode decomposition).

Test-first: ``odmd`` does not exist yet, so this file is RED until the spec is implemented. ODMD
estimates the ground-state energy from ONLY the complex survival amplitude
s_k = <phi_0|e^{-ik tau H}|phi_0> -- the first row of the overlap matrix QKSD already measures --
via Hankel-matrix DMD with an SVD-truncated pseudoinverse (arXiv:2306.01858). No Hamiltonian
matrix element is ever measured, so the lambda-scaled sampling noise that dominates KQD vanishes;
at a matched count of measured elements and shots per element, ODMD should beat KQD.

The KQD comparison arm is the validated one from msd.py (same centered frame, same per-element
Hadamard-test noise scales, pinned by tests/test_msd_sampling_spec.py). Errors use the MEDIAN over
noise seeds. PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np
import pytest
from pyscf import lib

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from msd import build_msd_problem, sample_ground_energy
from odmd import build_odmd_problem, odmd_energy, sample_odmd_energy
from reachability import TIGHT_SCF_CONV_TOL, _dense_hf_projection, reachable_mask

SYSTEMS = {
    "h2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "h4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "n2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
_CACHE = {}


def _case(key):
    if key not in _CACHE:
        mh = build_molecular_hamiltonian(**SYSTEMS[key])
        _CACHE[key] = (mh, build_odmd_problem(mh, n=20))
    return _CACHE[key]


def _median_abs_err(prob, shots, seeds=100, n=None, svd_threshold=None):
    errs = [abs(sample_odmd_energy(prob, shots, sd, n=n, svd_threshold=svd_threshold) - prob.ref)
            for sd in range(seeds)]
    return float(np.median(errs))


def test_G1_overlaps_alone_recover_fci():
    """Noiseless ODMD at K=20 matches FCI < 1e-5 Ha on H2, H4 and N2 CAS(6,6)."""
    for key in SYSTEMS:
        mh, prob = _case(key)
        energy, _ = odmd_energy(prob.s, prob.tau)
        assert abs(energy + prob.offset - mh.ground_state_energy()) < 1e-5, key


def test_G2_depth_convergence_and_no_variational_floor():
    """Error shrinks with depth on N2 -- and the K=8 estimate falls BELOW FCI (non-variational,
    unlike QKSD's Ritz values; the recorded boundary)."""
    _, prob = _case("n2")
    e8, _ = odmd_energy(prob.s[:8], prob.tau)
    e16, _ = odmd_energy(prob.s[:16], prob.tau)
    assert abs(e16 - prob.ref) < abs(e8 - prob.ref), (e8 - prob.ref, e16 - prob.ref)
    assert abs(e16 - prob.ref) < 1e-5, e16 - prob.ref
    assert e8 < prob.ref - 1e-4, e8 - prob.ref            # dips below the exact ground state


def test_G3_matched_budget_beats_kqd():
    """DEFINITION OF DONE: N2 CAS(6,6), matched measured-element count (ODMD K=16 overlaps vs
    KQD n=8 overlaps + 8 H elements), same shots per element: ODMD median error wins big."""
    mh, prob = _case("n2")
    kqd = build_msd_problem(mh, n=8, order=8, delta=0.38)
    assert abs(kqd.tau - prob.tau) < 1e-12                # same centered frame
    for shots, min_ratio in ((10_000, 15.0), (100_000, 4.0)):
        m_odmd = _median_abs_err(prob, shots, n=16)
        m_kqd = float(np.median([abs(sample_ground_energy(kqd, shots, "kqd", sd) - prob.ref)
                                 for sd in range(100)]))
        assert m_odmd < m_kqd, (shots, m_odmd, m_kqd)
        assert m_kqd / m_odmd > min_ratio, (shots, m_kqd / m_odmd)


def test_G4_svd_truncation_is_the_mechanism():
    """At 1e4 shots the noise-aware truncation keeps ODMD < 5 mHa median; removing it
    (threshold 1e-10) inflates the median by > 20x."""
    _, prob = _case("n2")
    m_on = _median_abs_err(prob, 10_000, n=16)
    m_off = _median_abs_err(prob, 10_000, n=16, svd_threshold=1e-10)
    assert m_on < 5e-3, m_on
    assert m_off / m_on > 20.0, (m_off, m_on)


def test_G5_ref_is_hf_reachable_not_global_fci_ground():
    """chem-1x3: none of G1-G4's systems (H2, H4 chain, N2 CAS(6,6)) separate the global FCI
    ground from the HF-reachable one -- HF reaches the true ground on all three, so
    ``build_odmd_problem``'s ``ref = reach[0] - mu`` (chem-3ds's fix, from the earlier
    ``w_eig[0] - mu``) is unexercised by every existing gate. Square H4 at a = 1.19 Angstrom with
    a symmetric SCF (``TIGHT_SCF_CONV_TOL``) is a witness where it is not: HF's determinant has
    (numerically) zero population on the global ground (symmetry-forbidden, same family as
    tests/test_eigenstate_reachability_spec.py's G9 control), so the survival-amplitude signal
    s_k = <phi_0|e^{-ik tau H}|phi_0> has no component at that energy and can never converge to
    it -- `ref` must be the lower bound the signal actually targets, the HF-reachable level.

    chem-7cx: the free-RHF branch on this near-degenerate geometry is nondeterministic across
    ambient OpenMP thread counts; pin pyscf to 1 thread for the reproducible residue (same fix as
    benchmark_nbn.py's `_tight_scf`, tests/test_nbn_scf_determinism_spec.py).
    """
    atom = "H 0 0 0; H 1.19 0 0; H 1.19 1.19 0; H 0 1.19 0"
    ambient = lib.num_threads()
    lib.num_threads(1)
    try:
        mh = build_molecular_hamiltonian(atom=atom, conv_tol=TIGHT_SCF_CONV_TOL)
    finally:
        lib.num_threads(ambient)

    w_eig, V, pops = _dense_hf_projection(mh)
    reach = w_eig[reachable_mask(mh, w_eig, V, pops, 1e-8)].real
    assert pops[0] < 1e-20, pops[0]                          # global ground: HF-unreachable
    assert reach.min() - w_eig[0] > 1e-2, reach.min() - w_eig[0]  # ...and it's a DIFFERENT level

    prob = build_odmd_problem(mh, n=20)
    assert prob.ref == pytest.approx(reach.min() - prob.mu, abs=1e-8)    # the fix
    assert prob.ref != pytest.approx(w_eig[0] - prob.mu, abs=1e-3)       # not the old formula

    # and the signal itself backs this: noiseless ODMD converges to ref, not the unreachable
    # global ground -- the fix matches what the observable can actually see, not just a relabeling.
    energy, _ = odmd_energy(prob.s, prob.tau)
    assert abs(energy - prob.ref) < 1e-5, (energy, prob.ref)
    assert abs(energy - (w_eig[0] - prob.mu)) > 1e-2, (energy, w_eig[0] - prob.mu)
