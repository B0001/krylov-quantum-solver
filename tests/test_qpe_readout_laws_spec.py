"""
Acceptance gates G1-G4 for specs/SPEC_qpe_readout_laws.md (QPE readout -- the precision law is a
staircase, and the state-prep law is overlap-independent).

Deliberately no new library code: the checks live entirely here, reusing `qpe_walk_readout.py`'s
existing `run_qpe`/`hartree_fock_vector` unmodified -- a genuine external verification of behavior
the module already exhibits (in its `__main__` printout) but never gates.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from qpe_walk_readout import run_qpe

T_RANGE = list(range(4, 16))


@pytest.fixture(scope="module")
def h2_system():
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, 2, 2)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), 2)
    from qubitization_blueprint import build_qubit_hamiltonian
    H, _n = build_qubit_hamiltonian(h1, eri, 2)
    _, Vk = np.linalg.eigh(H)
    ground = Vk[:, 0]
    return h1, eri, e_core, float(cas.e_tot), ground


def _precision_curve(h2_system):
    h1, eri, e_core, casci, ground = h2_system
    errs, lam = [], None
    for t in T_RANGE:
        E_est, lam, _ps, _olap = run_qpe(h1, eri, 2, e_core, ground, t)
        errs.append(abs(E_est - casci))
    return errs, lam


def test_G1_precision_error_monotonically_non_increasing(h2_system):
    """Adding phase bits never hurts the argmax point estimate -- err(t) is monotonically
    non-increasing, the true (staircase) shape of dyadic-grid refinement."""
    errs, _lam = _precision_curve(h2_system)
    for e_next, e_prev in zip(errs[1:], errs[:-1]):
        assert e_next <= e_prev + 1e-15, (errs)


def test_G2_precision_bounded_by_constant_times_lambda_over_2t(h2_system):
    """A real, provable upper bound -- err(t) <= 3 * lambda/2^t at every tested t -- replacing the
    docstring's bare '~lambda/2^t' with an actually-checked bound."""
    h1, eri, e_core, casci, ground = h2_system
    for t in T_RANGE:
        E_est, lam, _ps, _olap = run_qpe(h1, eri, 2, e_core, ground, t)
        err = abs(E_est - casci)
        assert err <= 3.0 * lam / 2 ** t, (t, err, lam / 2 ** t)


def _overlap_sweep_ratios(h2_system, t, seed=0):
    h1, eri, e_core, _casci, ground = h2_system
    rng = np.random.default_rng(seed)
    dim = ground.shape[0]
    ratios = []
    for p in (0.99, 0.9, 0.7, 0.5, 0.3, 0.1, 0.05, 0.02):
        r = rng.normal(size=dim) + 1j * rng.normal(size=dim)
        r = r - np.vdot(ground, r) * ground
        r = r / np.linalg.norm(r)
        trial = np.sqrt(p) * ground + np.sqrt(1 - p) * r
        trial = trial / np.linalg.norm(trial)
        _E_est, _lam, ps, olap = run_qpe(h1, eri, 2, e_core, trial, t)
        ratios.append(ps / olap)
    return ratios


@pytest.mark.parametrize("t", (8, 10, 12))
def test_G3_success_overlap_ratio_is_overlap_independent(h2_system, t):
    """THE FINDING / definition of done: across an overlap sweep spanning a 50x range (0.02 to
    0.99), the band width of p_success(window)/measured_overlap at FIXED t is < 0.05 -- the
    state-prep bottleneck's t-dependent prefactor does not depend on the overlap value itself."""
    ratios = _overlap_sweep_ratios(h2_system, t)
    assert max(ratios) - min(ratios) < 0.05, (t, ratios)


@pytest.mark.parametrize("t", (8, 10, 12))
def test_G4_ratio_bounded_but_not_monotonic_in_t(h2_system, t):
    """Boundary, recorded not smoothed over: the representative ratio stays in a sane regime
    [0.85, 1.05] at every tested t, but (checked across the three t values together, see the
    cross-check below) it is NOT monotonically increasing in t -- t=10 is measurably lower than
    t=8, proving G3's overlap-independence does not extend to t-independence or monotonicity."""
    ratios = _overlap_sweep_ratios(h2_system, t)
    assert all(0.85 <= r <= 1.05 for r in ratios), (t, ratios)


def test_G4_cross_check_ratio_is_not_monotonic_across_t(h2_system):
    """The explicit non-monotonicity claim: t=10's ratio is lower than t=8's, so the sequence
    across t=8,10,12 is not monotonically increasing."""
    r8 = min(_overlap_sweep_ratios(h2_system, 8))
    r10 = min(_overlap_sweep_ratios(h2_system, 10))
    r12 = min(_overlap_sweep_ratios(h2_system, 12))
    assert r10 < r8, (r8, r10, r12)
    assert r12 > r10, (r8, r10, r12)
