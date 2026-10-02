"""
Acceptance gates G0-G3 (+G2b) for specs/SPEC_centered_pds_shots.md (chem-u87): does centering
change how shot noise in <H^n> amplifies into PDS(K)?

Circuit-level noise: ONE random-Pauli classical-shadow record of |HF> feeds every moment, raw and
centered (symmetry-projected, see centered_pds_shots.py). The pre-registered run is 16 trials,
seeds 2000..2015, budgets 1e3..1e6 snapshots sampled directly. Thresholds were fixed before it.
G2 and G2b are KILLED and gated as measured (spec section 6).

The study (~1-2 min) runs in a spawned child with PYTHONHASHSEED=0 and BLAS pinned to 1 thread
(``pinned_n2_study``): raw-frame M at K>=6 is past 1/eps, so its last bits pick the result.
PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

from functools import lru_cache

import numpy as np
import scipy.linalg as sl
from qiskit.quantum_info import SparsePauliOp

import centered_pds_shots as cps
from classical_shadows import collect_classical_shadow, shadow_energy_samples
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from moment_expansion import hamiltonian_moments, pds_energy

TOL = 1.6e-3


@lru_cache(maxsize=None)
def _n2():
    return build_molecular_hamiltonian(**cps.N2)


@lru_cache(maxsize=None)
def _paths():
    return cps.pinned_n2_study(16, 2000)       # [base path, 8 term re-orderings]


def _study():
    return _paths()[0]


def _b(S):
    return cps.SAMPLED.index(S)


def test_G0a_noiseless_reproduces_centered_pds():
    mh = _n2()
    lam, V, hf, lo, hi = cps.sector_eigensystem(mh)
    d0 = cps.exact_weights(V, lo, hi, hf)
    ref, off = hamiltonian_moments(mh, 15)
    assert np.max(np.abs(cps.moments(lam, d0, 0.0, 15) - ref) / np.abs(ref)) < 1e-12
    res = _study()
    e0 = res["e0"]
    # sector E0 is the global minimum (scipy: numpy's complex eigh at 4096 is broken on this macOS)
    assert abs(sl.eigvalsh(mh.qubit_hamiltonian.to_matrix())[0] + off - e0) < 1e-9
    cen = cps.moments(lam, d0, res["mu"], 15)
    spec = {4: 1.076e-3, 5: 0.260e-3, 6: 0.0540e-3, 7: 0.00587e-3, 8: 0.00106e-3}
    for K, want in spec.items():
        got = pds_energy(cen, K, off + res["mu"]) - e0
        assert abs(got - want) < 0.01 * want, (K, got, want)


def test_G0b_same_estimator_as_classical_shadows_h2():
    """Sector-projected sum(lam*d) == classical_shadows' estimator of <P H P> on the same snapshots."""
    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    lam, V, hf, lo, hi = cps.sector_eigensystem(mh)
    n = mh.num_spatial_orbitals
    idx = (hi[:, None] << n | lo[None, :]).ravel()
    P = np.zeros((2 ** mh.num_qubits,) * 2)
    P[idx, idx] = 1.0
    op = SparsePauliOp.from_operator(P @ mh.qubit_hamiltonian.to_matrix() @ P)
    psi = np.asarray(mh.hf_state().data, dtype=complex)
    bases, signs = collect_classical_shadow(psi, mh.num_qubits, 3000, seed=7)
    ours = np.sum(lam * cps.shadow_sector_weights(bases, signs, V, lo, hi))
    theirs = shadow_energy_samples(bases, signs, op).mean()
    assert abs(ours - theirs) < 1e-10, (ours, theirs)


def test_G0c_fast_sampler_matches_collect_classical_shadow():
    mh = _n2()
    _, _, hf, _, _ = cps.sector_eigensystem(mh)
    psi = np.asarray(mh.hf_state().data, dtype=complex)
    for bases, signs in (collect_classical_shadow(psi, mh.num_qubits, 300, seed=3),
                         cps.sample_hf_shadow(hf, 300, np.random.default_rng(3))):
        z = bases == 2
        assert np.all(signs[z] == np.broadcast_to(1 - 2 * hf, signs.shape)[z])
        assert abs(signs[~z].mean()) < 4 / np.sqrt((~z).sum())


def test_G1_centering_does_not_change_noise_amplification():
    """Same data, shift-covariant functional: raw == centered wherever raw float64 is trustworthy,
    while the shot-noise error there (deviation from noiseless PDS(K)) is >= 100x that tolerance."""
    res = _study()
    for S in cps.SAMPLED:
        b = _b(S)
        raw, cen = res["raw"][b], res["cen"][b]
        ok = res["cond_raw"][b] < 1e10                       # (trials, K)
        assert ok[:, :3].any(axis=0).all(), S                # K=1..3 compared at every budget
        same = (np.abs(raw - cen) < 1e-6) | (np.isnan(raw) & np.isnan(cen))   # or no root in both
        assert same[ok].all(), (S, np.nanmax(np.abs(raw - cen)[ok]))
        dev = np.nan_to_num(np.abs(cen - res["exact"]), nan=np.inf)
        noise = [np.median(dev[ok[:, k], k]) for k in range(3)]
        assert min(noise) >= 1e-4, (S, noise)


def test_G2_verdict_rule_KILLED():
    """KILLED. Pre-registered (from the scout): survives at 1e6 and NOT at 1e3. Measured on the
    pinned path (PYTHONHASHSEED=0, 1 BLAS thread, macOS 27 Accelerate): survives at BOTH, each
    time at exactly the 0.25 margin at K=7 (1e6: centered 13/16 vs raw 9/16; 1e3: 5/16 vs 1/16).
    Across the 9 float64 paths the pair holds on 2. Raw's K=7 count moves with the BLAS build, so
    the asserted (portable) form of the kill is that the pair fails on some path."""
    pair = [cps.survives(r, _b(10**6)) and not cps.survives(r, _b(10**3)) for r in _paths()]
    assert not all(pair), pair


def test_G2b_verdict_depends_on_float64_path_KILLED():
    """KILLED. Registered at review: the G2 verdict is the same on all 9 float64 paths (H's own
    term order plus 8 re-orderings, same snapshots) at every sampled budget. Measured on macOS:
    paths surviving at 1e3 / 1e4 / 1e5 / 1e6 = 3 / 9 / 2 / 4 of 9. Centered hit counts are identical
    on every path. Raw's at K=7-8 are not (1e6, K=7: 8-13 of 16 vs centered 13), and the rule's
    4/16 margin sits inside that spread. Asserted in portable form: some budget splits, centered
    never moves."""
    paths = _paths()
    assert any(len({cps.survives(r, _b(S)) for r in paths}) == 2 for S in cps.SAMPLED)
    for S in cps.SAMPLED:
        cen = np.array([cps.hit_rate(r["cen"][_b(S)]) for r in paths])
        assert (cen == cen[0]).all(), S


def test_G3_independent_moment_noise_is_pessimistic():
    """Same marginals, independent per power: centered K=8 median error >= 10x the shared-shot one."""
    res = _study()
    b = _b(10**6)
    shared = np.median(np.nan_to_num(np.abs(res["cen"][b][:, 7]), nan=np.inf))
    indep = np.median(np.nan_to_num(np.abs(res["cen_ind"][b][:, 7]), nan=np.inf))
    assert indep >= 10 * shared, (indep, shared)
