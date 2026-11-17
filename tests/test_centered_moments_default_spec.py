"""
Acceptance gates G1-G5 for specs/SPEC_centered_moments_default.md (chem-nxn): does centering
``moment_expansion.hamiltonian_moments`` by mu=<H> by default -- with a conditioning-guarded
fallback to raw -- unlock K=8 PDS with no caller-side centering, while leaving every existing
closed-spec gate (SPEC_moment_pds, SPEC_pds_excited_roots, SPEC_centered_pds) numerically unchanged?

Thread pinning (same pattern as tests/test_centered_pds_spec.py): raw-frame M at K>=6 is past 1/eps
and BLAS thread count picks the floating-point path.

PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

from functools import lru_cache

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from moment_expansion import cmx2_energy, hamiltonian_moments, pds_energy

CHEM_ACC = 1.6e-3
MAX_K = 8

GEOMETRIES = {
    "H2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "H4": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    "N2": dict(atom="N 0 0 0; N 0 0 1.10", basis="sto3g", active_electrons=6, active_orbitals=6),
}
DEEP = ("H4", "LiH", "N2")  # > 8 HF-reachable eigenstates: K=8 does not exhaust the Krylov space


def _cond(moments, K):
    M = np.array([[moments[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])
    return np.linalg.cond(M)


@lru_cache(maxsize=None)
def _mh(name):
    return build_molecular_hamiltonian(**GEOMETRIES[name])


def test_G1_moment_pds_gates_reproduced_exactly():
    """SPEC_moment_pds G1-G4, re-run verbatim through the new default (no `centered` flag) --
    must reproduce the same pass/fail pattern as tests/test_moment_pds_spec.py."""
    cases = {n: _mh(n) for n in ("H2", "H4", "LiH")}

    for mh in cases.values():  # G1: PDS(1) == <H>
        mu, off = hamiltonian_moments(mh, 2)
        assert abs(pds_energy(mu, 1, off) - (mu[1] + off)) < 1e-9

    for mh in cases.values():  # G2: variational + chemical accuracy at K=4
        fci = mh.ground_state_energy()
        mu, off = hamiltonian_moments(mh, 7)
        for K in (1, 2, 3, 4):
            e = pds_energy(mu, K, off)
            assert e >= fci - 1e-9, (mh.num_qubits, K, e, fci)
        assert abs(pds_energy(mu, 4, off) - fci) < CHEM_ACC

    for mh in cases.values():  # G3: monotone tightening
        mu, off = hamiltonian_moments(mh, 7)
        es = [pds_energy(mu, K, off) for K in (1, 2, 3, 4)]
        assert all(es[k + 1] <= es[k] + 1e-9 for k in range(3)), es

    mh = cases["H2"]  # G4: CMX non-variational boundary
    fci = mh.ground_state_energy()
    mu, off = hamiltonian_moments(mh, 4)
    assert cmx2_energy(mu, off) < fci - 1e-4
    assert pds_energy(mu, 2, off) >= fci - 1e-9


def test_G2_K8_variational_monotone_with_no_caller_centering():
    """THE UNLOCK: pds_energy(*hamiltonian_moments(mh, 2*K), K) -- the plain default call, no
    centered_frame, no explicit flag -- is variational and monotone to K=8 and beats
    SPEC_centered_pds's own G3 per-system thresholds."""
    tight = {"H4": 1e-6, "LiH": CHEM_ACC, "N2": 5e-6}
    for name in DEEP:
        mh = _mh(name)
        fci = mh.ground_state_energy()
        mu, off = hamiltonian_moments(mh, 2 * MAX_K)
        es = [pds_energy(mu, K, off) for K in range(1, MAX_K + 1)]
        assert all(e >= fci - 1e-9 for e in es), (name, [e - fci for e in es])
        assert all(es[k + 1] <= es[k] + 1e-9 for k in range(MAX_K - 1)), (name, es)
        assert es[MAX_K - 1] - fci < tight[name], (name, es[MAX_K - 1] - fci)
        assert es[MAX_K - 1] - fci < 0.5 * (es[3] - fci), name


def test_G3_default_actually_centers_on_deep_systems():
    """The fallback picked the centered branch (not a coincidental raw pass): moments[1] ~ 0 and
    offset absorbs <H>."""
    for name in DEEP:
        mh = _mh(name)
        mu, off = hamiltonian_moments(mh, 2 * MAX_K)
        mu_raw, off_raw = hamiltonian_moments(mh, 2 * MAX_K, centered=False)
        assert abs(mu[1]) < 1e-8, (name, mu[1])
        assert abs(off - (off_raw + mu_raw[1])) < 1e-8, (name, off, off_raw, mu_raw[1])


def test_G4_fallback_engages_on_h2_boundary_without_regressing_it():
    """H2 at max_order=7 (SPEC_moment_pds G2's own call): past K=n_reachable=2, centering makes
    the Hankel matrix at the deepest supported K (4) MORE ill-conditioned than raw, so the fallback
    returns the exact raw moments -- reproducing SPEC_moment_pds's original numbers bit-for-bit."""
    mh = _mh("H2")
    mu_default, off_default = hamiltonian_moments(mh, 7)
    mu_raw, off_raw = hamiltonian_moments(mh, 7, centered=False)
    assert np.allclose(mu_default, mu_raw)
    assert off_default == off_raw

    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    pops = np.abs(V.conj().T @ np.asarray(mh.hf_state().data)) ** 2
    assert int(np.sum(pops > 1e-8)) == 2  # the n_reachable=2 boundary from SPEC_centered_pds G2b
    cond_raw_k4 = _cond(mu_raw, 4)
    assert cond_raw_k4 > 1e14  # already past the double-precision cliff


def test_G5_centered_false_reproduces_other_closed_specs_raw_numbers():
    """centered=False is the exact raw computation SPEC_pds_excited_roots and SPEC_centered_pds's
    own gates depend on, not an approximation -- reproduce their specific recorded numbers."""
    mh_st = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0")
    mu, off = hamiltonian_moments(mh_st, 16, centered=False)
    assert _cond(mu, 7) > 1e14  # SPEC_pds_excited_roots G1b

    mh_n2 = _mh("N2")
    mu_n2, off_n2 = hamiltonian_moments(mh_n2, 2 * MAX_K, centered=False)
    assert all(_cond(mu_n2, K) > 1e17 for K in (6, 7, 8))  # SPEC_centered_pds G4
