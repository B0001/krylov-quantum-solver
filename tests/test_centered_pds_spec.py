"""
Acceptance gates G1-G5 for specs/SPEC_centered_pds.md (chem-70c): is PDS's high-K "det M -> 0"
accuracy floor an artifact of the raw electronic moment frame?

PDS(K) is shift-covariant in exact arithmetic: building the moments of H - mu*I and adding mu back
gives the identical energy. In float64 it is not: raw electronic-frame moments grow like ||H||^n, so
the Hankel moment matrix M inherits a dynamic range of ~||H||^(2K-2). These gates check that the
centered frame (``device_odmd.centered_frame``) is the same functional where raw is trustworthy
(G1), that it removes >= 8 orders of cond(M) at K=8 (G2) -- except where M is *genuinely* singular
because K exceeds the number of HF-reachable eigenstates (G2b, the boundary) -- that it unlocks
variational, monotone PDS at K=7-8 (G3), that raw really does break there (G4), and that the
centering shift need not be known accurately: mu = <H> (free: it is mu_1) works (G5).

Reuse only: hamiltonian_moments, pds_energy, device_odmd.centered_frame. ``_shifted`` below is the
same dataclasses.replace shift centered_frame performs, at an arbitrary mu (test-only, for G5).

Thread pinning (same as tests/test_pds_excited_roots_spec.py): raw-frame M at K>=6 is past 1/eps
and BLAS thread count picks the floating-point path. Centered results were identical across 1/2/4
threads; raw results were not (recorded in the spec, which is why G4 asserts only the robust part).

PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import dataclasses
from functools import lru_cache

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from device_odmd import centered_frame
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from moment_expansion import hamiltonian_moments, pds_energy

CHEM_ACC = 1.6e-3
MAX_K = 8

GEOMETRIES = {
    "H2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "H4": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    "N2": dict(atom="N 0 0 0; N 0 0 1.10", basis="sto3g", active_electrons=6, active_orbitals=6),
}
DEEP = ("H4", "LiH", "N2")  # > 8 HF-reachable eigenstates: K=8 does not exhaust the Krylov space


def _shifted(mh, mu):
    """H -> H - mu*I with mu moved into energy_offset (centered_frame's shift, at an arbitrary mu)."""
    shifted = (mh.qubit_hamiltonian - SparsePauliOp("I" * mh.num_qubits, coeffs=[mu])).simplify()
    return dataclasses.replace(mh, qubit_hamiltonian=shifted, energy_offset=mh.energy_offset + mu)


def _cond(moments, K):
    M = np.array([[moments[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])
    return np.linalg.cond(M)


@lru_cache(maxsize=None)
def _system(name):
    """(mh, fci, mu_center, raw (moments, offset), centered (moments, offset))."""
    mh = build_molecular_hamiltonian(**GEOMETRIES[name])
    mh_c, _, mu_c = centered_frame(mh)
    return (mh, mh.ground_state_energy(), mu_c,
            hamiltonian_moments(mh, 2 * MAX_K), hamiltonian_moments(mh_c, 2 * MAX_K))


def test_G1_centered_equals_raw_where_raw_is_trustworthy():
    """Same functional: |PDS_cent - PDS_raw| < 1e-6 Ha at every K where raw cond(M) < 1e10."""
    for name in GEOMETRIES:
        _, _, _, (raw, off), (cen, off_c) = _system(name)
        checked = 0
        for K in range(1, MAX_K + 1):
            if _cond(raw, K) < 1e10:
                checked += 1
                diff = abs(pds_energy(cen, K, off_c) - pds_energy(raw, K, off))
                assert diff < 1e-6, (name, K, diff)
        assert checked >= 2, (name, checked)   # the comparison is never vacuous


def test_G2_centering_cuts_cond_by_8_orders_at_K8():
    """THE CLAIM: log10 cond_raw(M) - log10 cond_cent(M) >= 8 at K=8 (measured 13.0 / 11.3 / 17.4)."""
    for name in DEEP:
        _, _, _, (raw, _), (cen, _) = _system(name)
        drop = np.log10(_cond(raw, MAX_K)) - np.log10(_cond(cen, MAX_K))
        assert drop >= 8.0, (name, drop)


def test_G2b_boundary_krylov_exhaustion_is_fundamental():
    """THE BOUNDARY: H2 has only 2 HF-reachable eigenstates, so M is exactly singular for K >= 3 in
    EVERY frame (det M -> 0 here is real, not a frame artifact) and centering cannot rescue it."""
    mh, fci, _, (raw, off), (cen, off_c) = _system("H2")
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    pops = np.abs(V.conj().T @ np.asarray(mh.hf_state().data)) ** 2
    assert int(np.sum(pops > 1e-8)) == 2
    assert abs(pds_energy(cen, 2, off_c) - fci) < 1e-9     # already exact at K = n_reachable
    assert _cond(cen, 3) > 1e14                            # singular past it, centered or not
    assert np.log10(_cond(raw, MAX_K)) - np.log10(_cond(cen, MAX_K)) < 8.0


def test_G3_centered_pds_variational_and_monotone_to_K8():
    """The unlock: in the centered frame PDS(K) >= FCI and PDS(K+1) <= PDS(K) for K = 1..8, and
    K=8 is far tighter than the old K=4 cap (measured H4 0.06 uHa, LiH 0.091 mHa, N2 1.06 uHa)."""
    tight = {"H4": 1e-6, "LiH": CHEM_ACC, "N2": 5e-6}
    for name in DEEP:
        _, fci, _, _, (cen, off_c) = _system(name)
        es = [pds_energy(cen, K, off_c) for K in range(1, MAX_K + 1)]
        assert all(e >= fci - 1e-9 for e in es), (name, [e - fci for e in es])
        assert all(es[k + 1] <= es[k] + 1e-9 for k in range(MAX_K - 1)), (name, es)
        assert es[MAX_K - 1] - fci < tight[name], (name, es[MAX_K - 1] - fci)
        assert es[MAX_K - 1] - fci < 0.5 * (es[3] - fci), name   # K=8 beats the K=4 cap


def test_G4_raw_frame_breaks_past_the_cap_on_n2():
    """Reproduces the raw-frame failure on N2 CAS(6,6): cond(M) is past 1/eps from K=6 on, and raw
    PDS departs from the (variational, monotone) centered value by > 1e-5 Ha somewhere in K=6..8.
    HOW it breaks is BLAS-path dependent (non-monotone K=6->7 at 1-2 threads; hundreds of Ha below
    FCI at K=7 or 8 at 4 threads) -- see the spec; only the path-independent part is gated."""
    _, _, _, (raw, off), (cen, off_c) = _system("N2")
    assert all(_cond(raw, K) > 1e17 for K in (6, 7, 8))
    devs = [abs(pds_energy(raw, K, off) - pds_energy(cen, K, off_c)) for K in (6, 7, 8)]
    assert max(devs) > 1e-5, devs


def test_G5_insensitive_to_a_badly_estimated_center():
    """mu from dense diagonalization is a validation-scale cheat; it need not be accurate.
    (a) N2: shifting by mu_center + delta, |delta| <= 2 Ha (~0.84 of the reachable half-width),
        leaves PDS(7), PDS(8) unchanged to < 1e-6 Ha.
    (b) mu = <H> (the first moment itself: free, no diagonalization) reproduces centered PDS at every
        K <= 8 to < 1e-6 Ha and still cuts cond(M) at K=8 by >= 8 orders (measured 9.5 / 12.2 / 14.1)."""
    mh, _, mu_c, _, (cen, off_c) = _system("N2")
    for delta in (-2.0, -1.0, -0.5, 0.5, 1.0, 2.0):
        mom, off = hamiltonian_moments(_shifted(mh, mu_c + delta), 2 * MAX_K)
        for K in (7, 8):
            assert abs(pds_energy(mom, K, off) - pds_energy(cen, K, off_c)) < 1e-6, (delta, K)

    for name in DEEP:
        mh, fci, _, (raw, _), (cen, off_c) = _system(name)
        mom, off = hamiltonian_moments(_shifted(mh, raw[1]), 2 * MAX_K)
        for K in range(1, MAX_K + 1):
            assert abs(pds_energy(mom, K, off) - pds_energy(cen, K, off_c)) < 1e-6, (name, K)
        assert np.log10(_cond(raw, MAX_K)) - np.log10(_cond(mom, MAX_K)) >= 8.0, name
