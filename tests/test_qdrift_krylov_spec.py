"""Gates for specs/SPEC_qdrift_krylov.md -- qDRIFT real-time Krylov vs exact sector FCI.

All runs: H2O/STO-3G (O 0 0 0; H 0 +-0.757 0.586), Krylov dim M=6, R=20 realizations, seeds 0-2,
exact statevector simulation of the qDRIFT channel. Measured values are in the comments.
"""
import dataclasses
from functools import cache

import numpy as np
import pytest
import scipy.linalg
from qiskit.quantum_info import SparsePauliOp

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from qdrift_krylov import (
    check_reference,
    qdrift_krylov,
    qdrift_sample_count,
    sector_ground_energy,
)

GEO = "O 0 0 0; H 0 0.757 0.586; H 0 -0.757 0.586"
SEEDS = (0, 1, 2)
LS = (50, 400, 2000)
CHEM_ACC = 1.6e-3
CAS = [(2, 2), (4, 4)]


@cache
def _mh(ne, no, mu=0.0):
    mh = build_molecular_hamiltonian(GEO, "sto3g", active_electrons=ne, active_orbitals=no)
    if mu:
        n = mh.num_qubits  # H + mu * N_op: moves the Fock-space minimum into another sector
        N_op = SparsePauliOp.from_sparse_list(
            [("", [], n / 2)] + [("Z", [q], -0.5) for q in range(n)], n
        )
        mh = dataclasses.replace(mh, qubit_hamiltonian=(mh.qubit_hamiltonian + mu * N_op).simplify())
    return mh


@cache
def _runs(ne, no, L, postselect=False, mu=0.0):
    mh = _mh(ne, no, mu)
    return tuple(
        qdrift_krylov(mh, krylov_dim=6, samples_per_step=L, realizations=20,
                      postselect=postselect, seed=s)
        for s in SEEDS
    )


def _errs(ne, no, L, postselect=False, mu=0.0):
    ref = sector_ground_energy(_mh(ne, no, mu))
    return np.array([r.energy - ref for r in _runs(ne, no, L, postselect, mu)])


# G1 (definition of done) -- converges toward FCI as qDRIFT samples per step grow.
@pytest.mark.parametrize("ne,no,bound", [(2, 2, 3e-5), (4, 4, CHEM_ACC)])
def test_g1_converges_with_samples(ne, no, bound):
    med = [float(np.median(_errs(ne, no, L))) for L in LS]
    # measured CAS(2,2): 2.84e-4, 3.69e-5, 8.79e-6 ; CAS(4,4): 6.87e-3, 3.40e-3, 8.92e-4
    assert med[0] > med[1] > med[2], med
    assert med[2] < bound, med


# G2 -- realization-averaged (H, S) is a Gram pair: never below the floor (exact-sector FCI here).
@pytest.mark.parametrize("ne,no", CAS)
@pytest.mark.parametrize("postselect", [False, True])
def test_g2_variational_floor(ne, no, postselect):
    for L in LS:
        assert np.all(_errs(ne, no, L, postselect) >= -1e-9), (L, _errs(ne, no, L, postselect))


# G3 -- boundary: the samples needed for chemical accuracy grow with lambda.
def test_g3_lambda_boundary():
    assert 4 * _runs(2, 2, 400)[0].lam < _runs(4, 4, 400)[0].lam  # 2.44 vs 11.22
    assert np.median(_errs(2, 2, 400)) < CHEM_ACC          # 3.7e-5: CAS(2,2) is easy at L=400
    assert np.median(_errs(4, 4, 400)) > CHEM_ACC          # 3.4e-3: CAS(4,4) is not


# G4 -- the finding: single Pauli rotations leak electron number; the unconditional floor is the
# FOCK-space minimum, and postselecting onto the N-electron sector restores the N-electron floor.
def test_g4_particle_number_leak_and_floor():
    assert max(r.max_leakage for r in _runs(2, 2, 50)) < 1e-12        # measured 0
    leak50 = max(r.max_leakage for r in _runs(4, 4, 50))               # measured 0.15
    leak2k = max(r.max_leakage for r in _runs(4, 4, 2000))             # measured 5.3e-3
    assert leak50 > 0.05 and leak2k < leak50
    mh = _mh(4, 4, -2.0)
    assert mh.ground_state_energy() < sector_ground_energy(mh) - 1.0   # Fock min 2.4 Ha lower
    raw = _errs(4, 4, 400, False, -2.0)                                 # -10.0, -19.7, -42.1 mHa
    assert np.all(raw < -5e-3), raw                                     # BELOW sector FCI
    ps = _errs(4, 4, 400, True, -2.0)                                   # +3.4, +3.3, +3.2 mHa
    assert np.all(ps >= -1e-9), ps
    # postselection never hurts in the physical (mu=0) cases: measured ~10-25% lower error
    assert np.median(_errs(4, 4, 2000, True)) <= np.median(_errs(4, 4, 2000))


# G5 -- the pasted template's defects are caught, not silently computed.
def test_g5_template_defects_detected():
    mh = _mh(2, 2)
    vac = np.zeros(2**mh.num_qubits, complex)
    vac[0] = 1.0
    with pytest.raises(ValueError, match="sector"):
        check_reference(mh, vac)                                        # vacuum start refused
    with pytest.raises(ValueError, match="sector"):
        qdrift_krylov(mh, psi0=vac, seed=0)
    check_reference(mh, mh.hf_state().data)                             # HF accepted
    # template basis: full-H evolution of the vacuum gives phase copies -> rank-1 S
    Hd = mh.qubit_hamiltonian.to_matrix()
    states = [vac] + [scipy.linalg.expm(-1j * a * Hd) @ vac for a in (0.1, -0.1, 0.1)]
    V = np.array(states).T
    S, H = V.conj().T @ V, V.conj().T @ Hd @ V
    assert np.sum(np.linalg.eigvalsh(S) > 1e-8) == 1
    with pytest.raises(np.linalg.LinAlgError):
        scipy.linalg.eigh(H, S)                                          # unthresholded solve dies


# G6 -- qDRIFT cost is lambda^2 t^2 / eps, not lambda^2 / eps^2.
def test_g6_sample_count_scaling():
    n = qdrift_sample_count(10.0, 1.0, 1e-3)
    assert n == 200_000
    assert qdrift_sample_count(10.0, 2.0, 1e-3) == 4 * n
    assert qdrift_sample_count(10.0, 1.0, 5e-4) == 2 * n
