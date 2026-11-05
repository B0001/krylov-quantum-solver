"""Gates for specs/SPEC_rebuild_gauge_search.md, REVISED 2026-11-05 (chem-g14): _driver_orbitals
rebuilds through the builder's own SCF path (molecular_hamiltonian._driver_problem + map_canonical),
so the rebuild reproduces mh's operator bit-for-bit and the reproducibility check is strict again.
No SCF retry, no sign-gauge search: they existed to absorb independent-SCF-trajectory noise that a
shared, thread-pinned path cannot produce. The first version's G3-G5 tested that deleted machinery;
the spec's section 10 records it and why it went.

pyscf + qiskit-nature only (no block2): own process, per the test_*_spec.py convention.
"""
import dataclasses

import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp
from qiskit_nature.second_q.drivers import InitialGuess, PySCFDriver

import hybrid_quantum_solver.molecular_hamiltonian as mham
import reachability
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

H2 = "H 0 0 0; H 0 0 0.74"
LIN_H4 = "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"


def _sq(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


# The builder configurations test_eigenstate_reachability_spec.py hands to _driver_orbitals: its G2
# square-H4 grid, the G4 broken-symmetry witness (builder default conv_tol), a tight and a
# symmetric_scf build, an open shell, and both active-space shapes (int and (n_alpha, n_beta)).
SYSTEMS = [dict(atom=_sq(a), conv_tol=ct) for a in (1.05, 1.10, 1.19, 1.35) for ct in (1e-6, 1e-9)]
SYSTEMS += [dict(atom=_sq(1.20)), dict(atom=_sq(1.19), conv_tol=1e-13),
            dict(atom=_sq(1.19), symmetric_scf=True), dict(atom=H2), dict(atom=LIN_H4),
            dict(atom=LIN_H4, spin=2),
            dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
            dict(atom="O 0 0 0; H 0 0 0.97", spin=1, active_electrons=(3, 2), active_orbitals=4)]


def _residual(a, b):
    """Largest non-identity coefficient of a - b (the identity term is a frame shift)."""
    diff = (a - b).simplify(atol=0)
    non_identity = diff.paulis.x.any(axis=1) | diff.paulis.z.any(axis=1)
    return float(np.abs(diff.coeffs[non_identity]).max(initial=0.0))


def _count_scf(monkeypatch):
    calls = {"run_pyscf": 0, "to_problem": 0}
    for name in calls:
        real = getattr(PySCFDriver, name)

        def counted(self, *a, _real=real, _name=name, **kw):
            calls[_name] += 1
            return _real(self, *a, **kw)

        monkeypatch.setattr(PySCFDriver, name, counted)
    return calls


# --- G1: DEFINITION OF DONE -- the rebuild reproduces mh bit-for-bit --------------------------------

@pytest.mark.parametrize("kw", SYSTEMS, ids=lambda kw: ",".join(f"{k}={v}" for k, v in kw.items()))
def test_G1_rebuild_reproduces_the_operator_bit_for_bit(kw, monkeypatch):
    """With the reproducibility tolerance set to ZERO, every non-identity coefficient must match:
    the rebuild is the builder's computation, not a second SCF trajectory near it."""
    mh = build_molecular_hamiltonian(**kw)
    monkeypatch.setattr(reachability, "_SAME_OPERATOR_ATOL", 0.0)
    assert reachability._driver_orbitals(mh) is not None


# --- G2: the check is strict ----------------------------------------------------------------------

def test_G2_check_refuses_an_edit_smaller_than_scf_trajectory_noise():
    """1e-10 on one coefficient is below the smallest independent-SCF-trajectory residual measured
    (3.6e-10: the linear-H4 triplet at the default conv_tol, spec section 9). A tolerance loose
    enough to absorb a separate SCF run accepts this edit; the strict one refuses it. An identity
    shift is the same operator in a shifted frame and still verifies."""
    assert reachability._SAME_OPERATOR_ATOL <= 1e-12
    mh = build_molecular_hamiltonian(atom=_sq(1.35))
    n = mh.qubit_hamiltonian.num_qubits
    edited = dataclasses.replace(
        mh, qubit_hamiltonian=(mh.qubit_hamiltonian + SparsePauliOp("Z" + "I" * (n - 1), 1e-10))
        .simplify(atol=0))
    assert _residual(edited.qubit_hamiltonian, mh.qubit_hamiltonian) == pytest.approx(1e-10)
    assert reachability._driver_orbitals(edited) is None
    shifted = dataclasses.replace(
        mh, qubit_hamiltonian=(mh.qubit_hamiltonian + SparsePauliOp("I" * n, 0.37)).simplify())
    assert reachability._driver_orbitals(shifted) is not None


# --- G3: one SCF on every path ----------------------------------------------------------------------

def test_G3_match_and_refusal_each_cost_one_scf(monkeypatch):
    """Before this revision a refusal ran 5 SCFs and 5 x (2**(k-1) - 1) sign-flipped AO2MO
    transforms first (k = 4 here: 40 to_problem calls) -- work a deterministic rebuild cannot use."""
    fresh = build_molecular_hamiltonian(atom=_sq(1.19))
    stale = dataclasses.replace(
        fresh, qubit_hamiltonian=build_molecular_hamiltonian(atom=_sq(1.10)).qubit_hamiltonian)
    calls = _count_scf(monkeypatch)
    assert reachability._driver_orbitals(fresh) is not None
    assert calls == {"run_pyscf": 1, "to_problem": 1}, calls
    calls.update(run_pyscf=0, to_problem=0)
    assert reachability._driver_orbitals(stale) is None
    assert calls == {"run_pyscf": 1, "to_problem": 1}, calls


# --- G4: the rebuild follows the builder's SCF path -------------------------------------------------

def test_G4_rebuild_follows_a_change_to_the_builders_scf_path(monkeypatch):
    """Change the builder's SCF path (here: its initial guess, the kind of change a branch-selection
    fix such as chem-8ls might make) and the rebuild must still reproduce the operator exactly --
    because it calls the same path, not a copy of it. Non-vacuous: the changed path really does
    produce a different operator, so a rebuild through the old path would be refused."""
    atom = _sq(1.10)
    default = build_molecular_hamiltonian(atom=atom)

    class AtomGuess(mham.PySCFDriver):
        def __init__(self, *a, **kw):
            super().__init__(*a, init_guess=InitialGuess.ATOM, **kw)

    monkeypatch.setattr(mham, "PySCFDriver", AtomGuess)
    changed = build_molecular_hamiltonian(atom=atom)
    assert _residual(changed.qubit_hamiltonian, default.qubit_hamiltonian) \
        > reachability._SAME_OPERATOR_ATOL
    assert reachability._driver_orbitals(changed) is not None


# --- G5: build_args carries every SCF argument the builder took -------------------------------------

def test_G5_build_args_records_symmetric_scf():
    """The rebuild replays build_args through _driver_problem, so the builder must record every
    SCF argument _driver_problem takes -- symmetric_scf included, not just its floored conv_tol."""
    assert build_molecular_hamiltonian(atom=_sq(1.19), symmetric_scf=True) \
        .build_args["symmetric_scf"] is True
    assert build_molecular_hamiltonian(atom=H2).build_args["symmetric_scf"] is False
