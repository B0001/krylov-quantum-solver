"""Gates for specs/SPEC_rebuild_gauge_search.md (chem-g14): G3-G5, the retry/sign-search mechanics
inside reachability._driver_orbitals itself. G1/G2 (no regression / the chem-g14 failures clear)
are pinned, and revised, by tests/test_eigenstate_reachability_spec.py and the spec's own Status
section -- this file does not re-pin them."""
import dataclasses

from qiskit.quantum_info import SparsePauliOp

import reachability
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

H2 = "H 0 0 0; H 0 0 0.74"


def _sq(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


# --- G3: cheap when unneeded ------------------------------------------------------------------

def test_g3_zero_retry_zero_flip_match_costs_exactly_one_driver_call(monkeypatch):
    from qiskit_nature.second_q.drivers import PySCFDriver

    mh = build_molecular_hamiltonian(atom=H2)
    calls = []
    real_run_pyscf = PySCFDriver.run_pyscf

    def counting_run_pyscf(self, *a, **kw):
        calls.append(1)
        return real_run_pyscf(self, *a, **kw)

    monkeypatch.setattr(PySCFDriver, "run_pyscf", counting_run_pyscf)
    found = reachability._driver_orbitals(mh)
    assert found is not None
    assert len(calls) == 1, "a rebuild that already matches unflipped must not retry or search"


# --- G4: the energy pre-filter rejects a genuinely different SCF minimum ----------------------

def test_g4_energy_mismatch_triggers_a_retry_not_a_false_accept(monkeypatch):
    from qiskit_nature.second_q.drivers import PySCFDriver

    mh = build_molecular_hamiltonian(atom=H2)
    real_run_pyscf = PySCFDriver.run_pyscf
    calls = []

    def perturbing_run_pyscf(self, *a, **kw):
        calls.append(1)
        result = real_run_pyscf(self, *a, **kw)
        if len(calls) == 1:
            self._calc.e_tot += 0.01     # first attempt: a different, "worse" SCF minimum
        return result

    monkeypatch.setattr(PySCFDriver, "run_pyscf", perturbing_run_pyscf)
    found = reachability._driver_orbitals(mh)
    assert found is not None, "the genuine match on retry must still be found"
    assert len(calls) >= 2, "the energy-mismatched first attempt must not be accepted outright"


def test_g4_a_persistently_wrong_energy_exhausts_retries_and_refuses(monkeypatch):
    from qiskit_nature.second_q.drivers import PySCFDriver

    mh = build_molecular_hamiltonian(atom=H2)
    real_run_pyscf = PySCFDriver.run_pyscf

    def always_wrong_run_pyscf(self, *a, **kw):
        result = real_run_pyscf(self, *a, **kw)
        self._calc.e_tot += 0.01         # every attempt lands on the "wrong" minimum
        return result

    monkeypatch.setattr(PySCFDriver, "run_pyscf", always_wrong_run_pyscf)
    assert reachability._driver_orbitals(mh) is None


# --- G5: the sign search is honestly bounded ----------------------------------------------------

def test_g5_an_oversized_active_space_refuses_rather_than_searching_forever(monkeypatch):
    monkeypatch.setattr(reachability, "_MAX_SIGN_SEARCH_ORBITALS", 2)
    monkeypatch.setattr(reachability, "_MAX_REBUILD_RETRIES", 1)
    mh = build_molecular_hamiltonian(atom=_sq(1.1))        # 4 active orbitals: k - 1 = 3 > bound
    # an operator that can never match (same qubit count, but not the real Hamiltonian) forces the
    # unflipped comparison to fail and fall through to the sign-search bound check.
    wrong = SparsePauliOp("I" * mh.qubit_hamiltonian.num_qubits, 1.0)
    mismatched = dataclasses.replace(mh, qubit_hamiltonian=wrong)
    assert reachability._driver_orbitals(mismatched) is None


def test_g5_bound_is_inclusive_of_the_gated_systems(monkeypatch):
    """No system this module's own gate file exercises needs more than 6 active orbitals (module
    comment); confirm the default bound (10) does not clip any of them off by accident."""
    mh = build_molecular_hamiltonian(atom=_sq(1.1))
    k = mh.qubit_hamiltonian.num_qubits // 2
    assert k - 1 <= reachability._MAX_SIGN_SEARCH_ORBITALS
