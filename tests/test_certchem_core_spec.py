"""Spec gate for CertChem-M1 tasks 3-7: certified_energy over the validated pipeline.

Encodes specs/tasks/01-certchem-m1.md:
  Task 3 — H2/STO-3G returns a bracket containing the known FCI value.
  Task 4 — floor guard is the sole chokepoint: a sub-floor estimate raises
           FloorViolationError and no CertifiedResult is constructed.
  Task 5 — determinism: same inputs twice -> byte-identical serialized result.
  Task 6 — golden regression: H2, H4, LiH each (a) bracket contains FCI,
           (b) |estimate-FCI| <= tol, (c) floor passes.
  Task 7 — Mode.FAST returns a bare float, never a CertifiedResult.

Heavy (imports pyscf/qiskit); runs in the isolated spec process (Makefile test_*_spec glob).
"""

import dataclasses

import pytest

import certchem.core as core
from certchem import (
    CapExceededError,
    CertifiedResult,
    FloorViolationError,
    Mode,
    certified_energy,
    floor_guard,
)
from hybrid_quantum_solver import build_molecular_hamiltonian

# name -> (geometry, cas, tolerance in Ha). Small, exactly-referenceable systems.
GOLDEN = {
    "H2": ("H 0 0 0; H 0 0 0.735", (2, 2), 1e-6),
    "H4": ("H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7", (4, 4), 1e-6),
    "LiH": ("Li 0 0 0; H 0 0 1.6", (2, 5), 1e-4),
}


def _fci(mol, cas):
    mh = build_molecular_hamiltonian(
        atom=mol, basis="sto-3g", active_electrons=cas[0], active_orbitals=cas[1]
    )
    return mh.ground_state_energy()


# --- Task 3 & 6: golden regression gate --------------------------------------------------


@pytest.mark.parametrize("name", list(GOLDEN))
def test_golden_bracket_contains_fci_within_tol(name):
    mol, cas, tol = GOLDEN[name]
    result = certified_energy(mol, "sto-3g", cas)
    assert isinstance(result, CertifiedResult)
    fci = _fci(mol, cas)
    b = result.bracket
    # (a) bracket contains FCI
    assert b.lower_hartree <= fci <= b.upper_hartree, f"{name}: FCI {fci} outside {b}"
    # (b) estimate close to FCI
    assert abs(b.best_estimate_hartree - fci) <= tol, f"{name}: estimate off by > {tol}"
    # (c) floor passed (holding the object is the proof)
    assert result.certificate.floor_check == "pass"
    assert b.width >= 0.0


# --- Task 7: Mode.FAST ------------------------------------------------------------------


def test_fast_mode_returns_bare_float():
    mol, cas, _ = GOLDEN["H2"]
    out = certified_energy(mol, "sto-3g", cas, Mode.FAST)
    assert isinstance(out, float)
    assert not isinstance(out, CertifiedResult)
    # sanity: near FCI, not the double-offset garbage a wrong wiring would give
    assert abs(out - _fci(mol, cas)) < 1e-3


# --- Task 4: floor guard chokepoint -----------------------------------------------------


def test_floor_guard_unit_raises_below_floor():
    floor_guard(-1.0, -1.0)  # exactly at floor: ok
    floor_guard(-0.5, -1.0)  # above floor: ok
    with pytest.raises(FloorViolationError) as ei:
        floor_guard(-999.0, -1.1)
    assert ei.value.diagnostics["best_estimate_hartree"] == -999.0


def test_bad_estimate_raises_and_builds_no_result(monkeypatch):
    # Inject a broken estimator: best estimate far below a sane floor.
    monkeypatch.setattr(core, "_estimate", lambda mh, k, eps: (-999.0, -1.1, -1.0, 0.0))
    # Guard against a silent bypass: if a Bracket were ever built, fail loudly.
    monkeypatch.setattr(
        core,
        "Bracket",
        lambda *a, **k: pytest.fail("Bracket constructed despite floor violation"),
    )
    mol, cas, _ = GOLDEN["H2"]
    with pytest.raises(FloorViolationError):
        certified_energy(mol, "sto-3g", cas)


# --- Task 5: determinism ----------------------------------------------------------------


def test_determinism_byte_identical_serialization():
    mol, cas, _ = GOLDEN["H2"]
    r1 = certified_energy(mol, "sto-3g", cas)
    r2 = certified_energy(mol, "sto-3g", cas)
    assert dataclasses.asdict(r1) == dataclasses.asdict(r2)
    assert repr(r1) == repr(r2)


# --- caps propagate through the entry point ---------------------------------------------


def test_cap_exceeded_before_any_solve():
    with pytest.raises(CapExceededError) as ei:
        certified_energy("H 0 0 0; H 0 0 0.735", "cc-pVQZ", (2, 2))
    assert ei.value.cap == "basis"


# --- Task 8: worked example runs top-to-bottom ------------------------------------------


def test_n2_quickstart_example_runs():
    from certchem.examples import n2_quickstart

    n2_quickstart.main()  # raises if HF is not above the certified bracket, etc.
