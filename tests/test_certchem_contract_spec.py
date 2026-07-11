"""Spec gate for CertChem-M1 tasks 1-2: frozen contract types + envelope caps.

Encodes the acceptance criteria from specs/tasks/01-certchem-m1.md:
  Task 1 — contract types importable with ZERO solver deps; frozen; Bracket invariant;
           Mode.FAST distinct from CERTIFIED; the three exceptions exist.
  Task 2 — check_caps: one pass case, one fail per cap, error message names the cap.
"""

import dataclasses

import pytest

from certchem import (
    ALLOWED_BASES,
    MAX_SPIN_ORBITALS,
    Bracket,
    CapExceededError,
    Certificate,
    CertifiedResult,
    ConvergenceError,
    FloorViolationError,
    Mode,
    check_caps,
)


# --- Task 1: contract types --------------------------------------------------------------


def test_contract_imports_without_solver_deps():
    # The whole point of the contract layer: build/inspect results without pyscf/qiskit.
    # Must run in a FRESH interpreter — sys.modules is process-global and sibling spec tests
    # import pyscf, so an in-process sys.modules check would be meaningless.
    import subprocess
    import sys

    code = (
        "import certchem.contract, certchem.limits, sys;"
        "bad=[m for m in ('pyscf','qiskit','numpy','scipy') if m in sys.modules];"
        "assert not bad, bad; print('clean')"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, timeout=60
    )
    assert proc.returncode == 0, f"solver deps leaked into contract layer: {proc.stderr}"
    assert "clean" in proc.stdout


@pytest.mark.parametrize("cls", [Bracket, Certificate, CertifiedResult])
def test_contract_types_are_frozen(cls):
    assert dataclasses.is_dataclass(cls)
    assert cls.__dataclass_params__.frozen


def test_bracket_width_and_ordering_invariant():
    b = Bracket(lower_hartree=-1.5, upper_hartree=-1.0, best_estimate_hartree=-1.2)
    assert b.width == pytest.approx(0.5)
    # Invariant 2 enforced at construction: lower <= best <= upper.
    with pytest.raises(ValueError):
        Bracket(lower_hartree=-1.0, upper_hartree=-1.5, best_estimate_hartree=-1.2)
    with pytest.raises(ValueError):
        Bracket(lower_hartree=-1.5, upper_hartree=-1.0, best_estimate_hartree=-2.0)


def test_mode_fast_distinct_from_certified():
    assert Mode.FAST is not Mode.CERTIFIED
    assert Mode.FAST.value == "fast"
    assert Mode.CERTIFIED.value == "certified"


def test_exceptions_carry_their_vocabulary():
    with pytest.raises(FloorViolationError) as ei:
        raise FloorViolationError("below floor", diagnostics={"estimate": -999.0})
    assert ei.value.diagnostics["estimate"] == -999.0

    with pytest.raises(ConvergenceError) as ei:
        raise ConvergenceError("insufficient", partial={"krylov_dim": 3})
    assert ei.value.partial["krylov_dim"] == 3


def test_certified_result_composes():
    result = CertifiedResult(
        bracket=Bracket(-1.5, -1.0, -1.2),
        certificate=Certificate(
            method="temple_bound + variational_floor",
            floor_check="pass",
            krylov_dim=8,
            convergence="converged",
            solver_version="0.1.0",
            manifest=None,
        ),
    )
    assert result.certificate.floor_check == "pass"
    assert result.bracket.lower_hartree <= result.bracket.best_estimate_hartree


# --- Task 2: check_caps ------------------------------------------------------------------


def test_check_caps_pass_case():
    # H2 / sto-3g / CAS(2,2): 4 spin-orbitals, well inside the envelope.
    assert check_caps(molecule=object(), basis="STO-3G", cas=(2, 2)) is None


def test_check_caps_fail_basis_names_cap():
    with pytest.raises(CapExceededError) as ei:
        check_caps(molecule=object(), basis="cc-pVQZ", cas=(2, 2))
    assert ei.value.cap == "basis"
    assert "basis" in str(ei.value).lower()


def test_check_caps_fail_spin_orbitals_names_cap():
    # CAS(_, 9) -> 18 spin-orbitals > 16.
    with pytest.raises(CapExceededError) as ei:
        check_caps(molecule=object(), basis="sto-3g", cas=(10, 9))
    assert ei.value.cap == "spin_orbitals"
    assert str(MAX_SPIN_ORBITALS) in str(ei.value)


def test_check_caps_fail_electron_count_names_cap():
    with pytest.raises(CapExceededError) as ei:
        check_caps(molecule=object(), basis="sto-3g", cas=(99, 4))
    assert ei.value.cap == "electron_count"


def test_check_caps_fail_molecule_none_names_cap():
    with pytest.raises(CapExceededError) as ei:
        check_caps(molecule=None, basis="sto-3g", cas=(2, 2))
    assert ei.value.cap == "molecule"


def test_allowed_bases_are_lowercase_normalized():
    assert "sto-3g" in ALLOWED_BASES
    # case-insensitive acceptance
    assert check_caps(molecule=object(), basis="Sto-3g", cas=(2, 2)) is None
