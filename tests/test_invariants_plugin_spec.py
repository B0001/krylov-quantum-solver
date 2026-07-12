"""Spec gate for portfolio #20 — the pytest-invariants plugin (tasks 1-3, 5).

Red/green demo of the two decorators the solver actually needs (``lower_bound`` == the
variational-floor check, ``contains`` == bracket-brackets-reference), plus proof that the
failure messages are self-describing and that the plugin loads and emits its coverage report.

Deliberately self-contained: pure numbers, no pyscf/qiskit, so it is fast and safe to run in
the isolated ``test_*_spec.py`` process. The plugin/report checks shell out to a fresh pytest
so this module does not need the ``pytester`` fixture (which would require a root conftest).
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from pytest_invariants import InvariantViolation, invariant

# A stand-in for an FCI reference: H2/STO-3G ground state (Ha). Used as a floor and a bracket
# target, exactly as certchem.floor_guard and the golden suite use their FCI values.
FCI_H2 = -1.137270


# --- Task 3: lower_bound, green then red -------------------------------------------------


def test_lower_bound_passes_when_above_floor():
    @invariant.lower_bound(FCI_H2, provenance="FCI/STO-3G")
    def probe():
        return -1.100000  # a valid variational estimate: above the floor

    assert probe() is None  # returns None to pytest; no violation raised


def test_lower_bound_fails_below_floor_with_readable_message():
    # The "-hundreds of Hartree" pathology REFACTOR_PLAN.md documents: an estimate far below
    # the variational floor. The whole plugin exists to make this loud.
    @invariant.lower_bound(FCI_H2, provenance="FCI/STO-3G")
    def broken():
        return -813.0

    with pytest.raises(InvariantViolation) as ei:
        broken()
    msg = str(ei.value)
    # Message is legible without opening the plugin: relation, both values, provenance, margin.
    assert "Invariant `lower_bound` violated in broken" in msg
    assert "expected:  observed >= bound" in msg
    assert "observed:  -813" in msg
    assert "-1.13727" in msg
    assert "(provenance: FCI/STO-3G)" in msg
    assert "margin:" in msg


def test_lower_bound_key_extracts_from_object():
    class Bracket:
        best_estimate = -1.10

    @invariant.lower_bound(FCI_H2, key=lambda b: b.best_estimate, provenance="FCI")
    def probe():
        return Bracket()

    assert probe() is None


def test_lower_bound_accepts_callable_bound():
    @invariant.lower_bound(lambda: FCI_H2, provenance="FCI lookup")
    def probe():
        return -1.0

    assert probe() is None


# --- Task 3: contains, green then red ----------------------------------------------------


def test_contains_passes_when_reference_inside_interval():
    @invariant.contains(FCI_H2, provenance="FCI/STO-3G")
    def probe():
        return (-1.20, -1.05)  # (lower, upper) bracketing the reference

    assert probe() is None


def test_contains_fails_when_reference_outside_with_readable_message():
    @invariant.contains(FCI_H2, provenance="FCI/STO-3G")
    def bad_bracket():
        return (-2.00, -1.90)  # reference is above this interval

    with pytest.raises(InvariantViolation) as ei:
        bad_bracket()
    msg = str(ei.value)
    assert "Invariant `contains` violated in bad_bracket" in msg
    assert "expected:  lower <= reference <= upper" in msg
    assert "reference: -1.13727" in msg
    assert "(provenance: FCI/STO-3G)" in msg
    assert "interval:  [-2 Ha, -1.9 Ha]" in msg
    assert "ABOVE" in msg


def test_contains_reports_below_when_reference_under_interval():
    @invariant.contains(0.0, provenance="analytic")
    def bad_bracket():
        return (1.0, 2.0)

    with pytest.raises(InvariantViolation) as ei:
        bad_bracket()
    assert "BELOW" in str(ei.value)


def test_contains_key_maps_object_to_interval():
    class Bracket:
        lower, upper = -1.20, -1.05

    @invariant.contains(FCI_H2, key=lambda b: (b.lower, b.upper), provenance="FCI")
    def probe():
        return Bracket()

    assert probe() is None


# --- Task 1: stacking — one test, both invariants over the same observation ---------------


def test_stacked_invariants_both_run_and_report_metadata():
    @invariant.lower_bound(FCI_H2, key=lambda b: b["best"], provenance="FCI/STO-3G")
    @invariant.contains(FCI_H2, key=lambda b: (b["lo"], b["hi"]), provenance="FCI/STO-3G")
    def probe():
        return {"best": -1.10, "lo": -1.20, "hi": -1.05}

    assert probe() is None
    # Both specs recorded on the single composed wrapper (what the plugin reports over).
    names = {s.name for s in probe.__invariants__}
    assert names == {"lower_bound", "contains"}


def test_stacked_invariant_failure_names_the_broken_one():
    @invariant.lower_bound(FCI_H2, key=lambda b: b["best"], provenance="FCI")
    @invariant.contains(FCI_H2, key=lambda b: (b["lo"], b["hi"]), provenance="FCI")
    def probe():
        return {"best": -999.0, "lo": -1.20, "hi": -1.05}  # floor broken, bracket fine

    with pytest.raises(InvariantViolation) as ei:
        probe()
    assert "lower_bound" in str(ei.value)


# --- Task 2 & 5: the plugin loads and emits its coverage report ---------------------------

_INLINE_TEST = '''
from pytest_invariants import invariant

FCI = -1.137270

@invariant.contains(FCI, provenance="FCI/STO-3G")
def test_bracket_ok():
    return (-1.20, -1.05)

@invariant.lower_bound(FCI, provenance="FCI/STO-3G")
def test_floor_ok():
    return -1.10
'''


def _run_pytest(tmp_path: Path, *args: str) -> subprocess.CompletedProcess:
    test_file = tmp_path / "test_inline_invariants.py"
    test_file.write_text(textwrap.dedent(_INLINE_TEST))
    repo_root = Path(__file__).resolve().parent.parent
    return subprocess.run(
        [sys.executable, "-m", "pytest", str(test_file), "-p", "pytest_invariants", *args],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )


def test_plugin_loads_and_marker_registers(tmp_path):
    # `pytest -p pytest_invariants` loads clean and the auto-applied marker selects the tests.
    proc = _run_pytest(tmp_path, "-m", "invariants", "-v")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "2 passed" in proc.stdout
    # No "unknown marker" warning => the marker is registered by the plugin.
    assert "PytestUnknownMarkWarning" not in (proc.stdout + proc.stderr)


def test_invariant_report_emits_pass_counts(tmp_path):
    proc = _run_pytest(tmp_path, "--invariant-report")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    out = proc.stdout
    assert "invariant coverage report" in out
    # Both invariant kinds ran, each once, each passing.
    assert "lower_bound" in out
    assert "contains" in out
    assert "TOTAL" in out


# --- module-level decorated tests (so running THIS file under the plugin has real coverage) ---


@invariant.lower_bound(FCI_H2, provenance="FCI/STO-3G reference energy")
def test_estimate_respects_variational_floor():
    """A well-behaved estimator stays above the certified floor (cf. certchem.floor_guard)."""
    return -1.10  # stand-in Ritz estimate, above FCI_H2


@invariant.contains(FCI_H2, provenance="FCI/STO-3G reference energy")
def test_bracket_contains_fci():
    """The certified bracket brackets the FCI reference (cf. the golden-suite contains check)."""
    return (-1.20, -1.05)
