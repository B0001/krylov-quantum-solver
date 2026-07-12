"""Spec gate for ChemCheck M1 + Mode A (specs/SPEC_chemcheck.md, tasks 1-5,7,8)."""

import subprocess
import sys

import pytest
from jsonschema import Draft202012Validator

from chemcheck import (
    BENCHMARK_VERSION,
    ROUTING_OVERHEAD,
    TIERS,
    expected_total_error,
    headroom_factor,
    mode_b_energy_verdict,
    recompute_tier_reference,
    render_markdown,
    required_two_qubit_error,
    routing_overhead,
    score_mode_a,
)
from chemcheck.scorecard import scorecard_schema
from chemcheck.submission import SubmissionError, validate_submission

_SCOREABLE = [t for t in TIERS.values() if not t.aspirational]


# --- Gate 1: registry loads + frozen values verified ------------------------------------


def test_registry_loads_without_solver_deps():
    code = (
        "import chemcheck.tiers, chemcheck.budget, sys;"
        "bad=[m for m in ('pyscf','qiskit','numpy','scipy') if m in sys.modules];"
        "assert not bad, bad; print('clean')"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert "clean" in proc.stdout


def test_benchmark_version_and_t4_aspirational():
    assert BENCHMARK_VERSION == "chemcheck-2026.1"
    assert TIERS["T4"].aspirational is True
    assert TIERS["T4"].fci_reference_hartree is None
    assert TIERS["T4"].hamiltonian_sha256 is None


@pytest.mark.parametrize("name", ["T0", "T1", "T2", "T3"])
def test_frozen_tier_values_match_live_recompute(name):
    tier = TIERS[name]
    live = recompute_tier_reference(tier)
    assert live["hamiltonian_sha256"] == tier.hamiltonian_sha256
    assert live["two_qubit_gates_per_trotter_step"] == tier.two_qubit_gates_per_trotter_step
    assert live["hamiltonian_pauli_terms"] == tier.hamiltonian_pauli_terms
    assert abs(live["fci_reference_hartree"] - tier.fci_reference_hartree) < 1e-6


# --- Gate 2: submission validation ------------------------------------------------------

_VALID_DEVICE = {
    "benchmark_version": "chemcheck-2026.1",
    "device_spec": {
        "name": "TestQPU", "qubit_count": 27, "two_qubit_error": 5e-3,
        "connectivity": "heavy_hex", "native_gates": ["cx", "rz", "sx"],
    },
}


def test_valid_submissions_pass():
    validate_submission(_VALID_DEVICE)  # no runs -> Mode A, valid
    validate_submission({**_VALID_DEVICE, "device_spec": {
        **_VALID_DEVICE["device_spec"], "connectivity": "all_to_all", "t1_us": 100.0}})
    validate_submission({**_VALID_DEVICE, "device_spec": {
        **_VALID_DEVICE["device_spec"], "connectivity": "linear", "one_qubit_error": 1e-4}})


@pytest.mark.parametrize("mutate,bad_field", [
    (lambda s: s.pop("device_spec"), "device_spec"),
    (lambda s: s["device_spec"].pop("two_qubit_error"), "two_qubit_error"),
    (lambda s: s["device_spec"].update(two_qubit_error=2.0), "two_qubit_error"),
    (lambda s: s["device_spec"].update(connectivity="quantum_teleporter"), "connectivity"),
    (lambda s: s["device_spec"].update(qubit_count=0), "qubit_count"),
    (lambda s: s.update(benchmark_version="v1"), "benchmark_version"),
])
def test_invalid_submissions_rejected_with_pointer(mutate, bad_field):
    import copy
    sub = copy.deepcopy(_VALID_DEVICE)
    mutate(sub)
    with pytest.raises(SubmissionError) as ei:
        validate_submission(sub)
    assert bad_field in ei.value.path or bad_field in str(ei.value)


# --- Gate 3: routing overhead -----------------------------------------------------------


def test_routing_overhead_ordering():
    assert routing_overhead("all_to_all") == 1.0
    assert (ROUTING_OVERHEAD["heavy_hex"] > ROUTING_OVERHEAD["grid"]
            > ROUTING_OVERHEAD["all_to_all"])
    with pytest.raises(ValueError):
        routing_overhead("nonsense")


# --- Gate 4: error budget is a pure function --------------------------------------------


def test_expected_total_error_hand_cases():
    assert expected_total_error(10, 1.0, 0.0) == 0.0
    assert expected_total_error(1, 1.0, 1.0) == 1.0
    assert expected_total_error(1, 1.0, 0.1) == pytest.approx(0.1)
    assert expected_total_error(2, 1.0, 0.1) == pytest.approx(1 - 0.9**2)
    assert expected_total_error(3, 2.0, 0.01) == pytest.approx(1 - 0.99**6)


# --- Gate 5: headroom monotonicity ------------------------------------------------------


def test_headroom_threshold_and_monotonicity():
    req = required_two_qubit_error(1000, 1.0)  # p_required = 1/1000 = 1e-3
    assert req == pytest.approx(1e-3)
    assert headroom_factor(req, req) == pytest.approx(1.0)  # exactly at threshold
    assert headroom_factor(5e-4, req) < 1.0  # better error -> PASS
    assert headroom_factor(2e-3, req) > 1.0  # worse error -> FAIL
    # monotone decreasing in device quality
    hs = [headroom_factor(p, req) for p in (1e-2, 1e-3, 1e-4)]
    assert hs[0] > hs[1] > hs[2]


# --- Gate 6: floor detector -------------------------------------------------------------


@pytest.mark.parametrize("tier", _SCOREABLE)
def test_floor_detector_flags_known_bad(tier):
    # Old-codebase-style garbage: hundreds of Ha below the true ground state.
    v = mode_b_energy_verdict(tier.fci_reference_hartree - 500.0, tier)
    assert v["result"] == "UNPHYSICAL"
    assert v["floor_check"] == "violation"


@pytest.mark.parametrize("tier", _SCOREABLE)
def test_floor_detector_no_false_positive_on_golden(tier):
    v = mode_b_energy_verdict(tier.fci_reference_hartree, tier)  # exact reference
    assert v["result"] == "PASS"
    assert v["floor_check"] == "pass"


def test_energy_accuracy_bands():
    t = TIERS["T0"]
    f = t.fci_reference_hartree
    assert mode_b_energy_verdict(f + 1.0e-3, t)["result"] == "PASS"       # 1.0 mHa
    assert mode_b_energy_verdict(f + 10.0e-3, t)["result"] == "MARGINAL"  # 10 mHa
    assert mode_b_energy_verdict(f + 50.0e-3, t)["result"] == "FAIL"      # 50 mHa


# --- Gate 7: scorecard emitter ----------------------------------------------------------


def test_scorecard_validates_against_schema():
    card = score_mode_a(_VALID_DEVICE)
    Draft202012Validator(scorecard_schema()).validate(card)  # raises on invalid
    assert card["device_name"] == "TestQPU"
    assert {t["tier"] for t in card["tiers"]} == {"T0", "T1", "T2", "T3"}


def test_good_device_passes_t0_bad_device_fails_t3():
    good = score_mode_a({**_VALID_DEVICE, "device_spec": {
        **_VALID_DEVICE["device_spec"], "connectivity": "all_to_all", "two_qubit_error": 1e-5}})
    verdict = {t["tier"]: t["mode_a"]["result"] for t in good["tiers"]}
    assert verdict["T0"] == "PASS"

    bad = score_mode_a({**_VALID_DEVICE, "device_spec": {
        **_VALID_DEVICE["device_spec"], "connectivity": "linear", "two_qubit_error": 5e-2}})
    bad_verdict = {t["tier"]: t["mode_a"]["result"] for t in bad["tiers"]}
    assert bad_verdict["T3"] == "FAIL"


def test_render_markdown_carries_disclaimer():
    md = render_markdown(score_mode_a(_VALID_DEVICE))
    assert "not quantum advantage" in md
    assert "| T0 |" in md and "| T3 |" in md
