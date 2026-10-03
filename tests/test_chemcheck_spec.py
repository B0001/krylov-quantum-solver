"""Spec gate for ChemCheck M1 + Mode A (specs/SPEC_chemcheck.md, tasks 1-5,7,8)."""

import functools
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


@functools.lru_cache(maxsize=None)
def _live(name):
    return recompute_tier_reference(TIERS[name])


@pytest.mark.parametrize("name", ["T0", "T1", "T2", "T3"])
def test_frozen_tier_values_match_live_recompute(name):
    """Platform-independent part of the freeze: FCI energy and Pauli-term count."""
    tier, live = TIERS[name], _live(name)
    assert live["hamiltonian_pauli_terms"] == tier.hamiltonian_pauli_terms
    assert abs(live["fci_reference_hartree"] - tier.fci_reference_hartree) < 1e-6


# The CX count is portable for T0-T2 but NOT for T3: build_trotter_step orders terms by |coefficient|
# (trotter_krylov.canonical_term_order), a sign gauge cannot change that order but an MO rotation can.
# Measured on macOS 27 (Apple M3, scipy 1.15.3/Accelerate, 2026-10-02): T0/T1/T2 match; T3 live 6530
# vs frozen 6528, identical for PYTHONHASHSEED 0, 1 and 2. This contradicts the hand-off's "CX count
# passes everywhere" (the original single test asserted the hash first, so it never reached the CX).
@pytest.mark.parametrize("name", [
    "T0", "T1", "T2",
    pytest.param("T3", marks=pytest.mark.skipif(
        sys.platform != "linux",
        reason="T3 CX count is MO-gauge dependent (Trotter order is by |coeff|): macOS 27 gives "
               "6530 vs frozen 6528, deterministic over PYTHONHASHSEED 0/1/2. "
               "specs/SPEC_chemcheck.md platform note.")),
])
def test_frozen_tier_cx_count_matches_live_recompute(name):
    assert _live(name)["two_qubit_gates_per_trotter_step"] == TIERS[name].two_qubit_gates_per_trotter_step


# The canonical hash is a bitwise hash of SCF-derived integrals in whatever MO gauge the platform's
# eigensolver returns, so the frozen values only reproduce on the freezing platform (Linux). Measured
# on macOS 27 (same machine): T0 and T2 match; T1 live fff8a314e967 != frozen 2fe671eae4ae,
# reproduced exactly by flipping the sign of one MO (4 of the 16 sign patterns match); T3 live
# ccde9d59df68 != frozen 382de579ca32 and none of the 64 sign patterns matches -- its active space
# holds two degenerate pi pairs (eps = -0.57139 x2, 0.28019 x2), and its CX count differs too, so an
# in-pair rotation is inferred (not verified: needs the Linux orbitals). FCI energy and Pauli-term
# count above still match on every platform. specs/SPEC_chemcheck.md platform note.
@pytest.mark.skipif(
    sys.platform != "linux",
    reason="hash is MO-gauge dependent, frozen on Linux; macOS 27: T1 differs by an MO sign flip, "
           "T3 by more (no sign pattern of 64 matches; degenerate pi pairs), T0/T2 match. "
           "specs/SPEC_chemcheck.md platform note.",
)
@pytest.mark.parametrize("name", ["T0", "T1", "T2", "T3"])
def test_frozen_tier_hash_matches_live_recompute(name):
    assert _live(name)["hamiltonian_sha256"] == TIERS[name].hamiltonian_sha256


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
