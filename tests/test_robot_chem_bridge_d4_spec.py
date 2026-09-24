"""
Acceptance gates G1-G4 for specs/SPEC_robot_chem_bridge_d4.md.

Claim: a real, independently-checked certkit_bridge.Verdict -- the same object
certkit_bridge.py emits for an actual molecule -- maps deterministically onto D1's
TwoSidedClaim contract and drives the existing, unmodified D1/D2/D3 chain
(robot.planner -> robot.workspace -> robot.failclosed) to a real outcome. This closes
CADE epic chem-cjl's stated acceptance criterion at the actuation leg.

Needs the pinned checker (`uv pip install -e ".[certkit]"`); runs a real H2 computation
through certkit_bridge.run_case once (module-level cache). make gates runs this in its
own process, same isolation as test_certkit_regression_gate_spec.py.
"""
import importlib.util
import os
from pathlib import Path

import pytest

if importlib.util.find_spec("certkit") is None:
    if os.environ.get("CI"):
        raise RuntimeError(
            "certkit is not installed: the chem-bridge gate cannot run. "
            'Install the pinned checker with `uv pip install -e ".[certkit]"`.'
        )
    pytest.skip('certkit extra not installed (uv pip install -e ".[certkit]")',
                allow_module_level=True)

from certkit.interval import Iv

import certkit_bridge
from robot import chem_bridge
from robot.failclosed import MonitorAction, monitor_trajectory
from robot.planner import PlanRefusal, Trajectory
from robot.workspace import (
    Verdict,
    WorkspaceEnvelope,
    classify_trajectory,
)

ORIGIN = (0.0, 0.0, 0.0)
TARGET = (0.5, 0.5, 0.5)

# Big enough that a tight certificate's trajectory is certifiably INSIDE.
ENVELOPE = WorkspaceEnvelope(x=Iv(-1.0, 1.0), y=Iv(-1.0, 1.0), z=Iv(-1.0, 1.0))

_CACHE = {}


def _h2_verdicts():
    """The real, independently-checked certificates certkit_bridge emits for H2."""
    if "H2" not in _CACHE:
        out = Path(__file__).resolve().parent.parent / "certkit_out" / "H2"
        _, verdicts = certkit_bridge.run_case("H2", out)
        _CACHE["H2"] = {v.name: v for v in verdicts}
    return _CACHE["H2"]


def _by_rule_and_status(verdicts, name: str, expect_ok: bool):
    v = verdicts[name]
    assert v.ok is expect_ok, (
        f"fixture drifted: {name} expected ok={expect_ok}, got ok={v.ok} ({v.line})"
    )
    return v


def test_G1_a_verified_certificate_plans_and_continues():
    verdicts = _h2_verdicts()
    v = _by_rule_and_status(verdicts, "certificate_temple", expect_ok=True)
    assert v.hi - v.lo < 1e-6, f"fixture drifted: certificate_temple width {v.hi - v.lo!r}"

    claim = chem_bridge.verdict_to_claim(v)
    assert claim.abstained is False
    assert claim.lower == 0.0
    assert claim.upper == pytest.approx(v.hi - v.lo)

    result = chem_bridge.plan_from_verdict(v, ORIGIN, TARGET)
    assert isinstance(result, Trajectory)

    assert classify_trajectory(result, ENVELOPE) == Verdict.INSIDE
    assert monitor_trajectory(result, ENVELOPE).action == MonitorAction.CONTINUE


def test_G2_an_abstained_certificate_refuses_with_the_checkers_own_reason():
    verdicts = _h2_verdicts()
    v = _by_rule_and_status(verdicts, "certificate_sector", expect_ok=False)

    claim = chem_bridge.verdict_to_claim(v)
    assert claim.abstained is True
    assert claim.reason == v.line
    assert claim.lower is None
    assert claim.upper is None

    result = chem_bridge.plan_from_verdict(v, ORIGIN, TARGET)
    assert isinstance(result, PlanRefusal)
    assert result.reason == v.line


def test_G3_a_verified_but_loose_certificate_still_trips_the_failclosed_floor():
    verdicts = _h2_verdicts()
    v = _by_rule_and_status(verdicts, "certificate_gershgorin", expect_ok=True)
    width = v.hi - v.lo
    assert width > 0.15, (
        f"fixture drifted: certificate_gershgorin width {width!r} no longer exceeds "
        "SAFETY_FLOOR -- this gate needs a verified-but-loose certificate"
    )

    claim = chem_bridge.verdict_to_claim(v)
    assert claim.abstained is False

    # D1's planner does not know about D3's floor: a verified certificate always plans.
    result = chem_bridge.plan_from_verdict(v, ORIGIN, TARGET)
    assert isinstance(result, Trajectory)

    # D3's independent floor is what catches it.
    verdict = monitor_trajectory(result, ENVELOPE)
    assert verdict.action == MonitorAction.HALT, (
        "a VERIFIED-but-loose certificate must still trip the fail-closed floor -- "
        "verified is not the same thing as safe-to-actuate"
    )


def test_G4_the_bridge_trusts_ok_never_lo_hi_presence():
    """Regression gate for a real footgun found while probing.

    certkit_bridge.Verdict populates lo/hi from the numbers passed into emit() BEFORE the
    checker runs, so an ABSTAIN verdict can still carry a finite, deceptively tight bracket
    (H2's real certificate_sector Verdict does exactly this). A bridge that gated on
    "lo/hi are not None and finite" instead of `ok` would silently forward an unverified
    number into an actuation plan.
    """
    deceptive = certkit_bridge.Verdict(
        name="synthetic",
        rule="temple_inertia",
        ok=False,
        line="ABSTAIN   synthetic non-discharged gap",
        lo=-1.0,
        hi=-1.0,
    )

    claim = chem_bridge.verdict_to_claim(deceptive)
    assert claim.abstained is True
    assert claim.reason == deceptive.line

    result = chem_bridge.plan_from_verdict(deceptive, ORIGIN, TARGET)
    assert isinstance(result, PlanRefusal)
