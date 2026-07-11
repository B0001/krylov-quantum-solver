"""Mode A scorer + scorecard emitter (task 8) and the Mode B floor detector (task 7).

``score_mode_a`` maps a device spec sheet to a schema-valid scorecard; ``render_markdown``
renders it with the ``classically_simulable`` disclaimer. ``mode_b_energy_verdict`` is the
anti-fraud core — it reuses ``certchem.floor_guard`` so a sub-floor ("UNPHYSICAL") energy can
never be scored as a near miss.
"""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from certchem import FloorViolationError, floor_guard

from .budget import headroom_factor, required_two_qubit_error, routing_overhead
from .tiers import (
    ACCURACY_MARGINAL_MHA,
    ACCURACY_PASS_MHA,
    BENCHMARK_VERSION,
    TIERS,
    Tier,
)

_SCHEMA_PATH = (
    Path(__file__).resolve().parents[1]
    / "architecture" / "interfaces" / "chemcheck-scorecard.schema.json"
)
#: Below-floor energies within this many Ha of FCI are shot noise, not fraud (avoids false
#: positives on golden results, which sit AT the reference).
_FLOOR_TOLERANCE_HA = 1e-3


@lru_cache(maxsize=1)
def scoring_code_hash() -> str:
    """Version hash pinning the scoring code (reproducibility bundle). SHA-256 of the modules."""
    src = b""
    for name in ("tiers.py", "budget.py", "scorecard.py", "submission.py"):
        src += (Path(__file__).with_name(name)).read_bytes()
    return f"sha256:{hashlib.sha256(src).hexdigest()}"


def _score_tier_mode_a(tier: Tier, device: dict[str, Any]) -> dict[str, Any]:
    overhead = routing_overhead(device["connectivity"])
    required = required_two_qubit_error(tier.two_qubit_gates_per_trotter_step, overhead)
    headroom = headroom_factor(device["two_qubit_error"], required)
    return {
        "tier": tier.name,
        "classically_simulable": tier.classically_simulable,
        "mode_a": {
            "result": "PASS" if headroom <= 1.0 else "FAIL",
            "headroom_factor": headroom,
            "required_two_qubit_error": required,
            "model_uncertainty": "v1 depolarizing single-step model, uncalibrated (see SPEC)",
        },
    }


def score_mode_a(submission: dict[str, Any]) -> dict[str, Any]:
    """Score a submission's device spec against every scoreable tier (T0–T3).

    Returns a dict validating against chemcheck-scorecard.schema.json. Aspirational tiers (T4)
    are skipped. Does not itself validate the submission — call ``validate_submission`` first.
    """
    device = submission["device_spec"]
    tiers = [_score_tier_mode_a(t, device) for t in TIERS.values() if not t.aspirational]
    n_pass = sum(t["mode_a"]["result"] == "PASS" for t in tiers)
    highest = next(
        (t["tier"] for t in reversed(tiers) if t["mode_a"]["result"] == "PASS"), "none"
    )
    verdict = (
        f"{device['name']}: passes {n_pass}/{len(tiers)} tiers (Mode A, paper score); "
        f"highest passing tier = {highest}. Passing T0–T2 certifies stack correctness, "
        "not quantum advantage."
    )
    return {
        "benchmark_version": submission.get("benchmark_version", BENCHMARK_VERSION),
        "device_name": device["name"],
        "scoring_code_hash": scoring_code_hash(),
        "tiers": tiers,
        "verdict": verdict,
    }


def mode_b_energy_verdict(energy_hartree: float, tier: Tier) -> dict[str, Any]:
    """Score one measured energy against a tier (anti-fraud floor + accuracy bands).

    Floor check first: an energy below the variational floor (the exact FCI ground state) is
    ``UNPHYSICAL`` regardless of how close it looks — reusing ``certchem.floor_guard``.
    """
    floor = tier.fci_reference_hartree
    try:
        floor_guard(energy_hartree, floor, tol=_FLOOR_TOLERANCE_HA)
    except FloorViolationError:
        return {"result": "UNPHYSICAL", "floor_check": "violation",
                "energy_error_mha_raw": (energy_hartree - floor) * 1e3}

    err_mha = abs(energy_hartree - floor) * 1e3
    if err_mha <= ACCURACY_PASS_MHA:
        result = "PASS"
    elif err_mha <= ACCURACY_MARGINAL_MHA:
        result = "MARGINAL"
    else:
        result = "FAIL"
    return {"result": result, "floor_check": "pass", "energy_error_mha_raw": err_mha}


def render_markdown(scorecard: dict[str, Any]) -> str:
    """Render a scorecard as a GitHub-readable Markdown table with the honesty disclaimer."""
    lines = [
        f"# ChemCheck scorecard — {scorecard['device_name']}",
        "",
        f"Benchmark: `{scorecard['benchmark_version']}` · scoring code "
        f"`{scorecard['scoring_code_hash']}`",
        "",
        "| Tier | Classically simulable | Mode A | Headroom (×) | Req. 2q error |",
        "|------|----------------------|--------|--------------|---------------|",
    ]
    for t in scorecard["tiers"]:
        a = t["mode_a"]
        sim = "yes ⚠️" if t["classically_simulable"] else "no"
        lines.append(
            f"| {t['tier']} | {sim} | **{a['result']}** | {a['headroom_factor']:.2f} | "
            f"{a['required_two_qubit_error']:.2e} |"
        )
    lines += [
        "",
        f"**Verdict:** {scorecard['verdict']}",
        "",
        "> ⚠️ *Classically simulable (T0–T2): passing certifies stack correctness, not quantum "
        "advantage.* Headroom factors are order-of-magnitude (v1 uncalibrated model).",
    ]
    return "\n".join(lines)


def scorecard_schema() -> dict[str, Any]:
    """The output schema, for callers that want to validate a scorecard themselves."""
    return json.loads(_SCHEMA_PATH.read_text())
