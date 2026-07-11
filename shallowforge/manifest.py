"""Provenance manifest emitter (task 2), per architecture/interfaces/compiler-manifest.schema.json.

Every transform appends a :class:`TransformEntry` (name, parameters, predicted ε, lossless flag);
totals are assembled at export. A lossless transform emits ``predicted_epsilon_mha == 0``. The
manifest doubles as ChemCheck's compilation-audit artifact (ADR-0006/0007) and — per ADR-0007 —
a gate count is never reported without its ε (``cx_at_epsilon_claim``).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

import jsonschema
from jsonschema import Draft202012Validator

MANIFEST_VERSION = "shallowforge-1"

_SCHEMA_PATH = (
    Path(__file__).resolve().parents[1]
    / "architecture" / "interfaces" / "compiler-manifest.schema.json"
)

#: Transforms the schema permits (kept in sync with the schema enum).
_ALLOWED_TRANSFORMS = frozenset({
    "pauli_grouping", "term_ordering", "gate_cancellation", "trotter_order",
    "randomized_compiling", "thc_truncation", "qubitization", "peephole_optimization",
    "template_matching", "resynthesis",
})


@dataclass(frozen=True)
class TransformEntry:
    """One applied transform in the compilation stack."""

    transform: str
    parameters: dict[str, Any] = field(default_factory=dict)
    predicted_epsilon_mha: float = 0.0
    lossless: bool = True

    def __post_init__(self) -> None:
        if self.transform not in _ALLOWED_TRANSFORMS:
            raise ValueError(
                f"unknown transform {self.transform!r}; allowed: {sorted(_ALLOWED_TRANSFORMS)}"
            )
        if self.lossless and self.predicted_epsilon_mha != 0.0:
            raise ValueError("a lossless transform must emit predicted_epsilon_mha == 0")
        if self.predicted_epsilon_mha < 0:
            raise ValueError("predicted_epsilon_mha must be >= 0")

    def to_dict(self) -> dict[str, Any]:
        return {
            "transform": self.transform,
            "parameters": dict(self.parameters),
            "predicted_epsilon_mha": float(self.predicted_epsilon_mha),
            "lossless": bool(self.lossless),
        }


@lru_cache(maxsize=1)
def _validator() -> Draft202012Validator:
    return Draft202012Validator(json.loads(_SCHEMA_PATH.read_text()))


def validate_manifest(manifest: dict[str, Any]) -> None:
    """Raise ``jsonschema.ValidationError`` if the manifest violates the schema."""
    _validator().validate(manifest)


def build_manifest(
    hamiltonian_hash: str,
    stack: list[TransformEntry],
    *,
    cx_per_step: int,
    depth: int,
    ancillas: int,
    solver_version: str | None = None,
) -> dict[str, Any]:
    """Assemble a schema-valid manifest. Totals' predicted ε is the sum over the stack.

    ``cx_at_epsilon_claim`` is populated at the fixed ε = 1.6 mHa target (ADR-0007) so the gate
    count never travels without its ε. Validates before returning.
    """
    predicted_total = sum(e.predicted_epsilon_mha for e in stack)
    manifest: dict[str, Any] = {
        "manifest_version": MANIFEST_VERSION,
        "hamiltonian_hash": hamiltonian_hash,
        "stack": [e.to_dict() for e in stack],
        "totals": {
            "cx_per_step": int(cx_per_step),
            "depth": int(depth),
            "ancillas": int(ancillas),
            "predicted_epsilon_total_mha": float(predicted_total),
            "cx_at_epsilon_claim": {"epsilon_mha": 1.6, "cx_per_step": int(cx_per_step)},
        },
    }
    if solver_version is not None:
        manifest["solver_version"] = solver_version
    validate_manifest(manifest)
    return manifest


ValidationError = jsonschema.ValidationError
