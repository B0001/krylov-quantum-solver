#!/usr/bin/env python3
"""
senseforge.config -- one YAML file per sweep run.

Loader validates every field and names the offending one on failure (task 1, ``✓ Bad configs
fail with named field``); the resolved config is hashed (SHA-256 over canonical JSON) so every
output artifact can record exactly which config produced it (ADR-0008's cache-key idea, applied
to provenance rather than a persistent cache -- see senseforge/sweep.py module docstring for why
no persistent content-hash cache is implemented).
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import yaml

from nb3x8_gaps import NB3X8_LT_BULK

#: axis -> (grid_min, grid_max, grid_step) defaults, per the PRD (section 3).
_AXIS_DEFAULTS = {
    "strain": (-0.02, 0.02, 0.0025),
    "field": (0.0, 10.0, 0.5),
}


class ConfigError(ValueError):
    """A sweep config failed validation. ``.field`` names the offending key."""

    def __init__(self, message: str, *, field: str) -> None:
        super().__init__(message)
        self.field = field


@dataclass(frozen=True)
class SweepConfig:
    """A fully-resolved, validated SenseForge sweep run.

    ``cluster`` is fixed (not a sweep knob): this repo's only validated Nb3X8 model is the
    downfolded interlayer dimer (2-orbital generalized Hubbard cluster, 4 spin-orbitals) -- see
    senseforge/hamiltonian.py.

    NO ``krylov_dim`` FIELD (removed 2026-07-12). It used to sit here "for provenance/interface
    parity", claiming it "becomes the resolution knob for the Gate-1 convergence check". That was
    false on both counts: nothing read it (every ``Certificate`` here is built with
    ``krylov_dim=0``, because the gaps come from exact diagonalization and no Krylov subspace is
    ever constructed -- hamiltonian.py deviation (2)), and ``validation.py`` never referenced it.
    Worse, it was hashed into ``content_hash()`` and therefore STAMPED ON EVERY PUBLISHED DESIGN
    CARD as ``krylov_dim=12`` -- advertising a Krylov dimension for a calculation that used none.
    A dead field is untidy; a dead field that writes false provenance into a shipped artifact is
    a defect. See SPEC_senseforge.md.
    """

    halide: str                  # "Cl" | "Br" | "I"
    axis: str                    # "strain" | "field"
    grid_min: float
    grid_max: float
    grid_step: float
    output_dir: str = "results/senseforge"
    cluster: str = "Nb3X8 downfolded interlayer dimer (2-orbital generalized Hubbard cluster)"

    @property
    def system(self) -> str:
        return f"Nb3{self.halide}8"

    def grid(self) -> list:
        """The swept values, inclusive of ``grid_max`` (within float tolerance)."""
        n = round((self.grid_max - self.grid_min) / self.grid_step)
        return [round(self.grid_min + i * self.grid_step, 10) for i in range(n + 1)]

    def content_hash(self) -> str:
        """SHA-256 over the canonicalized (sorted-key) config, hex digest (ADR-0008 style key)."""
        payload = json.dumps(asdict(self), sort_keys=True)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def resolved_header(self) -> dict:
        """The dict every output artifact's header echoes (task 1: ``resolved config hash
        recorded``)."""
        return {**asdict(self), "system": self.system, "config_hash": self.content_hash()}


def _require(d: dict, key: str, cast, *, default=None) -> Any:
    if key not in d:
        if default is not None:
            return default
        raise ConfigError(f"missing required field {key!r}", field=key)
    try:
        return cast(d[key])
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"field {key!r}: {exc}", field=key) from exc


def load_config(path: str) -> SweepConfig:
    """Load and validate a YAML sweep config. Raises :class:`ConfigError` naming the bad field."""
    raw = yaml.safe_load(Path(path).read_text())
    if not isinstance(raw, dict):
        raise ConfigError("config root must be a mapping", field="<root>")

    halide = _require(raw, "halide", str)
    system = f"Nb3{halide}8"
    if system not in NB3X8_LT_BULK:
        raise ConfigError(
            f"halide {halide!r} -> {system!r} not in validated set "
            f"{sorted(h[3:-1] for h in NB3X8_LT_BULK)}",
            field="halide",
        )

    axis = _require(raw, "axis", str)
    if axis not in _AXIS_DEFAULTS:
        raise ConfigError(f"axis must be one of {sorted(_AXIS_DEFAULTS)}, got {axis!r}",
                          field="axis")
    default_min, default_max, default_step = _AXIS_DEFAULTS[axis]

    grid_min = _require(raw, "grid_min", float, default=default_min)
    grid_max = _require(raw, "grid_max", float, default=default_max)
    grid_step = _require(raw, "grid_step", float, default=default_step)
    if grid_step <= 0:
        raise ConfigError(f"grid_step must be positive, got {grid_step}", field="grid_step")
    if grid_max <= grid_min:
        raise ConfigError(f"grid_max ({grid_max}) must exceed grid_min ({grid_min})",
                          field="grid_max")

    if "krylov_dim" in raw:
        raise ConfigError(
            "krylov_dim was removed: SenseForge gaps come from exact diagonalization, no Krylov "
            "subspace is built, and the field only ever wrote false provenance into artifacts "
            "(see SweepConfig docstring / SPEC_senseforge.md)",
            field="krylov_dim",
        )

    output_dir = _require(raw, "output_dir", str, default=f"results/senseforge/{halide}")

    return SweepConfig(halide=halide, axis=axis, grid_min=grid_min, grid_max=grid_max,
                       grid_step=grid_step, output_dir=output_dir)
