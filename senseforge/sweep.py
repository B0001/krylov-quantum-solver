#!/usr/bin/env python3
"""
senseforge.sweep -- the grid -> certified-gap -> CSV driver (task 3).

Resumability (ADR-0008's cache-key idea, applied without a persistent service-layer cache):
ADR-0008 (architecture/adr/0008-content-hash-cache-and-async-jobs.md) describes a content-hash
cache keyed on canonicalized (geometry, basis, active space, ...) -- no such persistent cache is
implemented anywhere in this repo yet (checked: only ephemeral in-process dicts exist elsewhere,
e.g. sweep_hybrid_solver.py). "Resumable... for free" (task 3) does not hold literally. This
module gets the SAME observable behavior (kill mid-sweep, resume, identical final CSV) the direct
way instead: each grid point is a row keyed by its (rounded) x-value, written to the CSV
IMMEDIATELY after being computed (not batched), and a re-run reads the existing rows first and
skips any x already present. That is sufficient for the stated gate; a real persistent cache
(keyed on ``SweepConfig.content_hash()``, shared across runs/output dirs) is a separate,
un-built follow-up.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List

from certchem.contract import CertifiedResult
from senseforge.config import SweepConfig
from senseforge.hamiltonian import certified_strain_gap, zeeman_split_gap
from senseforge.headers import csv_header_lines

_FIELDNAMES = ["x", "gap", "lower", "upper", "bracket_width", "convergence"]


def compute_gap(config: SweepConfig, x: float) -> CertifiedResult:
    """The certified gap at grid value ``x`` for ``config.axis`` (see hamiltonian.py)."""
    if config.axis == "strain":
        return certified_strain_gap(config.system, x)
    if config.axis == "field":
        return zeeman_split_gap(config.system, x)
    raise ValueError(f"unknown axis {config.axis!r}")  # config.py already validates this


def csv_path(config: SweepConfig) -> Path:
    return Path(config.output_dir) / f"gap_vs_{config.axis}.csv"


def _read_existing(path: Path) -> Dict[str, dict]:
    if not path.exists():
        return {}
    with path.open(newline="") as f:
        data_lines = [ln for ln in f if not ln.startswith("#")]
    if not data_lines:
        return {}
    rows: Dict[str, dict] = {}
    for row in csv.DictReader(data_lines):
        rows[row["x"]] = row
    return rows


def run_sweep(config: SweepConfig, *, _limit: int = None) -> Path:
    """Run (or resume) the sweep, writing one row per grid point to ``gap_vs_{axis}.csv``.

    ``_limit``: internal test hook -- stop after computing this many NEW points (simulates a
    kill mid-sweep). Not part of the public interface.
    """
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = csv_path(config)

    existing = _read_existing(path)
    xs = config.grid()
    is_new_file = not path.exists()

    mode = "w" if is_new_file else "a"
    computed = 0
    with path.open(mode, newline="") as f:
        if is_new_file:
            for line in csv_header_lines(config.resolved_header(), config.cluster):
                f.write(line + "\n")
            writer = csv.DictWriter(f, fieldnames=_FIELDNAMES)
            writer.writeheader()
        else:
            writer = csv.DictWriter(f, fieldnames=_FIELDNAMES)

        for x in xs:
            key = repr(x)
            if key in existing:
                continue
            if _limit is not None and computed >= _limit:
                break
            result = compute_gap(config, x)
            row = {
                "x": key,
                "gap": result.bracket.best_estimate_hartree,
                "lower": result.bracket.lower_hartree,
                "upper": result.bracket.upper_hartree,
                "bracket_width": result.bracket.width,
                "convergence": result.certificate.convergence,
            }
            writer.writerow(row)
            f.flush()
            computed += 1

    return path


def read_sweep_csv(config: SweepConfig) -> List[dict]:
    """The full sweep result, ordered by ``config.grid()`` -- raises if any point is missing."""
    rows = _read_existing(csv_path(config))
    xs = config.grid()
    missing = [x for x in xs if repr(x) not in rows]
    if missing:
        raise ValueError(f"sweep incomplete: missing grid points {missing[:5]}...")
    return [rows[repr(x)] for x in xs]
