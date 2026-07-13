#!/usr/bin/env python3
"""
senseforge.pipeline -- config file -> sweep -> sensitivities -> candidates, one call.

The PRD's success criterion (sec 10): "Full pipeline reproducible from one config file + one
command." ``run_pipeline`` is that command's Python entry point; ``run_senseforge_sweep.py`` at
the repo root is the CLI wrapper.
"""
from __future__ import annotations

from pathlib import Path

from certchem.contract import Bracket, Certificate, CertifiedResult
from senseforge.candidates import build_candidates, write_candidate_report
from senseforge.config import SweepConfig, load_config
from senseforge.sensitivity import certified_central_differences
from senseforge.sweep import read_sweep_csv, run_sweep


def _row_to_result(row: dict) -> CertifiedResult:
    """A sweep CSV row, reconstituted as the CertifiedResult certified_central_differences
    expects -- the certificate fields are informational only (convergence is the one this
    pipeline reads back; the rest are round-tripped as placeholders since the CSV doesn't carry
    the full certificate)."""
    return CertifiedResult(
        bracket=Bracket(lower_hartree=float(row["lower"]), upper_hartree=float(row["upper"]),
                        best_estimate_hartree=float(row["gap"])),
        certificate=Certificate(method="from sweep CSV", floor_check="n/a",
                                krylov_dim=0, convergence=row["convergence"],
                                solver_version="n/a"),
    )


def run_pipeline(config: SweepConfig) -> Path:
    """Run the full SenseForge M1-M2 pipeline for one config. Returns the candidates.md path."""
    run_sweep(config)
    rows = read_sweep_csv(config)

    xs = [float(r["x"]) for r in rows]
    gaps = [float(r["gap"]) for r in rows]
    widths = [float(r["bracket_width"]) for r in rows]
    results = [_row_to_result(r) for r in rows]
    sensitivities = certified_central_differences(xs, results)

    candidates = build_candidates(config.halide, config.axis, xs, gaps, widths, sensitivities)
    return write_candidate_report(candidates, config.resolved_header(), config.cluster,
                                  config.output_dir)


def run_pipeline_from_file(config_path: str) -> Path:
    """CLI-facing entry point: load + validate the YAML config, then run the pipeline."""
    config = load_config(config_path)
    return run_pipeline(config)
