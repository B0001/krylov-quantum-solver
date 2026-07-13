#!/usr/bin/env python3
"""
run_senseforge_sweep.py -- SenseForge M1-M2 CLI: one YAML config -> one candidates.md.

    python run_senseforge_sweep.py senseforge/configs/nb3cl8_strain.yaml

See senseforge/hamiltonian.py for the two documented deviations from the literal task spec
(no ab-initio CIF-strain path; gaps are exact, not certified_gaps.py Krylov brackets) and
specs/tasks/04-senseforge.md / specs/full/spec-nb3x8-sensor-designer.md for the intended scope.
"""
from __future__ import annotations

import sys
import time

from senseforge.config import ConfigError
from senseforge.pipeline import run_pipeline_from_file


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <config.yaml>", file=sys.stderr)
        return 2
    config_path = sys.argv[1]

    t0 = time.time()
    try:
        report_path = run_pipeline_from_file(config_path)
    except ConfigError as exc:
        print(f"config error in field {exc.field!r}: {exc}", file=sys.stderr)
        return 1
    elapsed = time.time() - t0

    from senseforge.config import load_config
    n_points = len(load_config(config_path).grid())
    print(f"SenseForge sweep complete: {report_path}")
    print(f"  {n_points} grid points in {elapsed:.1f}s ({elapsed / n_points:.3f}s/point)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
