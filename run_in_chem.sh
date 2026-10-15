#!/usr/bin/env bash
#
# Validate the Krylov Quantum Solver with `uv`, the project's exclusive dependency manager.
#   Usage:  bash run_in_chem.sh
#
# Everything runs via `uv run`, which resolves/uses the project's `.venv` automatically --
# no env to create or activate first.
set -euo pipefail
cd "$(dirname "$0")"

echo "============================================================"
echo "[1/4] Package versions"
echo "============================================================"
uv run python - <<'PY'
import importlib
for pkg in ["qiskit", "qiskit_nature", "qiskit_aer", "pyscf", "ase", "scipy", "numpy"]:
    try:
        m = importlib.import_module(pkg)
        print(f"  {pkg:14s} {getattr(m, '__version__', '?')}")
    except Exception as e:
        print(f"  {pkg:14s} MISSING ({e.__class__.__name__})")
PY

echo
echo "============================================================"
echo "[2/4] Install / refresh dependencies"
echo "============================================================"
uv sync --extra dmrg --extra test --extra dev

echo
echo "============================================================"
echo "[3/4] Test suite  (expect: 32 passed, then 2 passed/1 skipped)"
echo "============================================================"
# block2 (DMRG) initialises its own OpenMP runtime and segfaults if it loads into a process that
# already imported pyscf/qiskit-aer. Run the block2/DMRG tests (test_dmrg_reference + every spec
# gate, test_*_spec.py) in their OWN process. The glob keeps this correct as new gates are added.
uv run pytest tests/ \
    --ignore-glob='tests/test_*_spec.py' --ignore=tests/test_dmrg_reference.py -q
uv run pytest tests/test_dmrg_reference.py tests/test_*_spec.py -q

echo
echo "============================================================"
echo "[4/4] Benchmarks"
echo "============================================================"
uv run python benchmark_krylov.py
echo
uv run python benchmark_resources.py

echo
echo "All steps completed via uv."
echo "Optional (slower): uv run python benchmark_n2.py"
