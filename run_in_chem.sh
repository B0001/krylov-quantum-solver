#!/usr/bin/env bash
#
# Validate the Krylov Quantum Solver in your local `chem` conda environment.
#   Usage:  bash run_in_chem.sh
#
# Everything runs via `conda run -n chem`, so you do not need to activate the env first.
set -euo pipefail
cd "$(dirname "$0")"
ENV="${CONDA_ENV:-chem}"

echo "============================================================"
echo "[1/4] Package versions in conda env '$ENV'"
echo "============================================================"
conda run -n "$ENV" python - <<'PY'
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
conda run -n "$ENV" python -m pip install -r requirements.txt

echo
echo "============================================================"
echo "[3/4] Test suite  (expect: 32 passed, then 2 passed/1 skipped)"
echo "============================================================"
# block2 (DMRG) initialises its own OpenMP runtime and segfaults if it loads into a process that
# already imported pyscf/qiskit-aer. Run the block2/DMRG tests in their OWN process to keep the
# suite green.
conda run -n "$ENV" python -m pytest tests/ \
    --ignore=tests/test_dmrg_reference.py --ignore=tests/test_hchain_extrapolation.py -q
conda run -n "$ENV" python -m pytest \
    tests/test_dmrg_reference.py tests/test_hchain_extrapolation.py -q

echo
echo "============================================================"
echo "[4/4] Benchmarks"
echo "============================================================"
conda run -n "$ENV" python benchmark_krylov.py
echo
conda run -n "$ENV" python benchmark_resources.py

echo
echo "All steps completed in conda env '$ENV'."
echo "Optional (slower): conda run -n $ENV python benchmark_n2.py"
