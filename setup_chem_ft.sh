#!/usr/bin/env bash
#
# setup_chem_ft.sh
# Creates an isolated FT-resource-estimation env by CLONING your working `chem` env,
# so `chem` is never touched. Adds only openfermion (+ pytest, which resource_estimates
# imports). Verifies the result and smoke-tests the estimator.
#
# Usage:   bash setup_chem_ft.sh
#
set -euo pipefail

SRC_ENV="${SRC_ENV:-chem}"
NEW_ENV="${NEW_ENV:-chem-ft}"

say() { printf '\n\033[1m[setup]\033[0m %s\n' "$*"; }
die() { printf '\n\033[31m[setup] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

# 0. sanity
command -v conda >/dev/null 2>&1 || die "conda not found on PATH. Open a shell where 'conda' works."
conda env list | awk '{print $1}' | grep -qx "$SRC_ENV" \
  || die "source env '$SRC_ENV' not found. Set SRC_ENV=<name> and retry."

if conda env list | awk '{print $1}' | grep -qx "$NEW_ENV"; then
  die "env '$NEW_ENV' already exists. Remove it (conda env remove -n $NEW_ENV) or set NEW_ENV=<name>."
fi

# 1. clone the known-good env (inherits your working pyscf -- the painful one on Mac)
say "Cloning '$SRC_ENV' -> '$NEW_ENV' (this copies your working pyscf/qiskit/sqd stack)..."
conda create -y -n "$NEW_ENV" --clone "$SRC_ENV"

# 2. add ONLY the FT extras, into the clone (never into $SRC_ENV)
say "Installing openfermion + pytest into '$NEW_ENV'..."
conda run -n "$NEW_ENV" python -m pip install --quiet openfermion pytest

# 3. verify nothing broke
say "Running 'pip check' in '$NEW_ENV'..."
if conda run -n "$NEW_ENV" python -m pip check; then
  say "pip check: no broken requirements."
else
  die "pip check reported conflicts in '$NEW_ENV'. Your '$SRC_ENV' is still untouched; \
inspect the clone or remove it with: conda env remove -n $NEW_ENV"
fi

# 4. confirm the combined stack imports together
say "Import test (qiskit + sqd + pyscf + resource_estimates)..."
conda run -n "$NEW_ENV" python - <<'PY'
import qiskit, qiskit_addon_sqd, pyscf
from openfermion.resource_estimates import df, sf
print(f"  OK: qiskit {qiskit.__version__} | pyscf {pyscf.__version__} | resource_estimates loaded")
PY

# 5. smoke-test the estimator if it's sitting next to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "$SCRIPT_DIR/ft_resource_estimator.py" ]]; then
  say "Smoke-testing ft_resource_estimator.py (N2/STO-3G FT cost)..."
  conda run -n "$NEW_ENV" python - <<PY
import sys; sys.path.insert(0, "$SCRIPT_DIR")
from pyscf import gto, scf
from ft_resource_estimator import estimate_from_mf
mf = scf.RHF(gto.M(atom="N 0 0 0; N 0 0 1.0977", basis="sto-3g")); mf.verbose=0; mf.kernel()
r = estimate_from_mf(mf)
print(f"  lambda_DF={r['lambda_DF']:.1f}  Toffoli={r['toffoli_total']:.2e}  logical_qubits={r['logical_qubits']}")
PY
else
  say "ft_resource_estimator.py not found next to this script -- skipping smoke test."
fi

say "Done. Activate with:  conda activate $NEW_ENV"
say "Your '$SRC_ENV' env was not modified."
