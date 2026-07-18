#!/usr/bin/env bash
# run_gates.sh — parallel, cached runner for spec acceptance gates.
#
# Each gate file runs in its OWN process (block2's OpenMP runtime segfaults if it
# shares an interpreter with pyscf/qiskit-aer — same constraint as the original
# sequential `make gates`). Parallelism is across processes, so isolation holds.
#
# Cache soundness: a gate is skipped only if it passed with an IDENTICAL global
# source state. The cache key = sha256(all tracked *.py) + sha256(gate file).
# This is deliberately conservative: ANY .py change invalidates ALL cached passes.
# It therefore only saves time on re-runs after no-code changes (docs, specs,
# repeated nightly runs) — but it can never wrongly skip a gate.
#
# Env:
#   GATE_RUN      command prefix (default: uv run)  e.g. GATE_RUN=""
#   GATE_JOBS     parallel processes (default: nproc)
#   GATE_NO_CACHE 1 = ignore cache (still records passes)
#   GATE_GLOB     gate file glob (default: tests/test_*_spec.py)
set -u

RUN=${GATE_RUN-uv run}
JOBS=${GATE_JOBS:-$(nproc 2>/dev/null || echo 4)}
GLOB=${GATE_GLOB:-tests/test_*_spec.py}
CACHE_DIR=.gates-cache
LOG_DIR=logs/gates
mkdir -p "$CACHE_DIR" "$LOG_DIR"

# shellcheck disable=SC2086
FILES=$(ls $GLOB 2>/dev/null || true)
if [ -z "$FILES" ]; then
    echo "no spec gate tests match: $GLOB"
    exit 0
fi

if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    SRC_DIGEST=$(git ls-files '*.py' | sort | xargs sha256sum 2>/dev/null | sha256sum | cut -d' ' -f1)
else
    # Outside git: no safe way to fingerprint sources -> disable caching.
    SRC_DIGEST="nogit-$(date +%s)"
fi
export RUN SRC_DIGEST CACHE_DIR LOG_DIR
export GATE_NO_CACHE=${GATE_NO_CACHE:-0}

run_gate() {
    f=$1
    name=$(basename "$f" .py)
    key=$( { echo "$SRC_DIGEST"; sha256sum "$f"; } | sha256sum | cut -d' ' -f1)
    cache_file="$CACHE_DIR/$name.$key"
    if [ -f "$cache_file" ] && [ "$GATE_NO_CACHE" != "1" ]; then
        echo "SKIP  $f  (cached pass, sources unchanged)"
        return 0
    fi
    start=$(date +%s)
    if $RUN python -m pytest "$f" -q >"$LOG_DIR/$name.log" 2>&1; then
        rm -f "$CACHE_DIR/$name".* 2>/dev/null
        touch "$cache_file"
        echo "PASS  $f  ($(( $(date +%s) - start ))s)"
        return 0
    else
        echo "FAIL  $f  ($(( $(date +%s) - start ))s)  log: $LOG_DIR/$name.log"
        return 1
    fi
}
export -f run_gate

echo "gates: $(echo "$FILES" | wc -w | tr -d ' ') files, $JOBS parallel processes, cache=$([ "$GATE_NO_CACHE" = "1" ] && echo off || echo on)"
echo "$FILES" | tr ' ' '\n' | xargs -P "$JOBS" -I{} bash -c 'run_gate "$@"' _ {}
STATUS=$?

if [ $STATUS -eq 0 ]; then
    echo "all spec gates passed"
else
    echo "GATE FAILURES — see logs above. Full output per gate in $LOG_DIR/"
    exit 1
fi
