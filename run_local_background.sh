#!/usr/bin/env bash
# =#############################################################################
# LOCAL BACKGROUND DRIVER: run_local_background.sh
# Bypasses Slurm requirements to run automated sweeps asynchronously on macOS.
# ##############################################################################

export LANG=en_US.UTF-8

# Generate dedicated tracking channels for your asynchronous task logs
mkdir -p logs
LOG_OUT="logs/local_sweep_run.log"

echo "================================================================================"
echo "INITIALIZING LOCAL ASYNCHRONOUS SWEEP OPERATIONS"
echo "================================================================================"
echo "Bypassing cluster managers. Task routing directly to internal hardware threads."
echo "Active Telemetry Stream Destination: ${LOG_OUT}"
echo "--------------------------------------------------------------------------------"

# Execute your pipeline sweep wrapper as an isolated background process fork
{
    echo "=== BACKGROUND LIFE CYCLE COMMENCED: $(date) ==="
    chmod +x run_sprint_benchmarks.sh
    ./run_sprint_benchmarks.sh
    
    echo -e "\n=== CHANNELS CONCLUDED. DRAWING VISUAL GRAPH METRICS ==="
    python plot_pes_curve.py
    echo "=== LIFECYCLE COMPLETED SUCCESSFULLY: $(date) ==="
} > "${LOG_OUT}" 2>&1 &

# Capture the unique background job Process ID assigned by your Zsh shell
BACKGROUND_PID=$!

echo "  -> Asynchronous process successfully initialized."
echo "  -> Assigned Local Task Thread PID: ${BACKGROUND_PID}"
echo ""
echo "You can check calculation progress at any time by monitoring the log tail:"
echo "    tail -f ${LOG_OUT}"
echo ""
echo "To terminate the running execution task mid-loop, run:"
echo "    kill ${BACKGROUND_PID}"
echo "================================================================================"