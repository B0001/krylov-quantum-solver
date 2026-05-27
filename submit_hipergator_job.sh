#!/usr/bin/env bash
# ##############################################################################
# SLURM BATCH SUBMIT HARNESS: submit_hipergator_job.sh
# Automates multi-node orchestration & geometric sweeps on HiPerGator clusters.
# ##############################################################################
#
#SBATCH --job-name=hybrid_q_shifter      # Unique job identifier in queue
#SBATCH --output=logs/sweep_%j.out       # Standard output log (%j appends JobID)
#SBATCH --error=logs/sweep_%j.err        # Standard error log
#SBATCH --mail-type=END,FAIL             # Status notifications update trigger
#SBATCH --mail-user=your.email@ufl.edu   # Target routing email destination
#
#SBATCH --nodes=1                        # Request exactly 1 compute node
#SBATCH --ntasks=1                       # Single parent execution task thread
#SBATCH --cpus-per-task=4                # Allocate 4 dedicated CPU cores for processing
#SBATCH --mem=16gb                       # Allocate 16 GB of system RAM
#SBATCH --time=02:00:00                  # Strict 2-hour wall-clock time limit
#SBATCH --partition=hpg                  # Direct job to standard HiPerGator partition
#
# ##############################################################################

# Fail immediately if any intermediate pipeline segment returns an error code
set -e

echo "================================################################================"
echo "INITIALIZING HIPERGATOR COMPUTE SLOT DISPATCH"
echo "================================================################================"
echo "Job Commenced on Node: $(hostname)"
echo "Current Allocation Directory: ${SLURM_SUBMIT_DIR}"
echo "--------------------------------------------------------------------------------"

# 1. Move into the designated cluster submission path context
cd "${SLURM_SUBMIT_DIR}"

# 2. Build explicit logs directory tracking matrix if not yet present
mkdir -p logs

# 3. Load native cluster module tools and switch to your specific conda workspace
echo "[STEP 1] Initializing system environments and application dependencies..."
module purge
module load conda/24.1.2  # Loads HiPerGator's base module package manager wrapper

# Activate your verified chemistry conda channel workspace
# If your environment is isolated within personal user home space paths:
conda activate chem

# 4. Enforce strict PYTHONPATH declaration 
# This ensures that headless nodes can instantly discover the edited package root
export PYTHONPATH="${SLURM_SUBMIT_DIR}:${PYTHONPATH}"

# 5. Launch the automated geometric Potential Energy Surface verification sweep
echo "[STEP 2] Launching automated potential energy surface calculation tracks..."
chmod +x run_sprint_benchmarks.sh
./run_sprint_benchmarks.sh

# 6. Execute data visualization rendering loops once data points compile
if [ -f "benchmark_results.csv" ]; then
    echo "[STEP 3] Exporting compiled datasets to publication graph metrics..."
    python plot_pes_curve.py
fi

echo ""
echo "================================================================================"
echo "[SUCCESS] Slurm processing lifecycle completed smoothly."
echo "================================================################================"