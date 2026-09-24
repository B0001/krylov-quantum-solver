#!/usr/bin/env bash
# ==============================================================================
# CERTCHEM GRAND CLOUD RUNNER: run_grand_sweep_on_cloud.sh
# Automates the complete GCP lifecycle for your 7-family, 113-point sweep.
# Utilizes EXIT traps to guarantee resource destruction under all conditions.
# ==============================================================================

set -euo pipefail

# Ensure target project ID is set before provisioning
export GCP_PROJECT_ID="${GCP_PROJECT_ID:-tlsmbl-compute}"
export GCP_IMAGE_URL="${GCP_IMAGE_URL:-gcr.io/tlsmbl-compute/certchem:latest}"
export GCP_REGION="${GCP_REGION:-us-central1}"

DEPLOY_DIR="architecture/deployment"
ORIGINAL_DIR="$PWD"

echo "================================================================================"
echo "INITIALIZING CERTCHEM GRAND GCP SWEEP RUNNER"
echo "================================================================================"
echo "Target GCP Project: ${GCP_PROJECT_ID}"
echo "Target Image:       ${GCP_IMAGE_URL}"
echo "--------------------------------------------------------------------------------"

# --- THE SAFETY NET: EXIT TRAP ------------------------------------------------
# This block is GUARANTEED to execute when this script finishes, fails, or is aborted.
cleanup() {
    echo ""
    echo "================================================================================"
    echo "[CLEANUP] TRIGGERING GUARANTEED GCP INFRASTRUCTURE TEARDOWN"
    echo "================================================================================"
    cd "${ORIGINAL_DIR}"
    if [ -d "${DEPLOY_DIR}" ]; then
        cd "${DEPLOY_DIR}"
        echo "Removing GCS results bucket from state to preserve your data..."
        terragrunt state rm google_storage_bucket.results || echo "Bucket already removed or not found in state."
        echo "Executing terragrunt destroy..."
        terragrunt destroy -auto-approve || echo "WARNING: Teardown failed or was already destroyed."
    fi
    echo "Teardown completed. Billing clock stopped."
}
trap cleanup EXIT
# ------------------------------------------------------------------------------

# 1. Provision GCP infrastructure
echo "[STEP 1] Provisioning VPC connector, Redis queue, and Cloud Run service..."
cd "${DEPLOY_DIR}"
echo "Synchronizing GCS results bucket state..."
terragrunt import google_storage_bucket.results "${GCP_PROJECT_ID}-certchem-results" || echo "Bucket already in state or skipping import..."
terragrunt apply -auto-approve

# 2. Extract API Endpoint URL
echo "[STEP 2] Acquiring active API gateway endpoint..."
API_URL=$(terragrunt output -raw api_url)
echo "Active endpoint: ${API_URL}"

# Return to root directory
cd "${ORIGINAL_DIR}"

# 3. Execute the high-throughput Python sweep orchestrator
echo "[STEP 3] Running the High-Throughput Grand Chemistry Sweep Suite..."
conda run -n chem python run_grand_sweep.py "${API_URL}"

echo ""
echo "Grand sweep execution finished successfully."
