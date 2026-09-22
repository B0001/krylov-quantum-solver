#!/usr/bin/env bash
# ==============================================================================
# CERTCHEM CLOUD RUNNER: run_on_cloud.sh
# Automates the entire GCP lifecycle: Spin-up -> Run -> Extract -> Teardown.
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
echo "INITIALIZING CERTCHEM GCP LIFECYCLE RUNNER"
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

# 3. Create a Python runner to submit and poll the heavy cc-pVQZ Be2 calculation
echo "[STEP 3] Running the Be2 cc-pVQZ Complete Basis Set (CBS) experiment..."
conda run -n chem python -c "
import urllib.request
import json
import time
import sys

api_url = '${API_URL}'.rstrip('/')
headers = {'Content-Type': 'application/json'}

# Define 13 Beryllium geometry sweep points (cc-pVQZ is extremely large for laptops!)
R_points = [2.0, 2.1, 2.2, 2.3, 2.4, 2.45, 2.6, 2.7, 2.8, 3.0, 4.0, 6.0, 8.0]
job_ids = []

print(f'Submitting {len(R_points)} geometry points to Cloud Run worker pool...')
for R in R_points:
    req_data = {
        'molecule': f'Be 0 0 0; Be 0 0 {R}',
        'basis': 'cc-pvdz',
        'active_space': [4, 8],
        'mode': 'certified',
        'krylov_dim': 10
    }
    
    req = urllib.request.Request(
        f'{api_url}/v1/energy',
        data=json.dumps(req_data).encode('utf-8'),
        headers=headers,
        method='POST'
    )
    
    try:
        with urllib.request.urlopen(req) as resp:
            resp_data = json.loads(resp.read().decode('utf-8'))
            job_id = resp_data['job_id']
            job_ids.append((R, job_id))
            print(f'  -> R={R:4.2f} A submitted | JobID: {job_id}')
    except Exception as e:
        print(f'  -> [ERROR] Failed to submit point R={R}: {e}')
        sys.exit(1)

# Trigger the Cloud Run queue worker job
import subprocess
print('\nTriggering the Cloud Run queue worker job to process the queued calculations...')
try:
    subprocess.run([
        'gcloud', 'run', 'jobs', 'execute', 'certchem-worker',
        '--region', 'us-central1', '--quiet'
    ], check=True, stdout=subprocess.DEVNULL)
    print('Queue worker job triggered successfully!')
except Exception as e:
    print(f'Warning: Failed to trigger queue worker job: {e}')

# Poll the API for async worker queue completion
completed_results = []
pending_jobs = list(job_ids)

print('\nEntering queue polling loop. Workers are executing calculations on Cloud Run...')
while pending_jobs:
    time.sleep(10)
    still_pending = []
    
    for R, job_id in pending_jobs:
        req = urllib.request.Request(f'{api_url}/v1/jobs/{job_id}', method='GET')
        try:
            with urllib.request.urlopen(req) as resp:
                job_status = json.loads(resp.read().decode('utf-8'))
                status = job_status.get('status')
                
                if status == 'completed':
                    res = job_status['result']
                    completed_results.append({
                        'R': R,
                        'best_estimate_hartree': res['best_estimate_hartree'],
                        'lower_bound_hartree': res['lower_bound_hartree'],
                        'upper_bound_hartree': res['upper_bound_hartree'],
                        'wall_time_s': job_status.get('wall_time_s', 0.0)
                    })
                    print(f'  -> [COMPLETE] R={R:4.2f} A | E={res[\"best_estimate_hartree\"]:.6f} Ha (in {job_status.get(\"wall_time_s\", 0.0)}s)')
                elif status == 'failed':
                    print(f'  -> [FAILED] R={R:4.2f} A | Error: {job_status.get(\"error\", {}).get(\"message\")}')
                else:
                    still_pending.append((R, job_id))
        except Exception as e:
            print(f'  -> [WARNING] Error checking job {job_id}: {e}')
            still_pending.append((R, job_id))
            
    pending_jobs = still_pending
    if pending_jobs:
        print(f'Checking queue... {len(pending_jobs)} calculations still running in cloud...')

# Save local CSV
if completed_results:
    output_path = 'data/be2_cbs_cloud_curve.csv'
    import csv
    import os
    os.makedirs('data', exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=completed_results[0].keys())
        w.writeheader()
        w.writerows(completed_results)
    print(f'\n[SUCCESS] All cloud calculations completed! Saved local curve data to: {output_path}')
else:
    print('\n[ERROR] No calculations completed successfully.')
"

echo ""
echo "Lifecycle script execution finished successfully."
