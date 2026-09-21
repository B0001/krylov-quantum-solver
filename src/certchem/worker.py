import os
import json
import logging
import time
import traceback
import redis
from google.cloud import storage

# Ensure thread pinning is applied immediately
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

# Import after setting env vars
from certchem import certified_energy, Mode, CertifiedResult
from certchem.contract import FloorViolationError, CapExceededError, ConvergenceError

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("certchem-worker")

# Connect to Redis
REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
redis_client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT, decode_responses=True)

# GCS storage client
RESULTS_BUCKET = os.getenv("RESULTS_BUCKET")
storage_client = storage.Client() if RESULTS_BUCKET else None


def upload_result_to_gcs(job_id: str, result_data: dict):
    """Upload completed job results to GCS bucket as a JSON blob."""
    if not RESULTS_BUCKET or not storage_client:
        logger.info("GCS results bucket not configured; skipping upload.")
        return

    try:
        bucket = storage_client.bucket(RESULTS_BUCKET)
        blob = bucket.blob(f"results/{job_id}.json")
        blob.upload_from_string(json.dumps(result_data), content_type="application/json")
        logger.info(f"Successfully uploaded result for job {job_id} to GCS.")
    except Exception as e:
        logger.error(f"Failed to upload result for job {job_id} to GCS: {e}")


def process_job(job_id: str):
    """Process a single chemistry energy job from the queue."""
    logger.info(f"Processing job {job_id}")

    # 1. Retrieve job details
    raw_job = redis_client.get(f"job:{job_id}")
    if not raw_job:
        logger.error(f"Job {job_id} details missing from Redis.")
        return

    job_data = json.loads(raw_job)
    request = job_data["request"]

    # 2. Update status to running
    job_data["status"] = "running"
    job_data["started_at"] = time.time()
    redis_client.set(f"job:{job_id}", json.dumps(job_data))

    try:
        # 3. Parse parameters
        molecule = request["molecule"]
        basis = request["basis"]
        active_space = tuple(request["active_space"])
        mode = Mode.FAST if request["mode"] == "fast" else Mode.CERTIFIED
        krylov_dim = request["krylov_dim"]

        # 4. Run the certified solver (pyscf, qiskit nature, temple bounds, etc.)
        t0 = time.time()
        result = certified_energy(
            molecule=molecule,
            basis=basis,
            cas=active_space,
            mode=mode,
            krylov_dim=krylov_dim
        )
        duration = time.time() - t0

        # 5. Format output.
        # Provenance first: the blob records WHAT WAS RUN, stated by the
        # producer. A consumer must never have to infer the geometry back out
        # of the energy -- doing that against a reference makes the reference
        # useless as a check (it can no longer disagree).
        provenance = {
            "job_id": job_id,
            "molecule": molecule,
            "basis": basis,
            "active_space": list(active_space),
            "mode": request["mode"],
            "krylov_dim": krylov_dim,
            "metadata": request.get("metadata"),
        }
        serialized_result = {}
        if isinstance(result, CertifiedResult):
            serialized_result = {
                "provenance": provenance,
                "best_estimate_hartree": result.bracket.best_estimate_hartree,
                "lower_bound_hartree": result.bracket.lower_hartree,
                "upper_bound_hartree": result.bracket.upper_hartree,
                "bracket_width_hartree": result.bracket.width,
                "certificate": {
                    "method": result.certificate.method,
                    "floor_check": result.certificate.floor_check,
                    "krylov_dim": result.certificate.krylov_dim,
                    "convergence": result.certificate.convergence,
                    "solver_version": result.certificate.solver_version,
                    "manifest": result.certificate.manifest
                }
            }
        else:
            # Mode.FAST returns a bare float
            serialized_result = {
                "provenance": provenance,
                "best_estimate_hartree": float(result),
                "certificate": {
                    "method": "fast_point_estimate",
                    "floor_check": "skipped",
                    "krylov_dim": krylov_dim
                }
            }

        # Save success metadata
        job_data["status"] = "completed"
        job_data["completed_at"] = time.time()
        job_data["wall_time_s"] = round(duration, 3)
        serialized_result["wall_time_s"] = round(duration, 3)
        job_data["result"] = serialized_result

        # Upload result blob to GCS
        upload_result_to_gcs(job_id, serialized_result)

    except FloorViolationError as err:
        logger.error(f"Variational floor violation in job {job_id}: {err}")
        job_data["status"] = "failed"
        job_data["error"] = {
            "type": "FloorViolationError",
            "message": str(err),
            "diagnostics": err.diagnostics
        }
    except CapExceededError as err:
        logger.error(f"Capability cap exceeded in job {job_id}: {err}")
        job_data["status"] = "failed"
        job_data["error"] = {
            "type": "CapExceededError",
            "message": str(err),
            "cap": err.cap
        }
    except ConvergenceError as err:
        logger.error(f"Solver convergence failed in job {job_id}: {err}")
        job_data["status"] = "failed"
        job_data["error"] = {
            "type": "ConvergenceError",
            "message": str(err),
            "partial": err.partial
        }
    except Exception as err:
        logger.error(f"Unexpected error in job {job_id}: {err}\n{traceback.format_exc()}")
        job_data["status"] = "failed"
        job_data["error"] = {
            "type": type(err).__name__,
            "message": str(err),
            "traceback": traceback.format_exc()
        }

    # Store back to Redis with a 24-hour TTL (86400 seconds)
    redis_client.setex(f"job:{job_id}", 86400, json.dumps(job_data))
    logger.info(f"Finished job {job_id} with status: {job_data['status']}")


def main():
    logger.info("Initializing CertChem Redis queue consumer worker...")
    logger.info(f"Connected to Redis host: {REDIS_HOST}:{REDIS_PORT}")
    if RESULTS_BUCKET:
        logger.info(f"GCS bucket for results: {RESULTS_BUCKET}")
    else:
        logger.warning("GCS bucket NOT configured; results will only reside in Redis memory.")

    # Start main BRPOP worker loop
    try:
        while True:
            try:
                # Block until a job is pushed onto the queue
                # BRPOP returns a tuple: (list_name, element_value)
                raw_element = redis_client.brpop("certchem:jobs", timeout=10)
                if raw_element:
                    _, job_id = raw_element
                    try:
                        process_job(job_id)
                    except Exception as e:
                        logger.error(f"Unhandled exception during job execution: {e}")
            except (redis.exceptions.TimeoutError, TimeoutError):
                # Standard empty queue timeout from network/Memorystore idle policies
                pass
            except redis.exceptions.ConnectionError as e:
                logger.warning(f"Redis connection interrupted: {e}. Reconnecting in 5s...")
                time.sleep(5)
            except Exception as e:
                logger.error(f"Queue retrieval error: {e}")
                time.sleep(2)
    except KeyboardInterrupt:
        logger.info("Worker terminated by user (KeyboardInterrupt). Exiting...")


if __name__ == "__main__":
    main()
