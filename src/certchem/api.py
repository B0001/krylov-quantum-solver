import os
import uuid
import json
import time
import hashlib
import logging
from typing import Any, Dict, Optional, Tuple
from fastapi import Depends, FastAPI, HTTPException, Security, status
from fastapi.security import APIKeyHeader
from pydantic import BaseModel, Field
import redis
from google.cloud import storage

from certchem import check_caps, CapExceededError

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("certchem-api")

app = FastAPI(
    title="CertChem API",
    description="Certified chemistry ground-state solver with rigorous two-sided brackets.",
    version="1.0.0"
)

# Connect to Redis
REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
redis_client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT, decode_responses=True)

# GCS storage client
RESULTS_BUCKET = os.getenv("RESULTS_BUCKET")
storage_client = storage.Client() if RESULTS_BUCKET else None

# ---------------------------------------------------------------------------
# API-key auth + per-key rate limiting (ADR-0009: app-layer auth, per-key
# limits from day one). Keys are compared as SHA-256 digests so a wrong key
# costs the same time as a right one, and so no raw key is held in a set that
# could surface in a traceback or repr.
# ---------------------------------------------------------------------------
def _digest(key: str) -> str:
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


ALLOWED_KEY_DIGESTS = frozenset(
    _digest(k.strip()) for k in os.getenv("ALLOWED_API_KEYS", "").split(",") if k.strip()
)
RATE_LIMIT_PER_MIN = int(os.getenv("RATE_LIMIT_PER_MIN", "100"))
# Caller metadata is opaque and echoed to storage, so it is bounded.
METADATA_MAX_BYTES = 1024
_api_key_header = APIKeyHeader(name="X-API-Key", auto_error=False)


def require_api_key(api_key: str | None = Security(_api_key_header)) -> str:
    """Authenticate the caller by X-API-Key, then charge them a request.

    Fails CLOSED: an unconfigured ALLOWED_API_KEYS serves 503 rather than
    silently running an open endpoint. Rate limiting is a fixed one-minute
    window per key (INCR + EXPIRE), which is one Redis round-trip and cannot
    grow unboundedly -- buckets expire themselves.
    """
    if not ALLOWED_KEY_DIGESTS:
        logger.error("ALLOWED_API_KEYS is unset; refusing to serve an unauthenticated API.")
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="API authentication is not configured on this deployment.",
        )

    digest = _digest(api_key) if api_key else ""
    if digest not in ALLOWED_KEY_DIGESTS:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing or invalid X-API-Key.",
            headers={"WWW-Authenticate": "X-API-Key"},
        )

    # Fixed-window counter, keyed by digest so no raw key lands in Redis.
    window = int(time.time() // 60)
    bucket = f"ratelimit:{digest[:16]}:{window}"
    pipe = redis_client.pipeline()
    pipe.incr(bucket)
    pipe.expire(bucket, 60)
    used = int(pipe.execute()[0])

    if used > RATE_LIMIT_PER_MIN:
        raise HTTPException(
            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
            detail=f"Rate limit exceeded: {RATE_LIMIT_PER_MIN} requests/minute.",
            headers={"Retry-After": str(60 - int(time.time() % 60))},
        )

    return digest


class EnergyRequest(BaseModel):
    molecule: str = Field(..., description="Molecule geometry in PySCF/XYZ atom format (e.g. 'H 0 0 0; H 0 0 0.74')")
    basis: str = Field("sto-3g", description="Atomic basis set (e.g. 'sto-3g', '6-31g')")
    active_space: Tuple[int, int] = Field(..., description="Active space as (electrons, orbitals)")
    mode: str = Field("certified", description="Execution mode: 'certified' (returns brackets) or 'fast' (bare point-estimate)")
    krylov_dim: int = Field(12, description="Krylov subspace dimension for diagonalization (default 12, must be >= 6 for self-certified)")
    metadata: Optional[Dict[str, Any]] = Field(
        None,
        description=(
            "Opaque caller-supplied provenance (e.g. {'family': 'n2_curve', "
            "'coordinate': 1.4}), echoed verbatim into the result blob. Lets a "
            "batch client label its own points instead of inferring them later."
        ),
    )


@app.get("/health")
def health_check():
    """Verify service health and Redis connectivity. Public: GCP health checks carry no key."""
    try:
        redis_client.ping()
        return {"status": "healthy", "redis": "connected"}
    except Exception as e:
        return {"status": "unhealthy", "redis_error": str(e)}


@app.get("/v1/limits")
def get_limits():
    """Return the hard-capped validation limits of the CertChem engine. Public."""
    from certchem.limits import ALLOWED_BASES, MAX_SPIN_ORBITALS
    return {
        "max_spin_orbitals": MAX_SPIN_ORBITALS,
        "validated_bases": sorted(ALLOWED_BASES),
        "notes": "Any request exceeding these caps will be rejected with HTTP 422 immediately."
    }


@app.post(
    "/v1/energy",
    status_code=status.HTTP_202_ACCEPTED,
    dependencies=[Depends(require_api_key)],
)
def submit_energy_job(req: EnergyRequest):
    """Submit a chemistry ground-state calculation job to the queue.

    Saves resource waste by running the capabilities check (check_caps) UPFRONT
    before the job is queued. Rejects invalid requests with HTTP 422 immediately.
    """
    # 1. Validate caps and mode upfront
    if req.mode not in ("certified", "fast"):
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail="mode must be either 'certified' or 'fast'"
        )

    if req.metadata is not None and len(json.dumps(req.metadata)) > METADATA_MAX_BYTES:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"metadata must serialize to at most {METADATA_MAX_BYTES} bytes",
        )

    try:
        check_caps(req.molecule, req.basis, req.active_space)
    except CapExceededError as err:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail={"error": "CapExceededError", "message": str(err), "cap": err.cap}
        )
    except Exception as err:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Validation failed: {str(err)}"
        )

    # 2. Generate a unique job ID
    job_id = uuid.uuid4().hex
    logger.info(f"Submitting job {job_id} to queue (mode: {req.mode}, basis: {req.basis})")

    # 3. Create job record
    job_data = {
        "job_id": job_id,
        "status": "queued",
        "request": {
            "molecule": req.molecule,
            "basis": req.basis,
            "active_space": req.active_space,
            "mode": req.mode,
            "krylov_dim": req.krylov_dim,
            "metadata": req.metadata,
        }
    }

    # Save initial status in Redis with a 24-hour expiration TTL
    redis_client.setex(f"job:{job_id}", 86400, json.dumps(job_data))
    # Push job_id onto the job queue
    redis_client.lpush("certchem:jobs", job_id)

    return {"job_id": job_id, "status": "queued"}


@app.get("/v1/jobs/{job_id}", dependencies=[Depends(require_api_key)])
def get_job_status(job_id: str):
    """Retrieve the status and results of a submitted job.

    Pulls the status from Redis. If the job is completed, it returns the result.
    If the result is in GCS, it streams it directly.
    """
    raw_job = redis_client.get(f"job:{job_id}")
    if not raw_job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Job not found or expired from Redis cache."
        )

    job_data = json.loads(raw_job)
    current_status = job_data.get("status")

    if current_status == "completed":
        # Check if the full result is stored directly in Redis
        if "result" in job_data:
            return job_data

        # Otherwise, retrieve the result from the GCS bucket
        if RESULTS_BUCKET and storage_client:
            try:
                bucket = storage_client.bucket(RESULTS_BUCKET)
                blob = bucket.blob(f"results/{job_id}.json")
                if blob.exists():
                    gcs_data = json.loads(blob.download_as_text())
                    job_data["result"] = gcs_data
                    return job_data
            except Exception as e:
                logger.error(f"Failed to fetch result for job {job_id} from GCS: {e}")
                raise HTTPException(
                    status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail="Job is marked as completed but its GCS result could not be retrieved."
                )

    return job_data
