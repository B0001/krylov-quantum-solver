# CertChem: one image, two entrypoints (api | worker) — ADR-0009
FROM python:3.12-slim AS base
ENV PYTHONUNBUFFERED=1 \
    # Pin math libraries to one thread per process; parallelism is process-level
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY src/ ./src/
COPY golden_suite/ ./golden_suite/

# Non-root runtime
RUN useradd --create-home appuser
USER appuser

# API entrypoint (default). Cloud Run injects $PORT; bind 0.0.0.0, never loopback.
CMD ["sh", "-c", "uvicorn src.certchem.api:app --host 0.0.0.0 --port ${PORT:-8080}"]

# Worker entrypoint (Cloud Run job overrides CMD):
#   python -m src.certchem.worker
