# syntax=docker/dockerfile:1
#
# Reproducible image for the hybrid quantum-classical chemistry pipeline.
#
# Python is pinned to 3.13, NOT latest: qiskit-aer publishes no cp314 wheel, so on 3.14 uv falls
# back to a scikit-build source build and the image fails to build at all. 3.13 also matches the
# local .venv, so container and laptop resolve the same lock.
#
#   docker build -t chem .
#   docker run --rm chem
#
# The DMRG (block2) and test extras are NOT installed by default. block2 initializes its own
# OpenMP runtime and segfaults if it loads into a process that already imported pyscf or
# qiskit-aer (see CLAUDE.md), so it belongs in a deliberately separate image or an explicit
# build-time --extra, not in the default runtime.

# --- Stage 1: resolve and build the environment ---
FROM python:3.13-slim-bookworm AS builder

# uv from the official distroless image -- pinned, no curl|sh bootstrap.
COPY --from=ghcr.io/astral-sh/uv:0.6.9 /uv /usr/local/bin/uv

ENV UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    UV_PYTHON_DOWNLOADS=never \
    UV_PYTHON=python3.13

WORKDIR /app

# Dependency layer first, so editing source does not invalidate the (slow) dependency solve.
# --no-install-project: the project itself is installed by the second sync, after the source copy.
COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-install-project --no-dev

# Then the source, then the project itself into the same venv.
# Explicit paths, never `COPY * .`: a multi-source glob copies each directory's CONTENTS into the
# destination, which would flatten hybrid_quantum_solver/ into /app and destroy the package.
COPY hybrid_quantum_solver/ ./hybrid_quantum_solver/
COPY tests/ ./tests/
COPY specs/ ./specs/
COPY *.py ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev

# --- Stage 2: lightweight runtime ---
FROM python:3.13-slim-bookworm AS runtime

WORKDIR /app

# The prepared venv, then the source it was installed from.
COPY --from=builder /app/.venv /app/.venv
COPY --from=builder /app/hybrid_quantum_solver/ ./hybrid_quantum_solver/
COPY --from=builder /app/tests/ ./tests/
COPY --from=builder /app/specs/ ./specs/
COPY --from=builder /app/*.py ./
COPY pyproject.toml ./

# Prepend the venv so `python` and console scripts resolve without activation.
ENV PATH="/app/.venv/bin:$PATH"

CMD ["python", "-c", "import hybrid_quantum_solver as h; print('ok:', len(h.__all__), 'exports')"]
