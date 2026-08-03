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
# This is also why the Kubernetes Jobs invoke bare `python`: on PATH it IS the
# venv interpreter, so the CLAUDE.md "always `uv run`" rule is already satisfied
# -- there is no second interpreter in the image to pick up by accident.
ENV PATH="/app/.venv/bin:$PATH"

# Bytecode is baked at build time (UV_COMPILE_BYTECODE); writing it at runtime
# would need a writable /app, which readOnlyRootFilesystem forbids.
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Thread count is deliberately NOT pinned here, unlike the sibling repos. These
# benchmarks are single-process and PySCF/BLAS-bound, so they genuinely want
# multiple threads -- but the count must match the pod's CPU limit, not the
# node's core count, or the runtime oversubscribes and thrashes. The Kubernetes
# Jobs set OMP_NUM_THREADS to their own cpu limit; this default keeps a plain
# `docker run` single-threaded and predictable rather than silently grabbing
# every core on the host.
ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# The account's home is /home/app, which does not exist and is on the read-only
# root filesystem anyway. Anything that resolves $HOME to write a dotfile --
# PySCF's ~/.pyscf_conf.py, matplotlib's cache -- fails or warns without this.
# /tmp is the emptyDir the Jobs mount.
ENV HOME=/tmp

# Headless plotting. Without MPLCONFIGDIR, matplotlib tries $HOME/.config on a
# read-only filesystem, prints a warning, and rebuilds its font cache on every
# single import -- measurably slow in a benchmark loop. Agg because there is no
# display; plot_benchmarks.py and plot_pes_curve.py both render to file.
ENV MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib

# Non-root runtime. Kubernetes pins runAsUser/fsGroup to this same 10001 so the
# results PVC is writable; keep the two in sync if either changes.
RUN groupadd --gid 10001 app \
    && useradd --uid 10001 --gid 10001 --no-create-home --shell /usr/sbin/nologin app

# Every benchmark writes to a relative `data/...` path resolved against /app
# (e.g. benchmark_n2.py's data/n2_dissociation.csv). .dockerignore excludes
# data/, so the directory does not otherwise exist -- create it as the mount
# point for the PVC. .dmrg_tmp is block2's scratch dir, mounted as an emptyDir.
RUN mkdir -p /app/data /app/results /app/.dmrg_tmp \
    && chown 10001:10001 /app/data /app/results /app/.dmrg_tmp

USER 10001:10001

# Import smoke check by default: cheap, and it fails loudly if the qiskit-nature
# import chain is broken (the scipy pin exists precisely because it can be).
CMD ["python", "-c", "import hybrid_quantum_solver as h; print('ok:', len(h.__all__), 'exports')"]
