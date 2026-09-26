# Worker image for chem's sandbox loop ([tool.sandbox] image = "claude-chem").
#
# block2 0.5.3 (the DMRG backend) ships no Linux arm64 wheel, and building it from
# source on arm64 fails (pybind11 + CMake/BLAS configuration). Its manylinux x86_64
# wheel installs cleanly, so chem's workers run as linux/amd64: Docker Desktop
# emulates them on Apple silicon (slower, but no compiling).
#
# These are the shared `claude` image's own steps (see `docker history claude`),
# rebuilt for amd64. Keep them in step with that image when it changes.
#   docker build -f sandbox.Dockerfile -t claude-chem .
# uv must be the amd64 build too: an arm64 uv runs natively in Docker's arm64 VM even
# inside an amd64 container, then fetches arm64 Python and asks for block2's missing wheel.
FROM --platform=linux/amd64 ghcr.io/astral-sh/uv:latest AS uv

FROM --platform=linux/amd64 node:22-bookworm-slim

RUN apt-get update && apt-get install -y git curl && rm -rf /var/lib/apt/lists/*
WORKDIR /workspace
RUN npm install -g @anthropic-ai/claude-code
COPY --from=uv /uv /uvx /usr/local/bin/

ARG BD_VERSION=1.1.2
RUN arch="$(dpkg --print-architecture)" \
    && f="beads_${BD_VERSION}_linux_${arch}.tar.gz" \
    && cd /tmp \
    && curl -fsSLO "https://github.com/gastownhall/beads/releases/download/v${BD_VERSION}/${f}" \
    && curl -fsSLO "https://github.com/gastownhall/beads/releases/download/v${BD_VERSION}/checksums.txt" \
    && grep " ${f}\$" checksums.txt | sha256sum -c - \
    && tar -xzf "$f" -C /usr/local/bin bd \
    && rm -f "$f" checksums.txt \
    && bd version

USER node
RUN bd metrics off
RUN git config --global beads.role contributor
ENTRYPOINT ["claude"]
