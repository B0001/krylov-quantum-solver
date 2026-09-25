#!/usr/bin/env bash
# Production driver for SPEC_hubbard_bethe §10 (bead chem-tjr). Resumable: finished rows are skipped.
set -u
for U in 4 2 8; do
  uv run python benchmark_hubbard_lieb_wu.py --route open --U $U --L 20,40,60,80,100 --bond-dims 200,400,800
  uv run python benchmark_hubbard_lieb_wu.py --route ring --U $U --L 16,20,24,28,32 --bond-dims 400,800,1600
done
echo ALL_DONE
