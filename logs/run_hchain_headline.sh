#!/bin/bash
cd /Users/benjaminhess/Downloads/chem
# pure D=400/800/1600 dataset: retire the mixed-schedule CSV (red flag 1)
if [ -f data/hchain_tdl.csv ] && [ ! -f data/hchain_tdl_preheadline.bak.csv ]; then
  mv data/hchain_tdl.csv data/hchain_tdl_preheadline.bak.csv
fi
echo "=== H-chain stage 1 (n=8..20, 4GB pool) start $(date) ==="
conda run -n chem python -u benchmark_hchain_tdl.py --protocol ramp \
  --ns 8,10,12,14,16,18,20 --bond-dims 400,800,1600 --stack-mem-gb 4 --threads 4
echo "=== H-chain stage 2 (n=22..26, 6GB pool) start $(date) ==="
conda run -n chem python -u benchmark_hchain_tdl.py --protocol ramp \
  --ns 22,24,26 --bond-dims 400,800,1600 --stack-mem-gb 6 --threads 4
echo "=== H-chain stage 3 (n=28,30, 8GB pool) start $(date) ==="
conda run -n chem python -u benchmark_hchain_tdl.py --protocol ramp \
  --ns 28,30 --bond-dims 400,800,1600 --stack-mem-gb 8 --threads 4
echo "=== HCHAIN PIPELINE DONE $(date) ==="
