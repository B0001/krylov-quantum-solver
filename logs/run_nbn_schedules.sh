#!/bin/bash
cd /Users/benjaminhess/Downloads/chem
echo "=== NbN schedule A (perD 400/800/1200) start $(date) ==="
conda run -n chem python -u logs/nbn_schedule.py A
echo "=== NbN schedule B (ramp 300/600/1200) start $(date) ==="
conda run -n chem python -u logs/nbn_schedule.py B
echo "=== NBN PIPELINE DONE $(date) ==="
