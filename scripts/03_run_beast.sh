#!/usr/bin/env bash
# 03_run_beast.sh
set -euo pipefail

# --- CONFIGURATION ---
XML="../02_beast/marburg_hybrid_model.xml"
THREADS=8
BEAGLE_SSE="true"

echo "Starting BEAST Run 1..."
beast -beagle_SSE -threads "$THREADS" -seed 1111 -overwrite -prefix "../02_beast/run1" "$XML" > "../02_beast/run1.log" 2>&1 &

echo "Starting BEAST Run 2..."
beast -beagle_SSE -threads "$THREADS" -seed 2222 -overwrite -prefix "../02_beast/run2" "$XML" > "../02_beast/run2.log" 2>&1 &

wait
echo "Both BEAST runs completed."

