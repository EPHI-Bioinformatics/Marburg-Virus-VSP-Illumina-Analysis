#!/usr/bin/env bash
# 04_postprocess.sh
set -euo pipefail

OUTDIR="../03_postprocessing"
RUN1_LOG="../02_beast/run1.log"
RUN2_LOG="../02_beast/run2.log"
RUN1_TREES="../02_beast/run1.trees"
RUN2_TREES="../02_beast/run2.trees"
BURNIN=10 # 10% burn-in (10 million states)

echo "[1/3] Checking ESS values using Tracer..."
# Running Tracer in CLI mode (assuming tracer CLI is installed)
# ESS > 200 indicates convergence
tracer-cli -burnin "$BURNIN" -log "$RUN1_LOG" -log "$RUN2_LOG" -out "$OUTDIR/ESS_summary.txt"

echo "[2/3] Combining posterior trees with LogCombiner..."
# 10% burn-in for each run
logcombiner -burnin 10000 -trees -reassociate \
    -o "$OUTDIR/marburg_combined.trees" "$RUN1_TREES" "$RUN2_TREES"

echo "[3/3] Generating MCC tree and tMRCA with TreeAnnotator..."
treeannotator -burnin 0 -heights mean \
    "$OUTDIR/marburg_combined.trees" \
    "$OUTDIR/marburg_MCC.tree"

echo "Done. MCC tree saved to $OUTDIR/marburg_MCC.tree"
echo "Extract tMRCA estimates (calendar years) from this tree."
