#!/usr/bin/env bash
# 01_temporal_signal.sh
set -euo pipefail

# --- CONFIGURATION ---
ALIGNMENT="../01_temporal_signal/MARV_global_aligned.fasta"
OUTGROUP="RAVV"
OUTDIR="../01_temporal_signal"
METADATA="../01_temporal_signal/dates.tsv" # Format: ID \t date (decimal)
THREADS=8

echo "[1/3] Building Maximum Likelihood tree (IQ-TREE)..."
# ModelFinder (BIC) + 1000 ultrafast bootstraps
iqtree2 -s "$ALIGNMENT" -m MFP -B 1000 -T "$THREADS" --prefix "$OUTDIR/MARV_ML"

echo "[2/3] Pruning outgroup ($OUTGROUP) using TreeTime/SeqKit..."
# Prune RAVV from the tree to restrict analysis to global MARV genomes
treetime prune --tree "$OUTDIR/MARV_ML.treefile" --outgroup "$OUTGROUP" --outtree "$OUTDIR/MARV_ML_pruned.treefile"

echo "[3/3] Preparing TempEst Root-to-Tip Regression..."
# TempEst requires a tree + metadata file.
# Note: TempEst is a GUI tool. For CLI, we generate the exact input file, 
# then open TempEst manually to optimize the root and extract the R^2.
cat << EOF > "$OUTDIR/TempEst_Instructions.txt"
--- TEMPEST INSTRUCTIONS ---
1. Open TempEst v1.5.3.
2. Load Tree: $OUTDIR/MARV_ML_pruned.treefile
3. Load Dates: $METADATA
4. Click "Parse Dates". 
5. Under the "Root-to-tip" tab, click "Find Best Root".
6. Record the Correlation Coefficient (R^2) and the root date.
7. Save the rerooted tree as: $OUTDIR/MARV_ML_rerooted.treefile
-----------------------------
EOF
echo "Done. Please follow instructions in $OUTDIR/TempEst_Instructions.txt"
