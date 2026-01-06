#!/bin/bash
set -euo pipefail

# -------------------------------------------------
# PATHS
# -------------------------------------------------
BASE_DIR="/media/betselotz/Expansion/MARV-Gen"

MSA="$BASE_DIR/results/10_msa/MARV.A.1/marburg_aligned.fasta"
TREE="$BASE_DIR/results/11_phylogeny/MARV.A.1/marburg_ml.treefile"
DATES="$BASE_DIR/metadata/all_seq_metadata.csv"
OUTDIR="$BASE_DIR/results/12_treetime/MARV.A.1"

# -------------------------------------------------
# CREATE OUTPUT DIRECTORY
# -------------------------------------------------
mkdir -p "$OUTDIR"

# -------------------------------------------------
# INITIALIZE CONDA
# -------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate treetime_env

# -------------------------------------------------
# SANITY CHECKS
# -------------------------------------------------
echo ">>> Checking input files..."
for f in "$MSA" "$TREE" "$DATES"; do
    if [ ! -s "$f" ]; then
        echo "❌ ERROR: File not found or empty: $f"
        exit 1
    fi
done
echo "✅ All input files found and valid."

# -------------------------------------------------
# RUN TREETIME
# -------------------------------------------------
echo ">>> Running TreeTime time-resolved phylogeny..."

treetime \
  --aln "$MSA" \
  --tree "$TREE" \
  --dates "$DATES" \
  --outdir "$OUTDIR" \
  --reroot oldest \
  --clock-rate 0.00022 \
  --branch-length-mode input \
  --clock-filter 0 \
  --gtr JC69

# -------------------------------------------------
# CLEAN UP
# -------------------------------------------------
conda deactivate

echo "------------------------------------------------"
echo "TreeTime analysis complete"
echo "Results saved in: $OUTDIR"
echo "------------------------------------------------"
