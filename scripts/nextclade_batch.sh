#!/bin/bash
set -euo pipefail

# ===============================
# Nextclade batch pipeline
# ===============================

# -------------------------------
# Auto-detect available CPU threads
# -------------------------------
if command -v nproc &>/dev/null; then
    THREADS=$(nproc)
else
    THREADS=4  # fallback if nproc is unavailable
fi

# -------------------------------
# Directories
# -------------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONSENSUS_DIR="$PROJECT_DIR/results/07_consensus"
COMPARE_DIRS="$PROJECT_DIR/reference_genomes/MARV_compare/human/cleaned $PROJECT_DIR/reference_genomes/MARV_compare/bat/cleaned"
DATASET_DIR="$PROJECT_DIR/database/nextclade_marburg_dataset"
OUTPUT_DIR="$PROJECT_DIR/results/09_nextclade"

# -------------------------------
# Activate conda environment
# -------------------------------
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate nextclade_env || { echo "ERROR: Could not activate conda environment 'nextclade_env'"; exit 1; }
else
    echo "ERROR: Conda not found. Please install Miniconda or Anaconda."
    exit 1
fi

# -------------------------------
# Info
# -------------------------------
echo "🧬 Running Nextclade batch pipeline"
echo "📁 Consensus dir  : $CONSENSUS_DIR"
echo "📁 Compare dirs   : $COMPARE_DIRS"
echo "📁 Local dataset  : $DATASET_DIR"
echo "📁 Output dir     : $OUTPUT_DIR"
echo "🧵 Using threads  : $THREADS"
echo "================================================="

# -------------------------------
# Collect all FASTA files as positional arguments
# -------------------------------
FASTA_FILES=$(find $CONSENSUS_DIR $COMPARE_DIRS -type f \( -name "*.fa*" -o -name "*.fasta" \))

echo "📥 Collecting FASTA sequences..."
NUM_FASTA=$(echo "$FASTA_FILES" | wc -w)
echo "✅ Found $NUM_FASTA FASTA sequences."
echo "================================================="

# -------------------------------
# Run Nextclade using positional arguments
# -------------------------------
nextclade run \
  -D "$DATASET_DIR" \
  -O "$OUTPUT_DIR" \
  --jobs "$THREADS" \
  $FASTA_FILES

echo "✅ Nextclade analysis complete. Results in $OUTPUT_DIR"

