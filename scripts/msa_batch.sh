#!/bin/bash
set -euo pipefail

# -------------------------------------------------
# CONFIGURATION
# -------------------------------------------------

# Using the absolute path to ensure accuracy regardless of execution directory
BASE_DIR="/media/betselotz/Expansion/MARV-Gen"
REF_DIR="$BASE_DIR/reference_genomes/MARV_compare"
CONSENSUS_DIR="$BASE_DIR/results/07_consensus"
OUTPUT_DIR="$BASE_DIR/results/10_msa"

COMBINED="$OUTPUT_DIR/combined_marburg.fasta"
ALIGNED_FINAL="$OUTPUT_DIR/marburg_aligned.fasta"

# MAFFT L-INS-i is highly accurate but memory-intensive.
# Capping threads prevents memory crashes on high-core systems.
MAX_MAFFT_THREADS=6

# -------------------------------------------------
# SETUP
# -------------------------------------------------

mkdir -p "$OUTPUT_DIR"

CPU_DETECTED=$(nproc)
THREADS=$(( CPU_DETECTED > MAX_MAFFT_THREADS ? MAX_MAFFT_THREADS : CPU_DETECTED ))

echo ">>> Detected CPU cores : $CPU_DETECTED"
echo ">>> MAFFT threads used : $THREADS (capped for stability)"

# Initialize Conda
CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate mafft_env

# -------------------------------------------------
# STEP 1: COMBINE FASTA FILES
# -------------------------------------------------

if [ -s "$COMBINED" ]; then
    echo ">>> Skipping FASTA combination (already exists): $COMBINED"
else
    echo ">>> Combining reference and consensus FASTA files..."
    
    # nullglob prevents errors if one extension (e.g., .fa) isn't found
    shopt -s nullglob
    cat "$REF_DIR"/*.{fasta,fa,fna} "$CONSENSUS_DIR"/*.{fasta,fa,fna} > "$COMBINED"
    shopt -u nullglob

    if [ ! -s "$COMBINED" ]; then
        echo "❌ ERROR: No FASTA files found in $REF_DIR or $CONSENSUS_DIR"
        exit 1
    fi
fi

# -------------------------------------------------
# STEP 2: MAFFT L-INS-i ALIGNMENT (Most Accurate)
# -------------------------------------------------

if [ -s "$ALIGNED_FINAL" ]; then
    echo ">>> Skipping MAFFT alignment (already exists): $ALIGNED_FINAL"
else
    echo ">>> Running MAFFT L-INS-i (Most Accurate Strategy)..."

    # --localpair + --maxiterate 1000 activates the L-INS-i algorithm
    mafft \
      --localpair \
      --maxiterate 1000 \
      --thread "$THREADS" \
      "$COMBINED" > "$ALIGNED_FINAL" || {
        echo "❌ ERROR: MAFFT L-INS-i failed."
        exit 1
      }

    if [ ! -s "$ALIGNED_FINAL" ]; then
        echo "❌ ERROR: Alignment output is empty!"
        exit 1
    fi
fi

conda deactivate

echo "------------------------------------------------"
echo "MSA PIPELINE COMPLETE"
echo "Final alignment: $ALIGNED_FINAL"
echo "------------------------------------------------"
