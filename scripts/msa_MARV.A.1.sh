#!/bin/bash
set -euo pipefail

# -------------------------------------------------
# CONFIGURATION
# -------------------------------------------------

BASE_DIR="/media/betselotz/Expansion/MARV-Gen"
REF_DIR="$BASE_DIR/reference_genomes/MARV_compare"
CONSENSUS_DIR="$BASE_DIR/results/07_consensus"
OUTPUT_DIR="$BASE_DIR/results/10_msa/MARV.A.1"

COMBINED="$OUTPUT_DIR/combined_marburg.fasta"
ALIGNED_FINAL="$OUTPUT_DIR/marburg_aligned.fasta"

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
# STEP 1: COMBINE FASTA FILES FOR SELECTED SAMPLES ONLY
# -------------------------------------------------

SAMPLES=(
    "JN408064.1"
    "JX458851.1"
    "KC545387.1"
    "KC545388.1"
    "JX458853.1"
    "JX458858.1"
    "EF446132.1"
    "ET_MARV_23"
    "ET_MARV_31"
    "ET_MARV_32"
    "ET_MARV_261"
    "ET_MARV_262"
    "ET_MARV_45"
    "ET_MARV_60"
)

if [ -s "$COMBINED" ]; then
    echo ">>> Skipping FASTA combination (already exists): $COMBINED"
else
    echo ">>> Combining selected samples from reference and consensus FASTA files..."
    
    > "$COMBINED"  # initialize empty file

    for SAMPLE in "${SAMPLES[@]}"; do
        # search in reference directory
        REF_FILE=$(find "$REF_DIR" -type f \( -name "${SAMPLE}.fasta" -o -name "${SAMPLE}.fa" -o -name "${SAMPLE}.fna" \) | head -n1)
        # search in consensus directory if not found
        if [ -z "$REF_FILE" ]; then
            REF_FILE=$(find "$CONSENSUS_DIR" -type f \( -name "${SAMPLE}.fasta" -o -name "${SAMPLE}.fa" -o -name "${SAMPLE}.fna" \) | head -n1)
        fi
        if [ -n "$REF_FILE" ]; then
            cat "$REF_FILE" >> "$COMBINED"
            echo "" >> "$COMBINED"  # ensure newline
        else
            echo "⚠️ WARNING: Sample not found: $SAMPLE"
        fi
    done

    if [ ! -s "$COMBINED" ]; then
        echo "❌ ERROR: No FASTA files found for selected samples."
        exit 1
    fi
fi

# -------------------------------------------------
# STEP 2: MAFFT L-INS-i ALIGNMENT
# -------------------------------------------------

if [ -s "$ALIGNED_FINAL" ]; then
    echo ">>> Skipping MAFFT alignment (already exists): $ALIGNED_FINAL"
else
    echo ">>> Running MAFFT L-INS-i (Most Accurate Strategy) on selected samples..."

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
