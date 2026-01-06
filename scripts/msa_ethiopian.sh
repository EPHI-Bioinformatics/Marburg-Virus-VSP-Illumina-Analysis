#!/bin/bash

# Define paths
INPUT_DIR="../results/07_consensus"
OUTPUT_DIR="../results/10_msa/Ethiopian_only"
COMBINED_FASTA="${OUTPUT_DIR}/temp_combined_ethiopian.fasta"
FINAL_MSA="${OUTPUT_DIR}/ethiopian_only_msa.fasta"

mkdir -p "$OUTPUT_DIR"

# --- Environment Activation ---
echo "Activating mafft_env..."
# Use 'source' or 'conda activate'. eval hook is more robust in scripts.
eval "$(conda shell.bash hook)"
conda activate mafft_env

echo "Step 1: Combining and cleaning headers..."
# Cleaning headers to ensure they are short and unique for clustering
cat "${INPUT_DIR}"/*.fasta "${INPUT_DIR}"/*.fa 2>/dev/null | sed 's/ .*//' > "$COMBINED_FASTA"

if [ ! -s "$COMBINED_FASTA" ]; then
    echo "Error: No consensus files found in ${INPUT_DIR}"
    conda deactivate
    exit 1
fi

echo "Step 2: Running high-accuracy MAFFT..."
# Using the recommended flags for L-INS-i accuracy:
# --localpair: Smith-Waterman pairwise alignment
# --maxiterate 1000: Iterative refinement for precision
# --reorder: Groups similar sequences together
mafft --localpair --maxiterate 1000 --reorder --thread -1 "$COMBINED_FASTA" > "$FINAL_MSA"

# --- Environment Deactivation ---
echo "Deactivating mafft_env..."
conda deactivate

# Clean up
rm "$COMBINED_FASTA"

echo "Process complete. Alignment saved to: ${FINAL_MSA}"
