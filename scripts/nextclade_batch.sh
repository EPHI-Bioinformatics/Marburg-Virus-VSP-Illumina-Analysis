#!/bin/bash
set -euo pipefail

# =========================================
# Configuration
# =========================================
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONSENSUS_DIR="$PROJECT_DIR/results/07_consensus"
REFERENCE_INPUT_DIR="$PROJECT_DIR/reference_genomes/MARV_compare"
OUTPUT_DIR="$PROJECT_DIR/results/09_nextclade"
THREADS=14

# Nextclade dataset
NEXTCLADE_DATASET="$PROJECT_DIR/database/nextclade_marburg_dataset"
NEXTCLADE_ENV="nextclade_env"

# Augur configuration
AUGUR_ENV="augur_env"
AUGUR_REFERENCE="$PROJECT_DIR/database/nextclade_marburg_dataset/reference.fasta"

# =========================================
# Find all FASTA files in both consensus and reference input directories
# =========================================
mapfile -t FASTA_INPUTS < <(
    find "$CONSENSUS_DIR" "$REFERENCE_INPUT_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \)
)

if [ ${#FASTA_INPUTS[@]} -eq 0 ]; then
    echo "❌ No FASTA files found in $CONSENSUS_DIR or $REFERENCE_INPUT_DIR. Exiting."
    exit 1
fi

echo -e "\n🧬 Starting Nextclade + Augur pipeline"
echo "📁 Consensus dir       : $CONSENSUS_DIR"
echo "📁 Reference input dir: $REFERENCE_INPUT_DIR"
echo "📁 Output dir          : $OUTPUT_DIR"
echo "🧵 Threads             : $THREADS"
echo "================================================="
echo "📥 Found ${#FASTA_INPUTS[@]} FASTA sequences"

# =========================================
# Run Nextclade
# =========================================
echo -e "\n🧬 Running Nextclade..."

# Skip if Nextclade output already exists
if [ -d "$OUTPUT_DIR/nextclade" ] && [ "$(ls -A $OUTPUT_DIR/nextclade)" ]; then
    echo " Skipping Nextclade, output already exists in $OUTPUT_DIR/nextclade"
else
    # Activate Nextclade environment
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate "$NEXTCLADE_ENV"

    mkdir -p "$OUTPUT_DIR/nextclade"

    # Run Nextclade with positional FASTA arguments
    nextclade run -D "$NEXTCLADE_DATASET" -O "$OUTPUT_DIR/nextclade" "${FASTA_INPUTS[@]}"

    echo "✅ Nextclade complete. Results in $OUTPUT_DIR/nextclade"

    conda deactivate
fi

# =========================================
# Run Augur phylogenetic analysis
# =========================================
echo -e "\n🌳 Starting Augur phylogenetic analysis..."

# Skip if Augur refined tree already exists
if [ -f "$OUTPUT_DIR/augur/refined_tree.nwk" ]; then
    echo " Skipping Augur, refined tree already exists at $OUTPUT_DIR/augur/refined_tree.nwk"
else
    # Check reference file exists
    if [ ! -f "$AUGUR_REFERENCE" ]; then
        echo "❌ Augur reference file not found: $AUGUR_REFERENCE. Exiting."
        exit 1
    fi

    # Activate Augur environment
    # >>> FIX: Re-source Conda init here to avoid 'CondaError: Run 'conda init''
    source ~/miniconda3/etc/profile.d/conda.sh
    # <<< END FIX
    conda activate "$AUGUR_ENV"

    mkdir -p "$OUTPUT_DIR/augur"

    # Align sequences
    augur align --sequences "${FASTA_INPUTS[@]}" \
                --reference-sequence "$AUGUR_REFERENCE" \
                --output "$OUTPUT_DIR/augur/aligned.fasta" \
                --nthreads "$THREADS"

    # Build phylogenetic tree
    augur tree --alignment "$OUTPUT_DIR/augur/aligned.fasta" \
               --output "$OUTPUT_DIR/augur/tree.nwk"

    # Refine tree
    augur refine --tree "$OUTPUT_DIR/augur/tree.nwk" \
                 --alignment "$OUTPUT_DIR/augur/aligned.fasta" \
                 --output-tree "$OUTPUT_DIR/augur/refined_tree.nwk"

    echo "✅ Augur analysis complete. Results in $OUTPUT_DIR/augur"

    # Deactivate Augur environment
    conda deactivate
fi

echo -e "\n🎉 Pipeline finished successfully!"
