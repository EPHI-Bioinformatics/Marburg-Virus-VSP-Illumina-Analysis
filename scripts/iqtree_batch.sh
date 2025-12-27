#!/bin/bash

# --- 1. SET PATHS ---
REF_DIR="../reference_genomes/MARV_compare"
CONSENSUS_DIR="../results/07_consensus"
MSA_DIR="../results/10_msa"
TREE_DIR="../results/11_phylogeny"
METADATA="../metadata/all_seq_metadata.csv"

# Filenames
COMBINED="$MSA_DIR/combined_marburg.fasta"
ALIGNED_MAX="$MSA_DIR/marburg_aligned.fasta"
ALIGNED_AUTO="$MSA_DIR/marburg_auto_aligned.fasta"
TREE_FILE="$TREE_DIR/marburg_ml.treefile"
TREE_CLEAN="$TREE_DIR/marburg_ml_clean.treefile"

# Create directories
mkdir -p "$MSA_DIR"
mkdir -p "$TREE_DIR"

# Initialize Conda
source "$(conda info --base)/etc/profile.d/conda.sh"

# --- 2. RUN MAFFT (MAX QUALITY) ---
if [ -f "$ALIGNED_MAX" ]; then
    echo ">>> Skipping MAFFT: Alignment file '$ALIGNED_MAX' already exists."
else
    echo ">>> Activating mafft_env..."
    conda activate mafft_env

    echo ">>> Combining FASTA files..."
    cat "$REF_DIR"/*.fasta "$CONSENSUS_DIR"/*.fa > "$COMBINED"

    echo ">>> Running High-Quality MSA (L-INS-i)..."
    mafft --localpair --maxiterate 1000 --thread -1 "$COMBINED" > "$ALIGNED_MAX"

    echo ">>> Running Auto MSA for comparison..."
    mafft --auto --thread -1 "$COMBINED" > "$ALIGNED_AUTO"
    
    conda deactivate
fi

# --- 3. RUN IQ-TREE ---
if [ -f "$TREE_FILE" ]; then
    echo ">>> Skipping IQ-TREE: Tree file '$TREE_FILE' already exists."
else
    echo ">>> Activating iqtree_env..."
    conda activate iqtree_env

    echo ">>> Running IQ-TREE Analysis..."
    iqtree -s "$ALIGNED_MAX" \
           -m MFP \
           -o DQ447649.1 \
           -B 1000 \
           -T AUTO \
           --prefix "$TREE_DIR/marburg_ml"
    
    conda deactivate
fi

# --- 4. CLEAN TREE TIP NAMES TO MATCH METADATA ---
if [ -f "$TREE_CLEAN" ]; then
    echo ">>> Skipping tree cleaning: '$TREE_CLEAN' already exists."
else
    echo ">>> Cleaning tree tip names to match metadata..."
    # Remove trailing underscores before branch lengths (common cause of TreeTime mismatches)
    sed -E 's/(_)+(:)/\2/g' "$TREE_FILE" > "$TREE_CLEAN"
    echo ">>> Cleaned tree saved as: $TREE_CLEAN"
fi

echo "------------------------------------------------"
echo "PROCESS CHECK COMPLETE"
echo "Alignment:      $ALIGNED_MAX"
echo "Tree File:      $TREE_FILE"
echo "Clean Tree:     $TREE_CLEAN"
echo "------------------------------------------------"

