#!/bin/bash
set -euo pipefail

# -------------------------------------------------
# PATHS
# -------------------------------------------------
BASE_DIR="/media/betselotz/Expansion/MARV-Gen"

MSA_DIR="$BASE_DIR/results/10_msa/MARV.A.1"
TREE_DIR="$BASE_DIR/results/11_phylogeny/MARV.A.1"
METADATA="$BASE_DIR/metadata/all_seq_metadata.csv"

ALIGNED_MAX="$MSA_DIR/marburg_aligned.fasta"
ALIGNED_AUTO="$MSA_DIR/marburg_auto_aligned.fasta"
TREE_FILE="$TREE_DIR/marburg_ml.treefile"
TREE_CLEAN="$TREE_DIR/marburg_ml_clean.treefile"

# -------------------------------------------------
# SETUP
# -------------------------------------------------
mkdir -p "$TREE_DIR"

# Initialize Conda
source "$(conda info --base)/etc/profile.d/conda.sh"

# -------------------------------------------------
# STEP 1: RUN AUTO MSA (Optional)
# -------------------------------------------------
if [ ! -f "$ALIGNED_AUTO" ]; then
    echo ">>> Activating mafft_env..."
    conda activate mafft_env

    echo ">>> Running Auto MSA on existing alignment..."
    mafft --auto --thread -1 "$ALIGNED_MAX" > "$ALIGNED_AUTO"

    conda deactivate
else
    echo ">>> Skipping Auto MSA (already exists): $ALIGNED_AUTO"
fi

# -------------------------------------------------
# STEP 2: RUN IQ-TREE WITH OUTGROUP
# -------------------------------------------------
if [ -f "$TREE_FILE" ]; then
    echo ">>> Skipping IQ-TREE: Tree file '$TREE_FILE' already exists."
else
    echo ">>> Activating iqtree_env..."
    conda activate iqtree_env

    echo ">>> Running IQ-TREE Analysis (JN408064.1 as outgroup)..."
    iqtree -s "$ALIGNED_MAX" \
           -m MFP \
           -o JN408064.1 \
           -B 1000 \
           -T AUTO \
           --prefix "$TREE_DIR/marburg_ml"

    conda deactivate
fi

# -------------------------------------------------
# STEP 3: CLEAN TREE TIP NAMES TO MATCH METADATA
# -------------------------------------------------
if [ -f "$TREE_CLEAN" ]; then
    echo ">>> Skipping tree cleaning: '$TREE_CLEAN' already exists."
else
    echo ">>> Cleaning tree tip names to match metadata..."
    sed -E 's/(_)+(:)/\2/g' "$TREE_FILE" > "$TREE_CLEAN"
    echo ">>> Cleaned tree saved as: $TREE_CLEAN"
fi

echo "------------------------------------------------"
echo "PHYLOGENY PIPELINE COMPLETE"
echo "Alignment:      $ALIGNED_MAX"
echo "Auto Alignment: $ALIGNED_AUTO"
echo "Tree File:      $TREE_FILE"
echo "Clean Tree:     $TREE_CLEAN"
echo "------------------------------------------------"
