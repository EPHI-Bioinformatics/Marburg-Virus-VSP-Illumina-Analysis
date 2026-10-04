#!/bin/bash
set -euo pipefail

# Suppress SSL warnings
export PYTHONWARNINGS="ignore:Unverified HTTPS request"

# -------------------------------------------------
# PATHS
# -------------------------------------------------
BASE_DIR="/media/betselotz/Expansion/MARV-Gen"
MSA_DIR="$BASE_DIR/results/10_msa/MARV.A.1"
TREE_DIR="$BASE_DIR/results/11_phylogeny/MARV.A.1"

ALIGNED_MAX="$MSA_DIR/marburg_aligned.fasta"
SNPS_ONLY="$MSA_DIR/marburg_snps_only.fasta"
TREE_PREFIX="$TREE_DIR/marburg_ml"

# -------------------------------------------------
# SETUP
# -------------------------------------------------
mkdir -p "$TREE_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate iqtree_env

# -------------------------------------------------
# STEP 1: EXTRACT SNP SITES (Strict Mode)
# -------------------------------------------------
if [ ! -f "$SNPS_ONLY" ]; then
    echo ">>> Extracting polymorphic sites (SNPs)..."
    snp-sites -m -c -o "$SNPS_ONLY" "$ALIGNED_MAX"
fi

# -------------------------------------------------
# STEP 2: RUN IQ-TREE
# -------------------------------------------------
echo ">>> Running IQ-TREE Analysis..."
iqtree -s "$SNPS_ONLY" \
       -m MFP \
       -o JN408064.1 \
       -B 1000 \
       -T AUTO \
       --prefix "$TREE_PREFIX"

# -------------------------------------------------
# STEP 3: FORMAT OUTPUTS (Cleaned & Newick)
# -------------------------------------------------
TREE_FILE="$TREE_PREFIX.treefile"
TREE_CLEAN="$TREE_DIR/marburg_ml_clean.treefile"
NEWICK_FINAL="$TREE_DIR/marburg_final.newick"

if [ -f "$TREE_FILE" ]; then
    echo ">>> Generating final Newick and cleaned files..."
    
    # 1. Create Cleaned version
    sed -E 's/(_)+(:)/\2/g' "$TREE_FILE" > "$TREE_CLEAN"
    
    # 2. Copy/Rename to .newick for downstream software
    cp "$TREE_CLEAN" "$NEWICK_FINAL"
    
    echo ">>> Cleaned tree: $TREE_CLEAN"
    echo ">>> Newick file:  $NEWICK_FINAL"
fi

conda deactivate

echo "------------------------------------------------"
echo "SNP-BASED PHYLOGENY PIPELINE COMPLETE"
echo "Final Newick File: $NEWICK_FINAL"
echo "------------------------------------------------"
