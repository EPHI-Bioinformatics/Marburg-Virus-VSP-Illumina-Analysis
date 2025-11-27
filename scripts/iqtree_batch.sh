#!/bin/bash
set -euo pipefail

# ========================================
# IQ-TREE Phylogeny Pipeline with Original Leaf Names
# ========================================

export THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
IQ_PREFIX="marv_phylogeny"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MSA_DIR="$PROJECT_DIR/results/09_msa2"
PHYLO_DIR="$PROJECT_DIR/results/10_phylogeny"
INTERMEDIATES_DIR="$PHYLO_DIR/intermediates"

ALN_TRIM="$MSA_DIR/all_sequences_aligned_renamed.fasta"
CLEANED_FASTA="$PHYLO_DIR/${IQ_PREFIX}_cleaned.fasta"
MAPPING_FILE="$PHYLO_DIR/${IQ_PREFIX}_header_mapping.tsv"
TREE_FILE="$PHYLO_DIR/${IQ_PREFIX}.treefile"
SUMMARY_FILE="$PHYLO_DIR/iqtree_summary.txt"
LOG_FILE="$PHYLO_DIR/${IQ_PREFIX}_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$PHYLO_DIR" "$INTERMEDIATES_DIR"

echo "⚙️ Running IQ-TREE (v3.0.1) with ${THREADS} threads." | tee "$LOG_FILE"

# --- Cleanup trap ---
cleanup_on_error() {
    echo "❌ IQ-TREE failed. Review $LOG_FILE" | tee -a "$LOG_FILE" >&2
    rm -f "$CLEANED_FASTA" 2>/dev/null || true
}
trap cleanup_on_error ERR

# --- 1. Check Input File ---
if [[ ! -s "$ALN_TRIM" ]]; then
    echo "❌ Trimmed MSA not found: $ALN_TRIM" | tee -a "$LOG_FILE"
    exit 1
fi

# --- 2. Clean FASTA headers ---
echo "🔹 Cleaning FASTA headers for IQ-TREE..." | tee -a "$LOG_FILE"
> "$CLEANED_FASTA"
> "$MAPPING_FILE"

while read -r line; do
    if [[ "$line" =~ ^">" ]]; then
        orig_header="${line#>}"
        cleaned_header=$(echo "$orig_header" | sed 's/ | /|/g; s/ /_/g')
        echo ">$cleaned_header" >> "$CLEANED_FASTA"
        echo -e "$cleaned_header\t$orig_header" >> "$MAPPING_FILE"
    else
        echo "$line" >> "$CLEANED_FASTA"
    fi
done < "$ALN_TRIM"

echo "✅ FASTA headers cleaned. Mapping file: $MAPPING_FILE" | tee -a "$LOG_FILE"

# --- 3. Skip if tree exists ---
if [[ -s "$TREE_FILE" ]]; then
    echo "✅ Treefile exists: $TREE_FILE" | tee -a "$LOG_FILE"
    exit 0
fi

# --- 4. Run IQ-TREE ---
IQTREE_EXEC=$(which iqtree3)
if [[ -z "$IQTREE_EXEC" ]]; then
    echo "❌ iqtree3 not found in PATH" | tee -a "$LOG_FILE"
    exit 1
fi

echo "🔹 Running IQ-TREE with cleaned FASTA..." | tee -a "$LOG_FILE"
START=$(date +%s)

"$IQTREE_EXEC" -s "$CLEANED_FASTA" -m MFP -bb 1000 -nt "$THREADS" 2>&1 | tee -a "$LOG_FILE"

END=$(date +%s)
RUNTIME=$((END-START))

# --- 5. Replace leaf names in treefile ---
IQTREE_RAW_TREE="$PHYLO_DIR/${IQ_PREFIX}.treefile"
IQTREE_ORIG_TREE="$PHYLO_DIR/${IQ_PREFIX}_orig_leafnames.treefile"

if [[ -s "$IQTREE_RAW_TREE" ]]; then
    echo "🔹 Replacing leaf names with original headers..." | tee -a "$LOG_FILE"

    # Build sed replacement commands from mapping file
    SED_SCRIPT=$(awk '{gsub(/[\[\]\/&]/,"\\&",$1); gsub(/[\[\]\/&]/,"\\&",$2); print "s/"$1"/"$2"/g"}' "$MAPPING_FILE")
    
    sed "$SED_SCRIPT" "$IQTREE_RAW_TREE" > "$IQTREE_ORIG_TREE"

    echo "✅ Tree with original leaf names saved to: $IQTREE_ORIG_TREE" | tee -a "$LOG_FILE"
else
    echo "❌ ERROR: IQ-TREE treefile not found!" | tee -a "$LOG_FILE"
    exit 1
fi

# --- 6. Move intermediates ---
find "$PHYLO_DIR" -maxdepth 1 -type f -name "$IQ_PREFIX.*" \
    ! -name "${IQ_PREFIX}.treefile" \
    ! -name "${IQ_PREFIX}_orig_leafnames.treefile" \
    -exec mv {} "$INTERMEDIATES_DIR/" \; 2>/dev/null || true

# --- 7. Summary ---
if [[ ! -f "$SUMMARY_FILE" ]]; then
    echo -e "Treefile\tRuntime_seconds\tBest_Model" > "$SUMMARY_FILE"
fi

BEST_MODEL=$(grep "ModelFinder best model:" "$INTERMEDIATES_DIR/$IQ_PREFIX.log" 2>/dev/null | awk '{print $NF}' || echo "N/A")
echo -e "$IQTREE_ORIG_TREE\t$RUNTIME\t$BEST_MODEL" >> "$SUMMARY_FILE"
echo "✅ Summary appended: $SUMMARY_FILE" | tee -a "$LOG_FILE"

