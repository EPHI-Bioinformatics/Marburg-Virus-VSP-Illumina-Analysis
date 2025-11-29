#!/bin/bash
set -euo pipefail

# ========================================
# IQ-TREE Phylogeny Pipeline (Improved)
# ========================================

# Reserve 2 cores for system, use rest
THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
echo "🧠 Using $THREADS threads for IQ-TREE"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MSA_DIR="$PROJECT_DIR/results/10_msa"
PHYLO_DIR="$PROJECT_DIR/results/11_phylogeny"

HOSTS=("human" "bat")
MASTER_LOG="$PHYLO_DIR/master_iqtree_summary.log"

mkdir -p "$PHYLO_DIR"
> "$MASTER_LOG"   # clear master log


# ========================================
# RUN IQTREE FUNCTION
# ========================================
run_iqtree() {

    local host=$1
    local msa_file="$MSA_DIR/$host/all_sequences_aligned.fasta"

    local phylo_host_dir="$PHYLO_DIR/$host"
    local inter_dir="$phylo_host_dir/intermediates"
    local log_file="$phylo_host_dir/iqtree_${host}_$(date +%Y%m%d_%H%M%S).log"

    mkdir -p "$phylo_host_dir" "$inter_dir"

    echo "----------------------------------------"
    echo "⚙️  Running IQ-TREE for host: $host"
    echo "📂 Input MSA: $msa_file"
    echo "📁 Output directory: $phylo_host_dir"
    echo "----------------------------------------"


    # =====================================================
    # 1. Validate MSA file exists
    # =====================================================
    if [[ ! -s "$msa_file" ]]; then
        echo "❌ ERROR: MSA file missing or empty for $host: $msa_file"
        echo -e "$host\tERROR: missing MSA" >> "$MASTER_LOG"
        return 1
    fi


    # =====================================================
    # 2. Count sequences BEFORE tree building
    # =====================================================
    local seq_count
    seq_count=$(grep -c "^>" "$msa_file")

    echo "📊 Sequence count for $host: $seq_count"

    if (( seq_count < 3 )); then
        echo "❌ ERROR: Not enough sequences for phylogeny (minimum is 3)."
        echo -e "$host\tERROR: insufficient sequences" >> "$MASTER_LOG"
        return 1
    fi


    # =====================================================
    # 3. Ensure IQ-TREE is installed
    # =====================================================
    IQTREE_EXEC=$(which iqtree3 || true)

    if [[ -z "$IQTREE_EXEC" ]]; then
        echo "❌ ERROR: iqtree3 not found in PATH"
        echo -e "$host\tERROR: iqtree3 missing" >> "$MASTER_LOG"
        return 1
    fi


    # =====================================================
    # 4. Run IQ-TREE
    # =====================================================
    echo "🚀 Starting IQ-TREE (ModelFinder + 1000 ultrafast bootstraps)..."

    "$IQTREE_EXEC" \
        -s "$msa_file" \
        -m MFP \
        -bb 1000 \
        -nt "$THREADS" \
        -pre "$phylo_host_dir/${host}" \
        2>&1 | tee "$log_file"

    echo "✅ IQ-TREE finished for $host"


    # =====================================================
    # 5. Validate tree output
    # =====================================================
    local tree_file="$phylo_host_dir/${host}.treefile"
    if [[ ! -s "$tree_file" ]]; then
        echo "❌ ERROR: Treefile missing for $host"
        echo -e "$host\tERROR: no treefile" >> "$MASTER_LOG"
        return 1
    fi


    # =====================================================
    # 6. Move intermediate files cleanly
    # =====================================================
    echo "📦 Organizing intermediate files..."

    for f in "$phylo_host_dir"/*; do
        [[ "$f" == "$tree_file" ]] && continue
        [[ "$f" == "$log_file" ]] && continue
        mv "$f" "$inter_dir/" 2>/dev/null || true
    done


    # =====================================================
    # 7. Write Host Summary
    # =====================================================
    local file_size
    file_size=$(stat -c%s "$tree_file")

    {
        echo "Host: $host"
        echo "Sequences: $seq_count"
        echo "Tree file: $tree_file"
        echo "Tree size (bytes): $file_size"
        echo "Log file: $log_file"
        echo "Date: $(date)"
    } > "$phylo_host_dir/${host}_summary.log"

    echo -e "$host\t$tree_file\t$file_size\t$log_file\t$seq_count" >> "$MASTER_LOG"

    echo "📝 Summary written for $host"
    echo "----------------------------------------"

}


# ========================================
# RUN ALL HOSTS IN PARALLEL
# ========================================
for host in "${HOSTS[@]}"; do
    run_iqtree "$host" &
done
wait


echo "🎉 All host-specific IQ-TREE phylogenies completed!"
echo "🗂 Master summary log: $MASTER_LOG"

