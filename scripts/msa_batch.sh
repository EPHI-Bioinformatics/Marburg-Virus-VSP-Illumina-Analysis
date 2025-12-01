#!/bin/bash
set -euo pipefail

# ========================================
# Fast MAFFT Batch Alignment Pipeline
# ========================================

# ------------------------
# Config
# ------------------------
THREADS=$(nproc)  # use all available CPU cores
BASE_DIR="/media/betselotz/Expansion/MARG2"
CONS_DIR="$BASE_DIR/results/07_consensus"
REF_DIR="$BASE_DIR/reference_genomes/MARV_compare"
MSA_BASE="$BASE_DIR/results/11_msa"
OUTGROUP="$BASE_DIR/reference_genomes/EF446131.1.fasta"

# Hosts to process
HOSTS=("bat" "human")

# Master summary log
MASTER_LOG="$MSA_BASE/master_summary.log"
mkdir -p "$MSA_BASE"
> "$MASTER_LOG"

echo "🧠 System Threads: $(nproc) | Using $THREADS threads for MAFFT"
echo "🎯 Starting fast batch MAFFT pipeline..."

# ------------------------
# Loop over hosts
# ------------------------
for HOST in "${HOSTS[@]}"; do
    echo "----------------------------------------"
    echo "⚙️ Processing host: $HOST"

    HOST_DIR="$MSA_BASE/$HOST"
    mkdir -p "$HOST_DIR"

    ALIGN_FILE="$HOST_DIR/all_sequences_aligned.fasta"
    HOST_LOG="$HOST_DIR/${HOST}_summary.log"
    > "$HOST_LOG"

    # ------------------------
    # Gather all sequences
    # ------------------------
    SEQ_LIST=$(find "$REF_DIR/$HOST/cleaned/" -name "*.fasta" -size +0 | sort)

    # Include all consensus sequences (fix for missing consensus)
    CONS_SEQ=$(find "$CONS_DIR" -name "*.fa" -size +0 | sort)

    ALL_INPUT="$HOST_DIR/all_sequences_raw.fasta"

    # Concatenate sequences + outgroup
    cat $SEQ_LIST $CONS_SEQ "$OUTGROUP" > "$ALL_INPUT"

    # Make headers unique
    awk '/^>/{print ">"$0"_"NR; next}{print}' "$ALL_INPUT" > "$ALL_INPUT.unique"
    mv "$ALL_INPUT.unique" "$ALL_INPUT"

    # ------------------------
    # Dry-run: show what will be aligned
    # ------------------------
    echo "📋 Dry-run: sequences to be aligned for $HOST:"
    grep "^>" "$ALL_INPUT" | sed 's/>//'
    TOTAL_SEQ=$(grep -c "^>" "$ALL_INPUT")
    echo "✅ Total sequences for $HOST: $TOTAL_SEQ"

    # ------------------------
    # Run MAFFT (fast & safe)
    # ------------------------
    echo "🔥 Running MAFFT --auto for $HOST..."
    mafft --thread "$THREADS" --auto "$ALL_INPUT" > "$ALIGN_FILE"
    echo "✅ MAFFT alignment complete: $ALIGN_FILE"

    # ------------------------
    # Host summary log
    # ------------------------
    NUM_SEQ_ALIGNED=$(grep -c "^>" "$ALIGN_FILE")
    FILE_SIZE=$(stat -c%s "$ALIGN_FILE")
    {
        echo "Host: $HOST"
        echo "Reference files: $(echo $SEQ_LIST | wc -w)"
        echo "Consensus files: $(echo $CONS_SEQ | wc -w)"
        echo "Outgroup: $(basename $OUTGROUP)"
        echo "Input sequences: $TOTAL_SEQ"
        echo "Aligned sequences: $NUM_SEQ_ALIGNED"
        echo "Output file size (bytes): $FILE_SIZE"
        echo "Output file: $ALIGN_FILE"
        echo "Date: $(date)"
    } > "$HOST_LOG"
    echo "📝 Host summary logged: $HOST_LOG"

    # Append to master summary
    echo -e "$HOST\t$NUM_SEQ_ALIGNED\t$FILE_SIZE\t$ALIGN_FILE" >> "$MASTER_LOG"

done

echo "🎉 All host-specific MAFFT alignments completed!"
echo "🗂 Master summary log: $MASTER_LOG"

