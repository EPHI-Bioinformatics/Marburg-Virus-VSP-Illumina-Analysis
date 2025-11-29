#!/bin/bash
set -euo pipefail

# ------------------------
# Config
# ------------------------
THREADS=$(nproc)
BASE_DIR="/media/betselotz/Expansion/MARG2"

CONS_DIR="$BASE_DIR/results/07_consensus"
REF_DIR="$BASE_DIR/reference_genomes/MARV_compare"

MSA_BASE="$BASE_DIR/results/10_msa"
OUTGROUP="$BASE_DIR/reference_genomes/EF446131.1.fasta"

HOSTS=("bat" "human")
MASTER_LOG="$MSA_BASE/master_summary.log"

mkdir -p "$MSA_BASE"

echo "🧠 System Threads: $(nproc) | Using $THREADS threads for MAFFT"
echo "🎯 Starting batch MAFFT pipeline..."

> "$MASTER_LOG"

# ------------------------
# Loop over hosts
# ------------------------
for HOST in "${HOSTS[@]}"; do

    echo "⚙️ Processing host: $HOST"

    HOST_DIR="$MSA_BASE/$HOST"
    mkdir -p "$HOST_DIR"

    ALIGN_FILE="$HOST_DIR/all_sequences_aligned.fasta"
    HOST_LOG="$HOST_DIR/${HOST}_summary.log"

    > "$HOST_LOG"

    # ------------------------
    # Concatenate FASTA sequences
    # ------------------------
    echo "🔗 Collecting sequences for $HOST..."

    SEQ_LIST=$(find "$REF_DIR/$HOST/cleaned/" -name "*.fasta" -size +0 | sort)
    CONS_SEQ=$(find "$CONS_DIR" -name "*$HOST*.fa" -size +0 | sort)

    TEMP_CONCAT="$HOST_DIR/all_sequences_raw.fasta"

    cat $SEQ_LIST $CONS_SEQ > "$TEMP_CONCAT"

    # ------------------------
    # Make headers unique
    # ------------------------
    awk '/^>/{print ">"$0"_"NR; next}{print}' "$TEMP_CONCAT" > "$TEMP_CONCAT.unique"
    mv "$TEMP_CONCAT.unique" "$TEMP_CONCAT"

    echo "✅ Headers made unique"

    # ------------------------
    # Additional Fix: Count sequences BEFORE MAFFT
    # ------------------------
    RAW_SEQ_COUNT=$(grep -c "^>" "$TEMP_CONCAT")
    NUM_REF_FILES=$(echo $SEQ_LIST | wc -w)
    NUM_CONS_FILES=$(echo $CONS_SEQ | wc -w)

    echo "📊 Sequence summary for $HOST:"
    echo "   Reference FASTA files: $NUM_REF_FILES"
    echo "   Consensus sequences:    $NUM_CONS_FILES"
    echo "   Total input sequences:  $RAW_SEQ_COUNT"

    # ------------------------
    # Run MAFFT (Memory-Safe Option 1: --auto)
    # ------------------------
    echo "🔥 Running MAFFT (--auto) for $HOST..."
    mafft --thread $THREADS --auto "$TEMP_CONCAT" > "$ALIGN_FILE"

    echo "✅ MAFFT alignment complete: $ALIGN_FILE"

    # ------------------------
    # Host Summary Log
    # ------------------------
    NUM_SEQ_ALIGNED=$(grep -c "^>" "$ALIGN_FILE")
    FILE_SIZE=$(stat -c%s "$ALIGN_FILE")

    {
        echo "Host: $HOST"
        echo "Reference files: $NUM_REF_FILES"
        echo "Consensus files: $NUM_CONS_FILES"
        echo "Input sequences: $RAW_SEQ_COUNT"
        echo "Aligned sequences: $NUM_SEQ_ALIGNED"
        echo "Output file size (bytes): $FILE_SIZE"
        echo "Output file: $ALIGN_FILE"
    } > "$HOST_LOG"

    echo "📝 Host summary logged: $HOST_LOG"

    echo -e "$HOST\t$NUM_SEQ_ALIGNED\t$FILE_SIZE\t$ALIGN_FILE" >> "$MASTER_LOG"

done

echo "🎉 All host-specific MAFFT alignments completed!"
echo "🗂 Master summary log: $MASTER_LOG"

