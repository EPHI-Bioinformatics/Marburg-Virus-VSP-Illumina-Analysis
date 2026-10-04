#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

###############################################################################
# 🧬 Host read removal using HOSTILE v2.0.2
###############################################################################

# -------------------------------
# Activate hostile environment
# -------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hostile_env

# -------------------------------
# Paths and threads
# -------------------------------
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$BASE_DIR/results/02_clean_reads"
OUT_DIR="$BASE_DIR/results/03_nonhuman_reads"
STATS_DIR="$OUT_DIR/alignment_stats"
HOST_INDEX="$BASE_DIR/reference_genomes/human_index"
THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))

mkdir -p "$OUT_DIR" "$STATS_DIR"

# -------------------------------
# Detect all R1 files safely
# -------------------------------
shopt -s nullglob
r1_files=("$RAW_DIR"/*_R1.trimmed.fastq.gz)
shopt -u nullglob

if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "❌ No R1 FASTQ files found in $RAW_DIR. Exiting."
    conda deactivate
    exit 1
fi

# -------------------------------
# CSV header
# -------------------------------
CSV_FILE="$STATS_DIR/host_removal_stats.csv"
echo "Sample,Original_Reads_PE,Nonhost_Reads_PE,Percent_Host_Removed" > "$CSV_FILE"

# -------------------------------
# Loop through samples
# -------------------------------
for R1 in "${r1_files[@]}"; do
    SAMPLE=$(basename "$R1" _R1.trimmed.fastq.gz)
    R2="$RAW_DIR/${SAMPLE}_R2.trimmed.fastq.gz"
    OUT1="$OUT_DIR/${SAMPLE}_nonhost_R1.fastq.gz"
    OUT2="$OUT_DIR/${SAMPLE}_nonhost_R2.fastq.gz"

    echo "============================================"
    echo "Processing sample: $SAMPLE"
    echo "============================================"

    # Skip if already processed
    if [[ -f "$OUT1" && -f "$OUT2" ]]; then
        echo "Skipping $SAMPLE (already has nonhost FASTQ)"
        continue
    fi

    # Check paired-end files exist
    if [[ ! -f "$R2" ]]; then
        echo "❌ Missing R2 for $SAMPLE. Skipping."
        continue
    fi

    # -------------------------------
    # Run hostile clean
    # -------------------------------
    # Note: Hostile 2.0+ appends '.clean_1.fastq.gz' to the input filename
    hostile clean \
        --index "$HOST_INDEX" \
        --fastq1 "$R1" \
        --fastq2 "$R2" \
        --threads "$THREADS" \
        -o "$OUT_DIR" \
        --force

    # -------------------------------
    # Handle Hostile's specific naming convention
    # -------------------------------
    # Input:  SAMPLE_R1.trimmed.fastq.gz
    # Output: SAMPLE_R1.trimmed.clean_1.fastq.gz
    TEMP1="$OUT_DIR/$(basename "$R1" .fastq.gz).clean_1.fastq.gz"
    TEMP2="$OUT_DIR/$(basename "$R2" .fastq.gz).clean_2.fastq.gz"

    if [[ -f "$TEMP1" && -f "$TEMP2" ]]; then
        mv "$TEMP1" "$OUT1"
        mv "$TEMP2" "$OUT2"
        echo "✅ $SAMPLE: renamed to nonhost FASTQs"
    else
        echo "❌ $SAMPLE: Expected files ($TEMP1) not found!"
        continue
    fi

    # -------------------------------
    # Calculate read counts for stats
    # -------------------------------
    # Using 'zgrep -c' on the '@' symbol header is faster than 'wc -l'
    ORIG_READS=$(zgrep -c "^@" "$R1")
    NONHOST_READS=$(zgrep -c "^@" "$OUT1")
    
    # Calculate percentage (using awk for floating point)
    PCT_REMOVED=$(awk "BEGIN {print ($ORIG_READS - $NONHOST_READS) / $ORIG_READS * 100}")

    echo "$SAMPLE,$ORIG_READS,$NONHOST_READS,${PCT_REMOVED}%" >> "$CSV_FILE"
    
    echo "📊 Stats: Original=$ORIG_READS | Clean=$NONHOST_READS | Removed=${PCT_REMOVED}%"
done

# -------------------------------
# Deactivate environment
# -------------------------------
conda deactivate
echo "============================================"
echo "✅ Host read removal completed."
echo "📝 Stats saved to: $CSV_FILE"
