#!/bin/bash
set -euo pipefail

###############################################################################
# Qualimap BAM QC Batch Script
# Description: Perform BAM QC for all BAM files in mapping directory
# Automatically activates conda environment, uses all available CPUs, and logs output
###############################################################################

# -------- CPU AUTO-DETECTION --------
TOTAL_THREADS=$(nproc)
SAFE_THREADS=$(( TOTAL_THREADS > 1 ? TOTAL_THREADS - 1 : 1 ))
echo "🧠 Using $SAFE_THREADS threads for Qualimap"

# -------- Project Directories --------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MAP_DIR="$PROJECT_DIR/results/04_mapped_bam"
QUALIMAP_DIR="$PROJECT_DIR/results/05_mapping_qc"
LOG_DIR="$QUALIMAP_DIR/logs"

mkdir -p "$QUALIMAP_DIR" "$LOG_DIR"

# -------- Activate Qualimap Conda Environment --------
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qualimap_env
set -u

# -------- Run Qualimap on each BAM with logging --------
for bam in "$MAP_DIR"/*.bam; do
    SAMPLE_NAME=$(basename "$bam" .bam)
    OUT_DIR="$QUALIMAP_DIR/$SAMPLE_NAME"
    mkdir -p "$OUT_DIR"
    
    LOG_FILE="$LOG_DIR/${SAMPLE_NAME}_qualimap.log"
    
    echo -e "🔹 Running Qualimap on $SAMPLE_NAME..."
    echo -e "   Log file: $LOG_FILE"

    # Redirect stdout and stderr to log file
    qualimap bamqc -bam "$bam" -outdir "$OUT_DIR" -nt "$SAFE_THREADS" \
        > "$LOG_FILE" 2>&1
    
    echo -e "✅ Completed Qualimap for $SAMPLE_NAME"
done

echo -e "🎉 All BAM QC completed! Results in $QUALIMAP_DIR"

