#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_DIR="$PROJECT_DIR/results/04_mapped_bam"
CONS_DIR="$PROJECT_DIR/results/06_consensus"
LOG_FILE="$PROJECT_DIR/pipeline_consensus_log_$(date +%Y%m%d_%H%M%S).txt"

mkdir -p "$CONS_DIR"
export THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))

echo "Starting consensus generation." | tee -a "$LOG_FILE"
echo "Using $THREADS threads (samtools + ivar)" | tee -a "$LOG_FILE"
echo "--------------------------------------------------------" | tee -a "$LOG_FILE"

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ivar_env
set -u

SAMTOOLS_EXEC=$(which samtools)
IVAR_EXEC=$(which ivar)

BAM_FILES=("$BAM_DIR"/*.sorted.bam)
if [[ ! -e "${BAM_FILES[0]}" ]]; then
    echo "FATAL: No *.sorted.bam files in $BAM_DIR" | tee -a "$LOG_FILE"
    conda deactivate
    exit 1
fi

for bam in "${BAM_FILES[@]}"; do
    START_TIME=$(date +%s)
    sample=$(basename "$bam" .sorted.bam)
    consensus_file="$CONS_DIR/${sample}.fa"

    if [[ -s "$consensus_file" ]]; then
        echo "Existing consensus for $sample, skipping." | tee -a "$LOG_FILE"
        continue
    fi

    echo "Processing $sample" | tee -a "$LOG_FILE"

    if ! (
        "$SAMTOOLS_EXEC" mpileup -A -Q 0 "$bam" |
        "$IVAR_EXEC" consensus -p "$CONS_DIR/$sample" -q 20 -t 0.7 -m 1
    ) 2>> "$LOG_FILE"; then
        echo "ERROR: Consensus failed for $sample" | tee -a "$LOG_FILE"
        rm -f "$consensus_file"
        continue
    fi

    if [[ -f "$CONS_DIR/${sample}.fa" ]]; then
        sed -i "1s/.*/>${sample} | Homo sapiens | Ethiopia | 2025/" "$CONS_DIR/${sample}.fa"
        echo "Consensus saved: $CONS_DIR/${sample}.fa" | tee -a "$LOG_FILE"
    else
        echo "WARNING: Consensus file missing for $sample" | tee -a "$LOG_FILE"
    fi

    END_TIME=$(date +%s)
    echo "Runtime for $sample: $((END_TIME - START_TIME))s" | tee -a "$LOG_FILE"
done

conda deactivate
echo "Environment deactivated." | tee -a "$LOG_FILE"

TOTAL=$(find "$CONS_DIR" -maxdepth 1 -name "*.fa" -type f | wc -l)
echo "Consensus generation complete." | tee -a "$LOG_FILE"
echo "Total FASTA files: $TOTAL" | tee -a "$LOG_FILE"
echo "Log saved: $LOG_FILE" | tee -a "$LOG_FILE"

