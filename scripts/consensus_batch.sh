#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_DIR="$PROJECT_DIR/results/04_mapping_results/bam"
CONS_DIR="$PROJECT_DIR/results/07_consensus"
LOG_FILE="$PROJECT_DIR/pipeline_consensus_log_$(date +%Y%m%d_%H%M%S).txt"

mkdir -p "$CONS_DIR"
export THREADS=$(( $(nproc) > 4 ? 4 : 1 ))

echo "Starting consensus generation (Accurate Headers & Unlimited Depth)." | tee -a "$LOG_FILE"
echo "--------------------------------------------------------" | tee -a "$LOG_FILE"

# Initialize Conda
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ivar_env
set -u

BAM_FILES=("$BAM_DIR"/*.sorted.bam)
if [[ ! -e "${BAM_FILES[0]}" ]]; then
    echo "FATAL: No *.sorted.bam files in $BAM_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

for bam in "${BAM_FILES[@]}"; do
    START_TIME=$(date +%s)
    
    # 1. Clean the sample name: replace spaces/dots with underscores, remove trailing underscores
    raw_sample=$(basename "$bam" .sorted.bam)
    sample=$(echo "$raw_sample" | tr '[:space:].' '_' | sed 's/_*$//')
    
    consensus_file="$CONS_DIR/${sample}.fa"

    if [[ -s "$consensus_file" ]]; then
        echo "Existing consensus for $sample, skipping." | tee -a "$LOG_FILE"
        continue
    fi

    echo "Processing $sample..." | tee -a "$LOG_FILE"

    # 2. Unlimited depth (-d 0) ensures high-coverage viral sites aren't capped
    samtools mpileup -A -d 0 -Q 20 "$bam" | \
    ivar consensus -p "$CONS_DIR/$sample" -q 20 -t 0.7 -m 1 -n N 2>> "$LOG_FILE"

    if [[ -f "$CONS_DIR/${sample}.fa" ]]; then
        # 3. Final header cleanup: apply sanitized name
        sed -i "1s/.*/>${sample}/" "$CONS_DIR/${sample}.fa"
        echo "Consensus saved: $CONS_DIR/${sample}.fa" | tee -a "$LOG_FILE"
    else
        echo "WARNING: Consensus failed for $sample" | tee -a "$LOG_FILE"
    fi

    END_TIME=$(date +%s)
    echo "Runtime: $((END_TIME - START_TIME))s" | tee -a "$LOG_FILE"
done

conda deactivate
echo "Done. Total files: $(ls "$CONS_DIR"/*.fa 2>/dev/null | wc -l)" | tee -a "$LOG_FILE"

