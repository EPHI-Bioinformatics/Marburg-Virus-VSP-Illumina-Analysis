#!/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NONHOST_DIR="$BASE_DIR/results/03_nonhuman_reads"
BAM_DIR="$BASE_DIR/results/04_mapping_results/bam"
REFERENCE="$BASE_DIR/reference_genomes/marburg_reference.fasta"

LOG_FILE="$BASE_DIR/pipeline_mapping_bwa_log_$(date +%Y%m%d_%H%M%S).txt"
SAMPLE_LOG_DIR="$BASE_DIR/logs/mapping"
mkdir -p "$BAM_DIR" "$SAMPLE_LOG_DIR"

echo "Starting high-confidence mapping with BWA (MAPQ ≥ 30)." | tee -a "$LOG_FILE"
echo -e "Sample\tInput_Reads\tMapped_Reads\tMapping_Rate" >> "$LOG_FILE"

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mapping_env
set -u

export THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
BWA_EXEC=$(which bwa)
SAMTOOLS_EXEC=$(which samtools)

[[ -z "$BWA_EXEC" || -z "$SAMTOOLS_EXEC" ]] && echo "Missing tools." && exit 1
[[ ! -f "$REFERENCE" ]] && echo "Missing reference." && exit 1

# Index reference if needed
[[ ! -f "${REFERENCE}.bwt" ]] && "$BWA_EXEC" index "$REFERENCE"
[[ ! -f "${REFERENCE}.fai" ]] && "$SAMTOOLS_EXEC" faidx "$REFERENCE"

for fq1 in "$NONHOST_DIR"/*_1.nonhost.fastq.gz; do
    [[ "$fq1" = "$NONHOST_DIR/*_1.nonhost.fastq.gz" ]] && break

    sample=$(basename "$fq1" _1.nonhost.fastq.gz)
    fq2="${fq1%_1.nonhost.fastq.gz}_2.nonhost.fastq.gz"
    outbam="$BAM_DIR/${sample}.sorted.bam"
    outbai="${outbam}.bai"
    SAMPLE_LOG_FILE="$SAMPLE_LOG_DIR/${sample}_bwa_mapping.log"

    [[ -s "$outbai" ]] && echo "Skipping $sample" | tee -a "$LOG_FILE" && continue
    [[ ! -f "$fq2" ]] && echo -e "${sample}\t0\t0\t0" >> "$LOG_FILE" && continue

    START=$(date +%s)
    READ_COUNT=$(zcat "$fq1" | wc -l | awk '{print $1/4}')
    [[ "$READ_COUNT" -lt 1 ]] && echo -e "${sample}\t${READ_COUNT}\t0\t0" >> "$LOG_FILE" && continue

    (
        "$BWA_EXEC" mem -t "$THREADS" "$REFERENCE" "$fq1" "$fq2" | \
        "$SAMTOOLS_EXEC" view -b -q 30 -F 4 - | \
        "$SAMTOOLS_EXEC" sort -@ "$THREADS" -o "$outbam" -
    ) 2> "$SAMPLE_LOG_FILE"

    set +u
    M1="${PIPESTATUS[0]}"
    M2="${PIPESTATUS[1]}"
    M3="${PIPESTATUS[2]}"
    set -u

    [[ "$M1" -ne 0 || "$M2" -ne 0 || "$M3" -ne 0 ]] && \
        echo -e "${sample}\t${READ_COUNT}\tFAIL\tFAIL" >> "$LOG_FILE" && \
        rm -f "$outbam" && continue

    "$SAMTOOLS_EXEC" index "$outbam" || {
        echo -e "${sample}\t${READ_COUNT}\tSORTED_BAM\tINDEX_FAIL" >> "$LOG_FILE"
        continue
    }

    MAPPED=$("$SAMTOOLS_EXEC" view -c "$outbam" || echo 0)
    RATE=$(awk "BEGIN { if ($READ_COUNT==0) print 0; else printf \"%.2f\", $MAPPED/$READ_COUNT*100 }")

    echo -e "${sample}\t${READ_COUNT}\t${MAPPED}\t${RATE}" >> "$LOG_FILE"

    END=$(date +%s)
    echo "Completed $sample in $((END - START))s" | tee -a "$LOG_FILE"

done

echo "All samples processed. Log: $LOG_FILE" | tee -a "$LOG_FILE"

