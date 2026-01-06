#!/bin/bash
set -euo pipefail

# CPU threads: reserve 2, min 1
THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
echo "Running with $THREADS threads (Total: $(nproc))"
echo "------------------------------------------------------"
echo "STEP: Host Removal and Viral Mapping"
echo "------------------------------------------------------"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLEAN_READS="$PROJECT_DIR/results/02_clean_reads"
NON_HOST_READS_DIR="$PROJECT_DIR/results/03_nonhuman_reads"

HUMAN_REF_DIR="$PROJECT_DIR/reference_genomes"
reference_genome="$HUMAN_REF_DIR/GCA_000001405.28_GRCh38.p13_genomic.fna"
index_prefix="$HUMAN_REF_DIR/human_index"

# Silence Conda initialization warnings
set +u
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null
conda activate host_removal_env 2>/dev/null
set -u

# Check Bowtie2 index
if [[ ! -f "${index_prefix}.1.bt2" ]]; then
    echo "ERROR: Bowtie2 index not found at ${index_prefix}.*.bt2"
    exit 1
fi

mkdir -p "$NON_HOST_READS_DIR"

PROCESSED_SAMPLES=0
SKIPPED_SAMPLES=0
FAILED_SAMPLES=0
SAMPLES_TO_PROCESS=()

shopt -s nullglob
r1_files=("$CLEAN_READS"/*_R1.trimmed.fastq.gz)
shopt -u nullglob

if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "No R1 files found in $CLEAN_READS"
else
    for forward_read in "${r1_files[@]}"; do
        sample=$(basename "$forward_read" _R1.trimmed.fastq.gz)
        out1="$NON_HOST_READS_DIR/${sample}_1.nonhost.fastq.gz"
        out2="$NON_HOST_READS_DIR/${sample}_2.nonhost.fastq.gz"

        if [[ -f "$out1" && -f "$out2" ]]; then
            echo "Skipping $sample (already processed)"
            SKIPPED_SAMPLES=$((SKIPPED_SAMPLES + 1))
        else
            SAMPLES_TO_PROCESS+=("$forward_read")
        fi
    done

    if [[ ${#SAMPLES_TO_PROCESS[@]} -gt 0 ]]; then
        echo "Starting host removal..."

        for forward_read in "${SAMPLES_TO_PROCESS[@]}"; do
            sample=$(basename "$forward_read" _R1.trimmed.fastq.gz)
            reverse_read="$CLEAN_READS/${sample}_R2.trimmed.fastq.gz"

            if [[ ! -f "$reverse_read" ]]; then
                echo "Missing R2 for $sample. Skipping."
                FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
                continue
            fi

            SAMPLE_ALIGNMENT_DIR="$NON_HOST_READS_DIR/alignment_stats/${sample}"
            mkdir -p "$SAMPLE_ALIGNMENT_DIR"

            bam="$SAMPLE_ALIGNMENT_DIR/${sample}.bam"
            flagstat="$SAMPLE_ALIGNMENT_DIR/${sample}.flagstat"
            out1="$NON_HOST_READS_DIR/${sample}_1.nonhost.fastq.gz"
            out2="$NON_HOST_READS_DIR/${sample}_2.nonhost.fastq.gz"
            tmp_prefix="$SAMPLE_ALIGNMENT_DIR/${sample}_tmp_unmapped"

            echo "Processing Sample $sample"

            # We redirect Bowtie2's stderr only to the log file, not the terminal
            bowtie2 \
                --very-sensitive-local \
                --score-min L,0,-0.6 \
                --ma 0 \
                -x "$index_prefix" \
                -1 "$forward_read" \
                -2 "$reverse_read" \
                -p "$THREADS" \
                --un-conc "$tmp_prefix" \
                -S - \
                2> "$SAMPLE_ALIGNMENT_DIR/${sample}_bowtie2.log" \
            | samtools view -@ "$THREADS" -bS - 2>/dev/null > "$bam"

            if [[ $? -ne 0 ]]; then
                echo "Error mapping $sample"
                rm -f "${tmp_prefix}.1" "${tmp_prefix}.2" || true
                FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
                continue
            fi

            gzip -c "${tmp_prefix}.1" > "$out1" 2>/dev/null
            gzip -c "${tmp_prefix}.2" > "$out2" 2>/dev/null
            rm -f "${tmp_prefix}.1" "${tmp_prefix}.2"

            samtools flagstat -@ "$THREADS" "$bam" > "$flagstat" 2>/dev/null

            PROCESSED_SAMPLES=$((PROCESSED_SAMPLES + 1))
        done
    fi
fi

echo "------------------------------------------------------"
echo "Host Removal Completed"
echo "Processed: $PROCESSED_SAMPLES"
echo "Skipped:   $SKIPPED_SAMPLES"
echo "Failed:    $FAILED_SAMPLES"
echo "Done."
