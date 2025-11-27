#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_DIR="$PROJECT_DIR/results/04_mapped_bam"
VAR_DIR="$PROJECT_DIR/results/05_variants"
MARBURG_REFERENCE="$PROJECT_DIR/reference_genomes/Marburg_reference.fasta"

THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
echo "Using $THREADS threads"
echo "----------------------------------------------------"

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ivar_env
set -u

mkdir -p "$VAR_DIR"

SAMTOOLS_EXEC=$(which samtools)
IVAR_EXEC=$(which ivar)

if [[ -z "$SAMTOOLS_EXEC" || -z "$IVAR_EXEC" ]]; then
    echo "FATAL: samtools or ivar not found."
    exit 1
fi

for sorted_bam in "$BAM_DIR"/*.sorted.bam; do
    [[ ! -f "$sorted_bam" ]] && echo "No .sorted.bam files found." && continue

    sample=$(basename "$sorted_bam" .sorted.bam)
    output_prefix="$VAR_DIR/${sample}_variants"
    output_file="${output_prefix}.tsv"

    [[ -s "$output_file" ]] && echo "Existing variants for $sample. Skipping." && continue

    echo "Calling variants for $sample"

    if ! "$SAMTOOLS_EXEC" mpileup -A -d 1000000 -B -Q 0 -f "$MARBURG_REFERENCE" "$sorted_bam" | \
        "$IVAR_EXEC" variants -r "$MARBURG_REFERENCE" -p "$output_prefix"; then
        echo "ERROR: Variant calling failed for $sample."
        continue
    fi

    echo "Done: $output_file"
    echo "----------------------------------------------------"
done

conda deactivate
echo "Variant calling complete."

