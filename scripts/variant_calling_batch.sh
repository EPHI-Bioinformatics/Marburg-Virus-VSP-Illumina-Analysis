#!/bin/bash
set -euo pipefail

# 1. Setup Directories
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_DIR="$PROJECT_DIR/results/04_mapping_results/bam"
VAR_DIR="$PROJECT_DIR/results/06_variants"
VCF_DIR="$PROJECT_DIR/results/07_vcf"
MARBURG_REFERENCE="$PROJECT_DIR/reference_genomes/Marburg_reference.fasta"

THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
echo "Using $THREADS threads"
echo "----------------------------------------------------"

# 2. Environment Setup
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ivar_env
set -u

mkdir -p "$VAR_DIR" "$VCF_DIR"

SAMTOOLS_EXEC=$(which samtools)
IVAR_EXEC=$(which ivar)
BCFTOOLS_EXEC=$(which bcftools)

if [[ -z "$SAMTOOLS_EXEC" || -z "$IVAR_EXEC" || -z "$BCFTOOLS_EXEC" ]]; then
    echo "FATAL: samtools, ivar, or bcftools not found. Ensure they are installed in ivar_env."
    exit 1
fi

# 3. Processing Loop
for sorted_bam in "$BAM_DIR"/*.sorted.bam; do
    [[ ! -f "$sorted_bam" ]] && echo "No .sorted.bam files found." && continue

    sample=$(basename "$sorted_bam" .sorted.bam)
    output_prefix="$VAR_DIR/${sample}_variants"
    output_tsv="${output_prefix}.tsv"
    output_vcf="$VCF_DIR/${sample}.vcf.gz"

    echo "Calling variants for $sample"

    # iVar: TSV calling for amplicon/viral depth analysis
    if ! "$SAMTOOLS_EXEC" mpileup -A -d 1000000 -B -Q 0 -f "$MARBURG_REFERENCE" "$sorted_bam" | \
        "$IVAR_EXEC" variants -r "$MARBURG_REFERENCE" -p "$output_prefix"; then
        echo "ERROR: iVar calling failed for $sample."
    fi

    # BCFtools: Haploid VCF calling (--ploidy 1)
    echo "Generating VCF for $sample (Haploid mode)"
    "$BCFTOOLS_EXEC" mpileup -Ou -f "$MARBURG_REFERENCE" "$sorted_bam" | \
    "$BCFTOOLS_EXEC" call -mv -Oz --ploidy 1 -o "$output_vcf"
    "$BCFTOOLS_EXEC" index "$output_vcf"

    echo "Done: $output_tsv and $output_vcf"
    echo "----------------------------------------------------"
done

conda deactivate
echo "Variant calling complete."
