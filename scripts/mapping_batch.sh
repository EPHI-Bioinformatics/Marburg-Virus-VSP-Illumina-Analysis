#!/bin/bash
set -euo pipefail

# --- 1. FORCE CONDA INITIALIZATION ---
# This ensures bwa/samtools/bcftools are found
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mapping_env

# Check that tools are actually found after activation
if ! command -v bwa &> /dev/null; then
    echo "ERROR: bwa not found. Is it installed in 'mapping_env'?"
    exit 1
fi

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NONHOST_DIR="$BASE_DIR/results/03_nonhuman_reads"
BAM_DIR="$BASE_DIR/results/04_mapping_results/bam"
VCF_DIR="$BASE_DIR/results/05_variants"
REFERENCE="$BASE_DIR/reference_genomes/marburg_reference.fasta"

mkdir -p "$BAM_DIR" "$VCF_DIR"

for fq1 in "$NONHOST_DIR"/*_nonhost_R1.fastq.gz; do
    [ -f "$fq1" ] || continue 
    
    sample=$(basename "$fq1" _nonhost_R1.fastq.gz)
    fq2="${NONHOST_DIR}/${sample}_nonhost_R2.fastq.gz"
    outbam="$BAM_DIR/${sample}.sorted.bam"
    outvcf="$VCF_DIR/${sample}.vcf.gz"

    echo "--- Mapping $sample ---"
    bwa mem -t 4 "$REFERENCE" "$fq1" "$fq2" | samtools view -b -q 30 -F 4 - | samtools sort -o "$outbam" -
    samtools index "$outbam"

    echo "--- Calling variants $sample ---"
    bcftools mpileup -Ou -f "$REFERENCE" "$outbam" | bcftools call -mv -Oz -o "$outvcf"
    bcftools index "$outvcf"
done

echo "Finished."
