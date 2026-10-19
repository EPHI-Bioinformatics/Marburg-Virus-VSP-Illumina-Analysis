#!/usr/bin/env bash
# 05_snps.sh
set -euo pipefail

# --- CONFIGURATION ---
ETH_ALIGNMENT="../04_snps/ethiopian_sequences.fasta"
REF_GENBANK="../04_snps/MARV_reference.gb"
OUTDIR="../04_snps"

echo "[1/4] Identifying SNPs using SNP-sites..."
snp-sites -v -o "$OUTDIR/ethiopian_SNPs.vcf" "$ETH_ALIGNMENT"

echo "[2/4] Mapping SNPs to GenBank annotations (Synonymous vs Non-synonymous)..."
# Using snpEff (requires a local database built from the GenBank file)
snpEff -v Marburg_virus "$OUTDIR/ethiopian_SNPs.vcf" > "$OUTDIR/ethiopian_SNPs_annotated.txt"

echo "[3/4] Calculating Nucleotide Diversity (Pi) and Watterson's Theta..."
# Using pixy or python scripts (assuming a wrapper script exists)
python3 calculate_diversity.py --fasta "$ETH_ALIGNMENT" --out "$OUTDIR/diversity_metrics.tsv"

echo "[4/4] Estimating gene-level dN/dS ratios..."
# Using HyPhy (BUSTED or SLAC)
hyphy slac --alignment "$ETH_ALIGNMENT" --tree "$OUTDIR/ethiopian_ML.treefile" --output "$OUTDIR/dNdS_results.json"

echo "SNP and diversity analysis complete."
