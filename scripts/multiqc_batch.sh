#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Input directories
FASTP_QC_DIR="$PROJECT_DIR/results/01_fastp"
QUALIMAP_DIR="$PROJECT_DIR/results/05_mapping_qc"

# Output directories
MULTIQC_DIR="$PROJECT_DIR/results/09_multiqc"
MULTIQC_FASTP_DIR="$MULTIQC_DIR/fastp"
MULTIQC_QUALIMAP_DIR="$MULTIQC_DIR/qualimap"
MULTIQC_COMBINED_DIR="$MULTIQC_DIR/combined"

# Create output directories
mkdir -p "$MULTIQC_FASTP_DIR" "$MULTIQC_QUALIMAP_DIR" "$MULTIQC_COMBINED_DIR"

# Activate Conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate multiqc_env || { echo "ERROR: Conda environment 'multiqc_env' not found."; exit 1; }

# ----------------------- Run MultiQC on Fastp -----------------------
REPORT_FILE="$MULTIQC_FASTP_DIR/multiqc_report.html"
if [[ -f "$REPORT_FILE" ]]; then
    echo "✅ Fastp MultiQC report already exists at $REPORT_FILE, skipping..."
else
    echo "🚀 Running MultiQC on fastp outputs..."
    multiqc "$FASTP_QC_DIR" -o "$MULTIQC_FASTP_DIR" --force
    echo "✅ Fastp MultiQC completed. Report saved at: $REPORT_FILE"
fi

# ----------------------- Run MultiQC on Qualimap --------------------
REPORT_FILE="$MULTIQC_QUALIMAP_DIR/multiqc_report.html"
if [[ -f "$REPORT_FILE" ]]; then
    echo "✅ Qualimap MultiQC report already exists at $REPORT_FILE, skipping..."
else
    echo "🚀 Running MultiQC on qualimap outputs..."
    multiqc "$QUALIMAP_DIR" -o "$MULTIQC_QUALIMAP_DIR" --force
    echo "✅ Qualimap MultiQC completed. Report saved at: $REPORT_FILE"
fi

# ----------------------- Run Combined MultiQC ---------------------
REPORT_FILE="$MULTIQC_COMBINED_DIR/multiqc_report.html"
if [[ -f "$REPORT_FILE" ]]; then
    echo "✅ Combined MultiQC report already exists at $REPORT_FILE, skipping..."
else
    echo "🚀 Running MultiQC on both Fastp + Qualimap outputs..."
    multiqc "$FASTP_QC_DIR" "$QUALIMAP_DIR" -o "$MULTIQC_COMBINED_DIR" --force
    echo "✅ Combined MultiQC completed. Report saved at: $REPORT_FILE"
fi

conda deactivate

