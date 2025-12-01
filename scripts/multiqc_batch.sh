#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Input directories
FASTP_QC_DIR="$PROJECT_DIR/results/01_fastp"
QUALIMAP_DIR="$PROJECT_DIR/results/05_mapping_qc"
NEXTCLADE_DIR="$PROJECT_DIR/results/12_nextclade"

# Output directories
MULTIQC_DIR="$PROJECT_DIR/results/10_multiqc"
MULTIQC_FASTP_DIR="$MULTIQC_DIR/fastp"
MULTIQC_QUALIMAP_DIR="$MULTIQC_DIR/qualimap"
MULTIQC_NEXTCLADE_DIR="$MULTIQC_DIR/nextclade"
MULTIQC_COMBINED_DIR="$MULTIQC_DIR/combined"

# Create output directories
mkdir -p "$MULTIQC_FASTP_DIR" "$MULTIQC_QUALIMAP_DIR" "$MULTIQC_NEXTCLADE_DIR" "$MULTIQC_COMBINED_DIR"

# Activate Conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate multiqc_env || { echo "ERROR: Conda environment 'multiqc_env' not found."; exit 1; }

# ----------------------- Run MultiQC on Fastp -----------------------
REPORT_FILE="$MULTIQC_FASTP_DIR/multiqc_report.html"
if [[ -f "$REPORT_FILE" ]]; then
    echo "✅ Fastp MultiQC report already exists at $REPORT_FILE, skipping..."
else
    echo "🚀 Running MultiQC on Fastp outputs..."
    multiqc "$FASTP_QC_DIR" -o "$MULTIQC_FASTP_DIR" --force
    echo "✅ Fastp MultiQC completed. Report saved at: $REPORT_FILE"
fi

# ----------------------- Run MultiQC on Qualimap --------------------
REPORT_FILE="$MULTIQC_QUALIMAP_DIR/multiqc_report.html"
if [[ -f "$REPORT_FILE" ]]; then
    echo "✅ Qualimap MultiQC report already exists at $REPORT_FILE, skipping..."
else
    echo "🚀 Running MultiQC on Qualimap outputs..."
    multiqc "$QUALIMAP_DIR" -o "$MULTIQC_QUALIMAP_DIR" --force
    echo "✅ Qualimap MultiQC completed. Report saved at: $REPORT_FILE"
fi

# ----------------------- Run MultiQC on Nextclade -------------------
REPORT_FILE="$MULTIQC_NEXTCLADE_DIR/multiqc_report.html"
if [[ -f "$REPORT_FILE" ]]; then
    echo "✅ Nextclade MultiQC report already exists at $REPORT_FILE, skipping..."
else
    echo "🚀 Running MultiQC on Nextclade outputs..."
    multiqc "$NEXTCLADE_DIR" -o "$MULTIQC_NEXTCLADE_DIR" --force
    echo "✅ Nextclade MultiQC completed. Report saved at: $REPORT_FILE"
fi

# ----------------------- Run Combined MultiQC for selected samples -----------------------
REPORT_FILE="$MULTIQC_COMBINED_DIR/multiqc_report.html"
if [[ -f "$REPORT_FILE" ]]; then
    echo "✅ Combined MultiQC report already exists at $REPORT_FILE, skipping..."
else
    echo "🚀 Running Combined MultiQC for selected samples..."

    # Define selected sample IDs
    SAMPLES=("ET_MARV_23" "ET_MARV_31" "ET_MARV_32" "ET_MARV_261" "ET_MARV_262")

    # Collect Fastp files
    FASTP_FILES=()
    for s in "${SAMPLES[@]}"; do
        FASTP_FILES+=($(find "$FASTP_QC_DIR" -type f -name "*$s*.json"))
    done

    # Collect Qualimap directories
    QUALIMAP_DIRS=()
    for s in "${SAMPLES[@]}"; do
        QUALIMAP_DIRS+=($(find "$QUALIMAP_DIR" -maxdepth 1 -type d -name "$s.sorted"))
    done

    # Collect Nextclade file (filtered CSV)
    NEXTCLADE_FILE="$NEXTCLADE_DIR/nextclade.csv"
    # MultiQC can read the CSV for all samples but we filter it to selected samples
    FILTERED_NEXTCLADE="$MULTIQC_COMBINED_DIR/nextclade_selected.csv"
    awk -v samples="$(IFS="|"; echo "${SAMPLES[*]}")" 'NR==1 || $1 ~ samples' "$NEXTCLADE_FILE" > "$FILTERED_NEXTCLADE"

    # Run MultiQC on selected files
    multiqc "${FASTP_FILES[@]}" "${QUALIMAP_DIRS[@]}" "$FILTERED_NEXTCLADE" -o "$MULTIQC_COMBINED_DIR" --force

    echo "✅ Combined MultiQC for selected samples completed. Report saved at: $REPORT_FILE"
fi

# Deactivate Conda
conda deactivate

