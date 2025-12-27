#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# -------------------------------
# Input directories
# -------------------------------
FASTP_QC_DIR="$PROJECT_DIR/results/01_fastp"
QUALIMAP_DIR="$PROJECT_DIR/results/05_mapping_qc"
NEXTCLADE_DIR="$PROJECT_DIR/results/09_nextclade/nextclade"

# -------------------------------
# Output directories
# -------------------------------
MULTIQC_DIR="$PROJECT_DIR/results/13_multiqc"
MULTIQC_FASTP_DIR="$MULTIQC_DIR/fastp"
MULTIQC_QUALIMAP_DIR="$MULTIQC_DIR/qualimap"
MULTIQC_NEXTCLADE_DIR="$MULTIQC_DIR/nextclade"
MULTIQC_COMBINED_DIR="$MULTIQC_DIR/combined"

mkdir -p \
  "$MULTIQC_FASTP_DIR" \
  "$MULTIQC_QUALIMAP_DIR" \
  "$MULTIQC_NEXTCLADE_DIR" \
  "$MULTIQC_COMBINED_DIR"

# -------------------------------
# Activate Conda
# -------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate multiqc_env || {
  echo "ERROR: Conda environment 'multiqc_env' not found."
  exit 1
}

# -------------------------------
# Sample list (single source of truth)
# -------------------------------
SAMPLES=(
  ET_MARV_23
  ET_MARV_31
  ET_MARV_32
  ET_MARV_45
  ET_MARV_60
  ET_MARV_261
  ET_MARV_262
)

EXPECTED_N=${#SAMPLES[@]}

# -------------------------------
# Integrity checks
# -------------------------------
echo "🔎 Performing sample integrity checks..."

for s in "${SAMPLES[@]}"; do
  [[ -f "$FASTP_QC_DIR/${s}.fastp.json" ]] || { echo "❌ Missing fastp JSON for $s"; exit 1; }
  [[ -d "$QUALIMAP_DIR/${s}.sorted" ]] || { echo "❌ Missing Qualimap directory for $s"; exit 1; }
done

echo "✅ All $EXPECTED_N samples found in Fastp and Qualimap"

# -------------------------------
# Filter Nextclade safely
# -------------------------------
NEXTCLADE_FILE="$NEXTCLADE_DIR/nextclade.csv"
FILTERED_NEXTCLADE="$MULTIQC_COMBINED_DIR/nextclade_selected.csv"

echo "🧬 Filtering Nextclade to selected samples..."

awk -F';' '
BEGIN {
  split("ET_MARV_23 ET_MARV_31 ET_MARV_32 ET_MARV_45 ET_MARV_60 ET_MARV_261 ET_MARV_262", s, " ")
  for (i in s) keep[s[i]] = 1
}
NR==1 { print; next }
($2 in keep) { print }
' "$NEXTCLADE_FILE" > "$FILTERED_NEXTCLADE"

# Validate Nextclade filtering
FOUND_NEXTCLADE=$(awk -F';' 'NR>1 {print $2}' "$FILTERED_NEXTCLADE" | wc -l)

if [[ "$FOUND_NEXTCLADE" -ne "$EXPECTED_N" ]]; then
  echo "❌ Nextclade filtering error: expected $EXPECTED_N samples, found $FOUND_NEXTCLADE"
  exit 1
fi

echo "✅ Nextclade filtering successful ($FOUND_NEXTCLADE samples)"

# -------------------------------
# Run MultiQC: Fastp
# -------------------------------
if [[ ! -f "$MULTIQC_FASTP_DIR/multiqc_report.html" ]]; then
  echo "🚀 Running Fastp MultiQC..."
  multiqc "$FASTP_QC_DIR" -o "$MULTIQC_FASTP_DIR" --force
fi

# -------------------------------
# Run MultiQC: Qualimap
# -------------------------------
if [[ ! -f "$MULTIQC_QUALIMAP_DIR/multiqc_report.html" ]]; then
  echo "🚀 Running Qualimap MultiQC..."
  multiqc "$QUALIMAP_DIR" -o "$MULTIQC_QUALIMAP_DIR" --force
fi

# -------------------------------
# Run MultiQC: Nextclade (filtered)
# -------------------------------
if [[ ! -f "$MULTIQC_NEXTCLADE_DIR/multiqc_report.html" ]]; then
  echo "🚀 Running Nextclade MultiQC (filtered)..."
  multiqc "$FILTERED_NEXTCLADE" -o "$MULTIQC_NEXTCLADE_DIR" --force
fi

# -------------------------------
# Run Combined MultiQC (FINAL)
# -------------------------------
if [[ ! -f "$MULTIQC_COMBINED_DIR/multiqc_report.html" ]]; then
  echo "🚀 Running Combined MultiQC..."

  FASTP_FILES=( "$FASTP_QC_DIR"/*.fastp.json )
  QUALIMAP_DIRS=( "$QUALIMAP_DIR"/*.sorted )

  multiqc \
    "${FASTP_FILES[@]}" \
    "${QUALIMAP_DIRS[@]}" \
    "$FILTERED_NEXTCLADE" \
    -o "$MULTIQC_COMBINED_DIR" \
    --force
fi

# -------------------------------
# Deactivate Conda
# -------------------------------
conda deactivate

echo "🎉 MultiQC completed successfully — all selected samples retained."

