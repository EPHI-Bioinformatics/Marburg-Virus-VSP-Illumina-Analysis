#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTP_OUT="$PROJECT_DIR/results/01_fastp"
MULTIQC_DIR="$PROJECT_DIR/results/08_multiqc"

mkdir -p "$MULTIQC_DIR"

if [[ -f "$MULTIQC_DIR/multiqc_report.html" ]]; then
    echo "MultiQC report exists, skipping..."
    exit 0
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate multiqc_env

multiqc "$FASTP_OUT" -o "$MULTIQC_DIR"

conda deactivate
echo "MultiQC completed. Report saved at: $MULTIQC_DIR/multiqc_report.html"

