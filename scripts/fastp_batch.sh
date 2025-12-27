#!/bin/bash
set -euo pipefail

# ----------------------- Configuration -----------------------
MIN_QUAL=20
MIN_LEN=50
SIZE_THRESHOLD=0.90
PROCESSED_SAMPLES=0
FAILED_SAMPLES=0

# ----------------------- Conda Environment -------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate fastp_env \
    || { echo "ERROR: Conda environment 'fastp_env' not found."; exit 1; }

command -v fastp >/dev/null \
    || { echo "ERROR: fastp not found in environment."; exit 1; }

# ----------------------- Thread Detection --------------------
AVAILABLE_THREADS=$(command -v nproc >/dev/null && nproc || echo 2)
THREADS=$(( AVAILABLE_THREADS < 1 ? 1 : AVAILABLE_THREADS ))

echo "Using $THREADS CPU threads (auto-detected)."
echo "--- fastp QC: Started ---"

# ----------------------- Directories --------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$PROJECT_DIR/raw_reads"
RESULTS_DIR="$PROJECT_DIR/results"
FASTP_QC_DIR="$RESULTS_DIR/01_fastp"
CLEAN_FASTQ_DIR="$RESULTS_DIR/02_clean_reads"

mkdir -p "$FASTP_QC_DIR" "$CLEAN_FASTQ_DIR"

# ----------------------- Main Loop ----------------------------
while read -r fq1; do
    sample=$(basename "$fq1" | sed -E 's/(_R1|_R2).*\.f(ast)?q(\.gz)?$//')
    fq2="${fq1/_R1/_R2}"

    # ---- Trimmed reads (UNCHANGED) ----
    out1="$CLEAN_FASTQ_DIR/${sample}_R1.trimmed.fastq.gz"
    out2="$CLEAN_FASTQ_DIR/${sample}_R2.trimmed.fastq.gz"

    # ---- fastp reports (FIXED for MultiQC) ----
    html="$FASTP_QC_DIR/${sample}.fastp.html"
    json="$FASTP_QC_DIR/${sample}.fastp.json"
    log="$FASTP_QC_DIR/${sample}.fastp.log"

    [[ -f "$fq2" ]] || { echo "Missing R2 for $sample. Skipping."; continue; }
    [[ -s "$fq1" && -s "$fq2" ]] || { echo "Empty file in $sample. Skipping."; continue; }

    fq1_size=$(stat -c%s "$fq1")
    fq2_size=$(stat -c%s "$fq2")

    if ! awk "BEGIN { exit ( $fq1_size < $fq2_size * $SIZE_THRESHOLD || $fq2_size < $fq1_size * $SIZE_THRESHOLD ) }"; then
        echo "Size mismatch for $sample"
    fi

    if [[ -s "$out1" && -s "$out2" ]]; then
        echo "Already processed: $sample"
        PROCESSED_SAMPLES=$((PROCESSED_SAMPLES + 1))
        continue
    fi

    echo "Processing $sample..."

    (
        fastp \
            -i "$fq1" -I "$fq2" \
            -o "$out1" -O "$out2" \
            -h "$html" -j "$json" \
            --report_title "$sample" \
            --thread "$THREADS" \
            --detect_adapter_for_pe \
            --length_required "$MIN_LEN" \
            --qualified_quality_phred "$MIN_QUAL" \
            --trim_poly_g --trim_poly_x \
            --low_complexity_filter \
            --cut_window_size 4 \
            --cut_mean_quality 20
    ) 2> "$log" \
        && PROCESSED_SAMPLES=$((PROCESSED_SAMPLES + 1)) \
        || { echo "fastp failed for $sample"; FAILED_SAMPLES=$((FAILED_SAMPLES + 1)); }

    echo "Done $sample."
done < <(find "$RAW_DIR" -maxdepth 1 -name "*_R1*.fastq*")

# ----------------------- Summary ------------------------------
echo "--- fastp QC: Completed ---"
echo "Processed: $PROCESSED_SAMPLES"
echo "Failed:    $FAILED_SAMPLES"

