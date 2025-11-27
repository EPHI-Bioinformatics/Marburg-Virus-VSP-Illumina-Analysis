#!/bin/bash
set -euo pipefail

# ========================================
# Full MSA Pipeline for Marburg Virus
# Supports long FASTA headers (e.g., "M23M | Homo sapiens | Ethiopia | 2025")
# Generates a CSV metadata file from headers
# ========================================

# -------- CPU AUTO-DETECTION --------
TOTAL_THREADS=$(nproc)
SAFE_THREADS=$(( TOTAL_THREADS > 4 ? TOTAL_THREADS - 2 : 2 ))
echo "🧠 System Threads: $TOTAL_THREADS  |  Safe Threads: $SAFE_THREADS"

# -------- Project Directories --------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONS_DIR="$PROJECT_DIR/results/06_consensus_selected"
REF_DIR="$PROJECT_DIR/reference_genomes/MARV_compare"
MSA_DIR="$PROJECT_DIR/results/09_msa"
OUTGROUP="$PROJECT_DIR/reference_genomes/EF446131.1.fasta"

mkdir -p "$MSA_DIR"

LOG_FILE="$MSA_DIR/$(date +%Y%m%d_%H%M%S)_msa_pipeline.log"

ALL_SEQ_INPUT="$MSA_DIR/all_sequences_input.fasta"
ALN_RAW="$MSA_DIR/all_sequences_raw_aligned.fasta"
ALN_TRIM="$MSA_DIR/all_sequences_aligned_trimmed.fasta"
METADATA_CSV="$MSA_DIR/msa_metadata.csv"

echo "⚙️ Starting full MSA Pipeline" | tee "$LOG_FILE"

# -------- Emergency Cleanup --------
cleanup() {
    echo "🧹 Cleanup triggered..." | tee -a "$LOG_FILE"
    rm -f "$ALL_SEQ_INPUT" "$ALN_RAW"
}
trap cleanup EXIT

# -------- Skip if final MSA exists --------
if [[ -s "$ALN_TRIM" ]]; then
    echo "✅ MSA already exists: $ALN_TRIM" | tee -a "$LOG_FILE"
    exit 0
fi

# -------- Validate Inputs --------
CONS_INPUT_FILES=("$CONS_DIR"/*.fa)

for f in "${CONS_INPUT_FILES[@]}" "$REF_DIR"/*.fasta "$OUTGROUP"; do
    if [[ ! -s "$f" ]]; then
        echo "❌ FATAL: Missing or empty FASTA: $f" | tee -a "$LOG_FILE"
        exit 1
    fi
    SEQ_LEN=$(grep -v ">" "$f" | tr -d '\n' | wc -c)
    if [[ $SEQ_LEN -lt 10000 ]]; then
        echo "❌ FATAL: Sequence too short: $f (${SEQ_LEN} bp)" | tee -a "$LOG_FILE"
        exit 1
    fi
done

echo "✅ All sequences validated." | tee -a "$LOG_FILE"

# -------- Combine Sequences --------
cat "${CONS_INPUT_FILES[@]}" "$REF_DIR"/*.fasta "$OUTGROUP" > "$ALL_SEQ_INPUT"

# -------- Ensure Unique Headers --------
awk '
/^>/ {
    full=$0
    count[full]++
    if (count[full] > 1) {
        $0 = full "_" count[full]
    }
}
{ print }
' "$ALL_SEQ_INPUT" > "$ALL_SEQ_INPUT.tmp"

mv "$ALL_SEQ_INPUT.tmp" "$ALL_SEQ_INPUT"
echo "✅ Headers preserved and made unique." | tee -a "$LOG_FILE"

# -------- Activate MAFFT environment --------
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mafft_env
set -u

# -------- MAFFT HIGH-ACCURACY --------
MAFFT_OPTS=(--globalpair --maxiterate 1000)
MAFFT_THREADS=2  # safe for 32 GB RAM

echo "🔥 Running MAFFT High-Accuracy Mode" | tee -a "$LOG_FILE"
if ! mafft "${MAFFT_OPTS[@]}" --adjustdirection --thread "$MAFFT_THREADS" \
      "$ALL_SEQ_INPUT" > "$ALN_RAW" 2>>"$LOG_FILE"; then
    echo "❌ MAFFT failed, retrying with safer mode..." | tee -a "$LOG_FILE"
    MAFFT_OPTS=(--localpair --maxiterate 1000)
    if ! mafft "${MAFFT_OPTS[@]}" --adjustdirection --thread "$SAFE_THREADS" \
          "$ALL_SEQ_INPUT" > "$ALN_RAW" 2>>"$LOG_FILE"; then
        echo "❌ FATAL: MAFFT failed in fallback mode." | tee -a "$LOG_FILE"
        exit 1
    fi
fi
echo "✅ MAFFT alignment complete." | tee -a "$LOG_FILE"

# -------- TRIMAL --------
echo "✂️ Trimming with trimAl..." | tee -a "$LOG_FILE"
if ! trimal -in "$ALN_RAW" -out "$ALN_TRIM" -fasta -automated1 2>>"$LOG_FILE"; then
    echo "❌ FATAL: trimAl failed" | tee -a "$LOG_FILE"
    exit 1
fi

# -------- Generate Metadata CSV --------
echo "Accession,Host,Country,Year" > "$METADATA_CSV"
grep "^>" "$ALN_TRIM" | sed 's/^>//' | while IFS='|' read -r acc host country year; do
    #

