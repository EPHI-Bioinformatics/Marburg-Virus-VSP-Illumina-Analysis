#!/bin/bash
set -euo pipefail

# ---------------------------
# CONFIGURATION
# ---------------------------
WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOWNLOAD_DIR="$WORKDIR/reference_genomes/MARV_downloads"
COMPARE_DIR="$WORKDIR/reference_genomes/MARV_compare"
SUMMARY_FILE="$COMPARE_DIR/selection_summary.tsv"
TEMP_FILE="$COMPARE_DIR/selection_temp.tsv"

MIN_LENGTH=18000  

mkdir -p "$COMPARE_DIR"

echo -e "Genome\tLength(bp)\tIncluded" > "$TEMP_FILE"
echo "Filtering Marburg virus genomes by length (≥ $MIN_LENGTH bp)..."

for fasta in "$DOWNLOAD_DIR"/*.fasta; do
    if [[ ! -f "$fasta" ]]; then
        echo "No FASTA files found in $DOWNLOAD_DIR"
        exit 1
    fi

    # Calculate genome length (ignore FASTA headers)
    length=$(grep -v "^>" "$fasta" | tr -d '\n' | wc -c)

    if (( length >= MIN_LENGTH )); then
        cp "$fasta" "$COMPARE_DIR/"
        echo -e "$(basename "$fasta")\t$length\tYES" >> "$TEMP_FILE"
        echo "✅ Copied $(basename "$fasta") (length: $length bp)"
    else
        echo -e "$(basename "$fasta")\t$length\tNO" >> "$TEMP_FILE"
        echo "⚠️ Skipped $(basename "$fasta") (length: $length bp)"
    fi
done

# Sort summary by length descending and save final summary
(head -n 1 "$TEMP_FILE" && tail -n +2 "$TEMP_FILE" | sort -k2,2nr) > "$SUMMARY_FILE"
rm "$TEMP_FILE"

echo "Filtering complete. Selected sequences are in $COMPARE_DIR"
echo "Summary written to $SUMMARY_FILE (sorted by length descending)"

