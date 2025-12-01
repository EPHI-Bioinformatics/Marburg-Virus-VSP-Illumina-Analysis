#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_DIR="$PROJECT_DIR/results/04_mapped_bam"
CONS_DIR="$PROJECT_DIR/results/07_consensus"
COVERAGE_DIR="$PROJECT_DIR/results/08_coverage"
mkdir -p "$COVERAGE_DIR"

THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))

log() { echo -e "$(date '+%Y-%m-%d %H:%M:%S') | $*"; }
log_sample() { echo -e "  → $*"; }

count_N_bases() {
    local fasta="$1"
    [[ ! -f "$fasta" ]] && echo "0 0 0.00" && return
    local total n_count percent
    total=$(grep -v "^>" "$fasta" | tr -d '\n' | wc -c)
    n_count=$(grep -v "^>" "$fasta" | tr -d '\n' | tr 'a-z' 'A-Z' | tr -cd 'N' | wc -c)
    percent=$(awk -v n="$n_count" -v len="$total" 'BEGIN{printf "%.2f", (n/len)*100}')
    echo "$n_count $total $percent"
}

calculate_coverage_stats() {
    local depth="$1" glen="$2"
    awk -v glen="$glen" '
    BEGIN{min=1e9; max=0; sum=0; sumsq=0; count=0; c1=c10=c20=c30=c50=0;}
    {d=$3; sum+=d; sumsq+=d*d; count++; if(d<min)min=d; if(d>max)max=d;
     if(d>=1)c1++; if(d>=10)c10++; if(d>=20)c20++; if(d>=30)c30++; if(d>=50)c50++;}
    END{
        mean=(count>0)?sum/count:0;
        stddev=(count>1)?sqrt((sumsq-(sum*sum)/count)/(count-1)):0;
        cov1=(c1/glen)*100; cov10=(c10/glen)*100; cov20=(c20/glen)*100;
        cov30=(c30/glen)*100; cov50=(c50/glen)*100;
        printf("%.2f\t%.0f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", mean, mean, min, max, stddev, cov1, cov10, cov20, cov30, cov50);
    }' "$depth"
}

main() {
    summary="$COVERAGE_DIR/depth_summary.tsv"
    echo -e "Sample\tMeanDepth\tMedianDepth\tMinDepth\tMaxDepth\tStdDevDepth\tCoverage1x\tCoverage10x\tCoverage20x\tCoverage30x\tCoverage50x\tN_count\tGenome_length\tPercent_N" > "$summary"

    for bam in "$BAM_DIR"/*.sorted.bam; do
        sample=$(basename "$bam" .sorted.bam)
        depth_file="$COVERAGE_DIR/${sample}.depth"
        fasta="$CONS_DIR/${sample}.fa"

        [[ -s "$depth_file" && -s "$summary" && $(grep -c "^${sample}\t" "$summary") -gt 0 ]] && log_sample "$sample already computed, skipping." && continue
        [[ ! -f "$bam" ]] && log_sample "BAM not found, skipping $sample." && continue
        [[ ! -f "$fasta" ]] && log_sample "FASTA not found, skipping $sample." && continue

        log_sample "Processing $sample..."
        samtools depth -a -@ "$THREADS" "$bam" > "$depth_file"
        genome_len=$(grep -v "^>" "$fasta" | tr -d '\n' | wc -c)
        stats=$(calculate_coverage_stats "$depth_file" "$genome_len")
        read n_count genome_len percent_n <<< $(count_N_bases "$fasta")

        echo -e "${sample}\t${stats}\t${n_count}\t${genome_len}\t${percent_n}" >> "$summary"

        mean_depth=$(echo "$stats" | cut -f1)
        cov1=$(echo "$stats" | cut -f6)
        cov30=$(echo "$stats" | cut -f9)
        log_sample "Done: Mean=${mean_depth}x, >=1x=${cov1}%, >=30x=${cov30}%, N_count=${n_count}, Percent_N=${percent_n}%"
    done

    echo -e "\n✅ Coverage summary at $summary"
}

for cmd in samtools awk grep tr wc; do
    command -v $cmd >/dev/null 2>&1 || { echo "$cmd not found"; exit 1; }
done

main

