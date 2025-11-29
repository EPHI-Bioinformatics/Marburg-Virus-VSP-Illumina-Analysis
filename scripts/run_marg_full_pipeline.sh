#!/bin/bash
set -euo pipefail

###############################################################################
# 🧬 MARG2 Viral Genomics Pipeline
# Author: Betselot Zerihun Ayano
#
# Description:
#   Full automated pipeline for Illumina-based viral genome analysis.
#   Runs sequential modules including QC, host removal, mapping,
#   BAM QC (Qualimap), variant calling, consensus generation, coverage,
#   and MultiQC reporting.
#
# Features:
#   ✔ Clean interactive terminal output
#   ✔ Spinners for long-running steps
#   ✔ Reproducible directory structure
#   ✔ Modular design (scripts/ subdirectory)
###############################################################################

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

# Directories
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS_DIR="$PROJECT_DIR/scripts"

FASTP_DIR="$PROJECT_DIR/results/01_fastp"
HOSTREM_DIR="$PROJECT_DIR/results/02_clean_reads"
NONHOST_DIR="$PROJECT_DIR/results/03_nonhuman_reads"
MAP_DIR="$PROJECT_DIR/results/04_mapped_bam"
MAPPING_QC_DIR="$PROJECT_DIR/results/05_mapping_qc"
VAR_DIR="$PROJECT_DIR/results/06_variants"
CONS_DIR="$PROJECT_DIR/results/07_consensus"
COV_DIR="$PROJECT_DIR/results/08_coverage"
MULTIQC_DIR="$PROJECT_DIR/results/09_multiqc"

mkdir -p "$FASTP_DIR" "$HOSTREM_DIR" "$NONHOST_DIR" "$MAP_DIR" "$MAPPING_QC_DIR" \
         "$VAR_DIR" "$CONS_DIR" "$COV_DIR" "$MULTIQC_DIR"

# Spinner
spinner() {
    local pid=$1
    local delay=0.1
    local spinstr='|/-\'
    while ps -p $pid > /dev/null 2>&1; do
        for i in $(seq 0 3); do
            printf "\r${CYAN}Processing... %s${NC}" "${spinstr:$i:1}"
            sleep $delay
        done
    done
    printf "\r"
}

# Section header
print_step() {
    echo -e "${YELLOW}========================================${NC}"
    echo -e "${YELLOW} $1 ${NC}"
    echo -e "${YELLOW}========================================${NC}"
}

echo -e "🚀 ${CYAN}Starting full MARG2 pipeline at $(date)${NC}\n"

# STEP 1: Fastp QC and trimming
print_step "STEP 1: Fastp QC and trimming"
bash "$SCRIPTS_DIR/fastp_batch.sh"
echo ""

# STEP 2: Host removal
print_step "STEP 2: Host removal"
bash "$SCRIPTS_DIR/host_removal_batch.sh" &
HOST_PID=$!
spinner $HOST_PID
wait $HOST_PID
echo -e "\r${GREEN}✅ Host removal completed.${NC}\n"

# STEP 3: Mapping non-host reads
print_step "STEP 3: Mapping non-host reads"
bash "$SCRIPTS_DIR/mapping_batch.sh"
echo -e "${GREEN}✅ All non-host reads mapped.${NC}\n"

# STEP 4: Qualimap BAM QC
print_step "STEP 4: Qualimap BAM QC"
bash "$SCRIPTS_DIR/qualimap_qc.sh"
echo -e "${GREEN}✅ Qualimap BAM QC completed.${NC}\n"

# STEP 5: Variant calling (iVar)
print_step "STEP 5: Variant calling"
bash "$SCRIPTS_DIR/variant_calling_batch.sh"
echo -e "${GREEN}✅ Variant calling completed.${NC}\n"

# STEP 6: Consensus generation (iVar)
print_step "STEP 6: Consensus generation"
bash "$SCRIPTS_DIR/consensus_batch.sh"
echo -e "${GREEN}✅ Consensus generation completed.${NC}\n"

# STEP 7: Coverage calculation
print_step "STEP 7: Coverage calculation"
bash "$SCRIPTS_DIR/coverage_batch.sh" &
COV_PID=$!
spinner $COV_PID
wait $COV_PID
echo -e "\r${GREEN}✅ Coverage calculation completed.${NC}\n"

# STEP 8: MultiQC (Fastp + Qualimap + Combined)
print_step "STEP 8: MultiQC reporting"
bash "$SCRIPTS_DIR/multiqc_batch.sh"
echo -e "${GREEN}✅ MultiQC reporting completed.${NC}\n"

echo -e "🎉 ${MAGENTA}All steps completed successfully!${NC}"

