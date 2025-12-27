#!/bin/bash
set -euo pipefail

###############################################################################
# 🧬 MARV-Gen: Marburg Viral Genomics & Phylodynamics Pipeline
###############################################################################

# ============================== CONFIG =======================================

# 1. Get the directory where THIS script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 2. Intelligent Root Detection
if [[ "$(basename "$SCRIPT_DIR")" == "scripts" ]]; then
    PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
else
    PROJECT_DIR="$SCRIPT_DIR"
fi

CPU=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 4)
export CPU

RUN_STEPS="${RUN_STEPS:-all}"
RUN_STEP12="${RUN_STEP12:-false}"

LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/marv_gen_link_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ============================== COLORS =======================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
BLUE='\033[0;34m'
NC='\033[0m' 
BOLD='\033[1m'

# ============================ DIRECTORIES ====================================
SCRIPTS_DIR="$PROJECT_DIR/scripts"
ENVS_DIR="$PROJECT_DIR/envs"
RAW_DATA="$PROJECT_DIR/raw_reads"
REF_DIR="$PROJECT_DIR/reference_genomes"

# Result Directories
FASTP_DIR="$PROJECT_DIR/results/01_fastp"
HOSTREM_DIR="$PROJECT_DIR/results/02_clean_reads"
NONHOST_DIR="$PROJECT_DIR/results/03_nonhuman_reads"
MAP_DIR="$PROJECT_DIR/results/04_mapping_results"
MAPPING_QC_DIR="$PROJECT_DIR/results/05_mapping_qc"
VAR_DIR="$PROJECT_DIR/results/06_variants"
CONS_DIR="$PROJECT_DIR/results/07_consensus"
COV_DIR="$PROJECT_DIR/results/08_coverage"
NEXTCLADE_DIR="$PROJECT_DIR/results/09_nextclade"
MSA_DIR="$PROJECT_DIR/results/10_msa"
IQTREE_DIR="$PROJECT_DIR/results/11_phylogeny"
TREETIME_DIR="$PROJECT_DIR/results/12_treetime"
MULTIQC_DIR="$PROJECT_DIR/results/13_multiqc"

# Create directories
mkdir -p "$LOG_DIR" "$ENVS_DIR" "$FASTP_DIR" "$HOSTREM_DIR" "$NONHOST_DIR" \
         "$MAP_DIR" "$MAPPING_QC_DIR" "$VAR_DIR" "$CONS_DIR" "$COV_DIR" \
         "$NEXTCLADE_DIR" "$MSA_DIR" "$IQTREE_DIR" \
         "$TREETIME_DIR" "$MULTIQC_DIR"

# ============================== FUNCTIONS ====================================

check_env() {
    local env=$1
    local yml_file="$ENVS_DIR/${env}.yml"

    echo -n -e "  📦 Checking ${BOLD}$env${NC}... "

    set +u 
    if conda info --envs | grep -q "$env"; then
        echo -e "${GREEN}✅ READY${NC}"
        set -u
        return 0
    else
        echo -e "${YELLOW}⚠️  MISSING${NC}"
        if [[ -f "$yml_file" ]]; then
            echo -e "      🚀 Attempting Auto-Install from $yml_file..."
            if conda env create -n "$env" -f "$yml_file" > /dev/null 2>&1; then
                echo -e "      ${GREEN}✨ $env installed successfully!${NC}"
                set -u
                return 0
            else
                echo -e "      ${RED}❌ Installation FAILED for $env.${NC}"
                set -u
                return 1
            fi
        else
            echo -e "      ${RED}❌ FAILED: $yml_file not found.${NC}"
            set -u
            return 1
        fi
    fi
}

print_step() {
    echo -e "\n${BLUE}┌──────────────────────────────────────────────────────────┐${NC}"
    echo -e "${BLUE}│${NC}  ${BOLD}${MAGENTA}▶ $1${NC}"
    echo -e "${BLUE}└──────────────────────────────────────────────────────────┘${NC}"
}

run_if_needed() {
    local step_name=$1
    local command=$2

    if [[ "$RUN_STEPS" != "all" && ! "$RUN_STEPS" =~ $step_name ]]; then
        echo -e "  ${CYAN}⏩ Skipping: $step_name (User Excluded)${NC}"
        return
    fi

    print_step "$step_name"
    echo -e "  ${CYAN}🔍 Scanning for missing samples...${NC}"
    
    if eval "$command"; then
         echo -e "  ${GREEN}✨ $step_name processing segment complete.${NC}"
    else
         echo -e "  ${RED}❌ Error in $step_name. Check logs.${NC}"
         exit 1
    fi
}

# ============================== PIPELINE START ===============================

clear
echo -e "${CYAN}${BOLD}"
echo "    __  ___ ___     ____ _    __      ______ ______ _   __ "
echo "   /  |/  //   |   / __ \ |  / /     / ____// ____// | / / "
echo "  / /|_/ // /| |  / /_/ / | / /_____/ / __ / __/  /  |/ /  "
echo " / /  / // ___ | / _, _/| |/ /_____/ /_/ // /___ / /|  /   "
echo "/_/  /_//_/  |_|/_/ |_| |___/       \____//_____//_/ |_/    "
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "  🧬 ${BOLD}Pipeline:${NC}  MARV-Gen"
echo -e "  📅 ${BOLD}Started:${NC}    $(date)"
echo -e "  📂 ${BOLD}Base Dir:${NC}  ${PROJECT_DIR}"
echo -e "  📄 ${BOLD}Log File:${NC}  ${YELLOW}${LOG_FILE}${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

# --- 1. DATA VALIDATION ---
echo -e "${YELLOW}${BOLD}🔍 STAGE 1: Data Validation${NC}"
if [[ ! -d "$RAW_DATA" ]] || [[ -z "$(ls -A "$RAW_DATA"/*.fastq.gz 2>/dev/null)" ]]; then
    echo -e "  ${RED}❌ ERROR: No .fastq.gz files found in $RAW_DATA${NC}"
    exit 1
fi
echo -e "  ${GREEN}✅ Data validation passed.${NC}"

# --- 2. ENVIRONMENT HEALTH CHECK ---
echo -e "\n${YELLOW}${BOLD}🔍 STAGE 2: Environment Health Check${NC}"
ENV_ERRORS=0
# Note: Changed msa_env to mafft_env here
envs=("fastp_env" "host_removal_env" "mapping_env" "qualimap_env" "ivar_env" "nextclade_env" "mafft_env" "iqtree_env" "treetime_env" "multiqc_env")

for e in "${envs[@]}"; do
    check_env "$e" || ENV_ERRORS=$((ENV_ERRORS+1))
done

if [ $ENV_ERRORS -gt 0 ]; then
    echo -e "\n${RED}❌ FATAL: $ENV_ERRORS environments could not be initialized.${NC}"
    exit 1
fi
echo -e "  ${GREEN}✅ Environment check passed.${NC}"

# ============================== EXECUTION ====================================

run_if_needed "STEP 1: Fastp QC" "bash $SCRIPTS_DIR/fastp_batch.sh"
run_if_needed "STEP 2: Host Removal" "bash $SCRIPTS_DIR/host_removal_batch.sh"
run_if_needed "STEP 3: Mapping Non-Host Reads" "bash $SCRIPTS_DIR/mapping_batch.sh"
run_if_needed "STEP 4: Qualimap BAM QC" "bash $SCRIPTS_DIR/qualimap_qc.sh"
run_if_needed "STEP 5: Variant Calling" "bash $SCRIPTS_DIR/variant_calling_batch.sh"
run_if_needed "STEP 6: Consensus Generation" "bash $SCRIPTS_DIR/consensus_batch.sh"
run_if_needed "STEP 7: Coverage Calculation" "bash $SCRIPTS_DIR/coverage_batch.sh"
run_if_needed "STEP 8: Nextclade" "bash $SCRIPTS_DIR/nextclade_batch.sh"
run_if_needed "STEP 10: MSA" "bash $SCRIPTS_DIR/msa_batch.sh"
run_if_needed "STEP 11: Phylogeny" "bash $SCRIPTS_DIR/iqtree_batch.sh"
run_if_needed "STEP 12: TreeTime" "bash $SCRIPTS_DIR/treetime_batch.sh"
run_if_needed "STEP 13: MultiQC" "bash $SCRIPTS_DIR/multiqc_batch.sh"
 

echo -e "\n${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "  🎉 ${MAGENTA}${BOLD}MARV-Gen pipeline completed!${NC}"
echo -e "  📂 Final results in: ${PROJECT_DIR}/results/"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
