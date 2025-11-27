# Marburg Virus VSP Illumina Analysis Pipeline

This repository contains a  analysis pipeline for processing Illumina sequencing data of Marburg virus (MARV) samples. The workflow includes raw data quality control, host read removal, mapping, variant calling, consensus generation, phylogenetic analysis, and coverage assessment. The pipeline is modular, reproducible, and suitable for batch processing of multiple samples.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Steps](#pipeline-steps)
- [Software Requirements](#software-requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Logging](#logging)
- [Authors](#authors)

---

## Overview

This pipeline automates the analysis of Illumina sequencing reads for Marburg virus. It is optimized for:

- High-throughput batch processing
- Flexible use of computational resources (multi-threaded)
- Comprehensive logging for reproducibility
- Integration with downstream analysis (phylogenetics, consensus sequences, coverage statistics)

---

## Pipeline Steps

The workflow is organized into the following steps, each implemented as a batch script:

| Step | Script                      | Description                                                                                               |
| ---- | --------------------------- | --------------------------------------------------------------------------------------------------------- |
| 1    | `copy_marv_fasta.sh`        | Filters downloaded Marburg virus genomes ≥18,000 bp and prepares reference files.                         |
| 2    | `fastp_batch.sh`            | Quality trimming, adapter removal, and filtering of raw FASTQ reads.                                      |
| 3    | `host_removal_batch.sh`     | Removes host reads to retain viral sequences.                                                             |
| 4    | `mapping_batch.sh`          | Maps reads to the Marburg reference genome using Minimap2.                                                |
| 5    | `variant_calling_batch.sh`  | Calls variants using samtools mpileup and iVar.                                                           |
| 6    | `consensus_batch.sh`        | Generates consensus sequences from BAM files using iVar.                                                  |
| 7    | `coverage_batch.sh`         | Computes per-base coverage, genome coverage, and ambiguous bases.                                         |
| 8    | `msa_batch.sh`              | Combines reference and consensus sequences, aligns with MAFFT, trims with trimAl, generates CSV metadata. |
| 9    | `iqtree_batch.sh`           | Builds phylogenetic trees from MSA using IQ-TREE, restores original leaf names, and creates summary.      |
| 10   | `multiqc_batch.sh`          | Aggregates QC reports across all samples.                                                                 |
| 11   | `run_marg_full_pipeline.sh` | Launches the entire workflow in sequence, handling intermediate directories and logging.                  |

---

## Software Requirements

The pipeline uses the following software tools:

- [Git](https://git-scm.com/)
- [Conda](https://docs.conda.io/en/latest/)
- [fastp](https://github.com/OpenGene/fastp) – Quality control of FASTQ reads
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) – Host read removal
- [samtools](http://www.htslib.org/) – BAM/SAM manipulation
- [iVar](https://andersen-lab.github.io/ivar/html/) – Variant calling and consensus generation
- [minimap2](https://github.com/lh3/minimap2) – Read mapping
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) – Multiple sequence alignment
- [trimAl](http://trimal.cgenomics.org/) – Alignment trimming
- [IQ-TREE](http://www.iqtree.org/) – Phylogenetic analysis
- [MultiQC](https://multiqc.info/) – Aggregated QC reporting
- Standard UNIX utilities: `awk`, `grep`, `tr`, `wc`


**Conda environments** are used to manage dependencies:

---

## Installation

1. Clone this repository:
```bash
git clone https://github.com/betselotz/Marburg-Virus-VSP-Illumina-Analysis.git
cd Marburg-Virus-VSP-Illumina-Analysis
```

2. Create required Conda environments:
 1. FastQ Quality Control (fastp)
```bash
conda create -n fastp_env -c bioconda -y fastp && conda run -n fastp_env fastp --version
```
 2. Host Read Removal (Bowtie2)
```bash
conda create -n host_env -c bioconda -y bowtie2 && conda run -n host_env bowtie2 --version
```
 3. Mapping (Minimap2 + Samtools)
```bash
conda create -n mapping_env -c bioconda -y minimap2 samtools && conda run -n mapping_env minimap2 --version && conda run -n mapping_env samtools --version
```
 4. Variant Calling & Consensus (iVar + Samtools)
```bash
conda create -n ivar_env -c bioconda -y samtools ivar && conda run -n ivar_env samtools --version && conda run -n ivar_env ivar --help
```
 5. Multiple Sequence Alignment (MAFFT + trimAl)
```bash
conda create -n mafft_env -c bioconda -y mafft trimal && conda run -n mafft_env mafft --version && conda run -n mafft_env trimal -version
```
```bash
# Phylogenetic Analysis (IQ-TREE)
```bash
conda create -n iqtree_env -c bioconda -y iqtree && conda run -n iqtree_env iqtree3 --version
```
7. MultiQC Reporting
```bash
conda create -n multiqc_env -c bioconda -y multiqc && conda run -n multiqc_env multiqc --version
```
8. MSA with trimAl (MAFFT + trimAl)
```bash
conda create -n mafft_env -c bioconda -y mafft trimal && conda run -n mafft_env mafft --version && conda run -n mafft_env trimal -version
```
9. IQ-TREE (redundant if already done in #6, but included for clarity)
```bash
conda create -n iqtree_env -c bioconda -y iqtree && conda run -n iqtree_env iqtree3 --version
```

3. Prepare input directories:
```bash
raw_reads/                  # Raw FASTQ files
reference_genomes/MARV_downloads/   # Downloaded MARV genomes
reference_genomes/MARV_compare/     # Filtered reference genomes
```

## Usage
1. Place raw FASTQ files in a designated directory (e.g., raw_reads/).
2. Prepare the reference genome in reference_genomes/Marburg_reference.fasta.
3. Run the full pipeline:
```bash
bash scripts/run_marg_full_pipeline.sh
```
4. Alternatively, run individual steps as needed:
```bash
bash scripts/fastp_batch.sh
```
```bash
bash scripts/host_removal_batch.sh
```
```bash
bash scripts/mapping_batch.sh
```
```bash
bash scripts/variant_calling_batch.sh
```
```bash
bash scripts/consensus_batch.sh
```
```bash
bash scripts/coverage_batch.sh
```
```bash
bash scripts/msa_batch.sh
```
```bash
bash scripts/iqtree_batch.sh
```
```bash
bash scripts/multiqc_batch.sh
```

## Directory Structure
```bash
Marburg-Virus-VSP-Illumina-Analysis/
├── raw_reads/
│   ├── MARV_X_1.fastq.gz
│   └── MARV_X_2.fastq.gz                 
├── reference_genomes/
│   ├── NC_001608.4.fasta (reference)
│   ├── EF446131.1.fasta (ourgroup)
│   ├── MARV_downloads/
│   └── MARV_compare/
├── results/
│   ├── 01_fastp
│   ├── 02_clean_reads
│   ├── 03_nonhuman_reads/
│   ├── 04_mapped_bam/
│   ├── 05_variants/
│   ├── 06_consensus/
│   ├── 07_coverage/
│   ├── 09_msa/
│   └── 10_phylogeny/
├── scripts/                    
├── logs/                       
└── README.md                   
```

## Logging
Each script writes per-sample logs in logs/.
A summary table is generated for mapping, variant calling, consensus, and coverage statistics.
Logs capture runtime, errors, and pipeline decisions.

## Authors

Betselot Zerihun Ayano – Principal developer
GitHub: @betselotz

## License

This repository is open for academic and research use.


