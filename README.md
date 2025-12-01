# Marburg Virus VSP Illumina Analysis Pipeline

This repository contains an analysis pipeline for processing Illumina sequencing data of Marburg virus (MARV) samples. The workflow includes raw data quality control, host read removal, mapping, BAM QC, variant calling, consensus generation, coverage assessment, clade assignment and MultiQC reporting. The pipeline is modular, reproducible, and suitable for batch processing of multiple samples.
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
| 5    | `qualimap_batch.sh`         | Performs BAM QC using Qualimap on mapped BAM files.                                                       |
| 6    | `variant_calling_batch.sh`  | Calls variants using samtools mpileup and iVar.                                                           |
| 7    | `consensus_batch.sh`        | Generates consensus sequences from BAM files using iVar.                                                  |
| 8    | `coverage_batch.sh`         | Computes per-base coverage, genome coverage, and ambiguous bases.                                         |
| 9    | `nextclade_batch.sh`        | Analyze clades using its own database.                                                                    |
| 9    | `multiqc_batch.sh`          | Aggregates QC reports from Fastp, and Qualimap and also clade from nextclade.                             |
| 10   | `msa_batch.sh`              | Combines reference and consensus sequences, aligns with MAFFT, trims with trimAl, generates CSV metadata. |
| 11   | `iqtree_batch.sh`           | Builds phylogenetic trees from MSA using IQ-TREE, restores original leaf names, and creates summary.      |
| 12   | `run_marg_full_pipeline.sh` | Launches the entire workflow in sequence, handling intermediate directories and logging.                  |

 !Note: Steps 10 and 11 (msa_batch.sh and iqtree_batch.sh) are not included in the automated full pipeline; they should be run separately after step 9 if phylogenetic analysis is needed.
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
- [Qualimap](http://qualimap.bioinfo.cipf.es) – BAM QC
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) – Multiple sequence alignment
- [trimAl](http://trimal.cgenomics.org/) – Alignment trimming
- [Qualimap] 
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
conda create -n fastp_env -c bioconda -y fastp
```
 2. Host Read Removal (Bowtie2)
```bash
conda create -n host_env -c bioconda -y bowtie2
```
 3. Mapping (Minimap2 + Samtools)
```bash
conda create -n mapping_env -c bioconda -y minimap2 samtools
```
4. BAM QC (Qualimap)
```bash
conda create -n qualimap_env -c bioconda -y qualimap
```
5. Variant Calling & Consensus (iVar + Samtools)
```bash
conda create -n ivar_env -c bioconda -y samtools ivar
```
 6. Multiple Sequence Alignment (MAFFT + trimAl)
```bash
conda create -n mafft_env -c bioconda -y mafft trimal
```
7.  Clade identification
```bash
conda create -n nextclade_env -c bioconda -y nextclade
```
8. Phylogenetic Analysis (IQ-TREE)
```bash
conda create -n iqtree_env -c bioconda -y iqtree
```
9. MultiQC Reporting
```bash
conda create -n multiqc_env -c bioconda -y multiqc
```

3. Prepare input directories:
```bash
raw_reads/                  # Raw FASTQ files
reference_genomes/MARV_downloads/   # Downloaded MARV genomes
reference_genomes/MARV_compare/     # Filtered reference genomes from downloaded MARV genomes
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
bash scripts/qualimap_batch.sh
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
bash scripts/nextclade_batch.sh
```
```bash
bash scripts/multiqc_batch.sh
```
```bash
bash scripts/msa_batch.sh
```
```bash
bash scripts/iqtree_batch.sh
```


## Directory Structure
```bash
Marburg-Virus-VSP-Illumina-Analysis/
├── raw_reads/
│   ├── MARV_X_1.fastq.gz
│   └── MARV_X_2.fastq.gz
├── reference_genomes/
│   ├── NC_001608.4.fasta
│   ├── EF446131.1.fasta
│   ├── MARV_downloads/
│   └── MARV_compare/
├── database/nextclade_marburg_dataset/
│   ├── marburg_tree.json
│   ├── pathogen.json
│   ├── reference.fasta
│   ├── examples.fasta
│   ├── genome_annotation.gff3
│   ├── CHANGELOG.md
│   └── README.md
├── results/
│   ├── 01_fastp/
│   ├── 02_clean_reads/
│   ├── 03_nonhuman_reads/
│   ├── 04_mapped_bam/
│   ├── 05_mapping_qc/
│   ├── 06_variants/
│   ├── 07_consensus/
│   ├── 08_coverage/
│   ├── 09_nextclade/
│   ├── 10_multiqc/
│   │   ├── fastp/
│   │   ├── nextclade/
│   │   ├── qualimap/
│   │   └── combined/
│   ├── 11_msa/
│   └── 12_iqtree/
├── scripts/
├── logs/
└── README.md
             
```

## Logging
Each script writes per-sample logs in logs/.
A summary table is generated for mapping, variant calling, consensus, and coverage statistics.
Logs capture runtime, errors, and pipeline decisions.

## Authors

Betselot Zerihun Ayano
GitHub: @betselotz

Melak Getu Bire
GitHub: MelakG13

## License

This repository is open for academic and research use.


