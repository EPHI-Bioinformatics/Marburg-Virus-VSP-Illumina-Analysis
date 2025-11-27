# Marburg Virus VSP Illumina Analysis Pipeline

This repository contains a full analysis pipeline for processing Illumina sequencing data of Marburg virus (MARV) samples. The workflow includes raw data quality control, host read removal, mapping, variant calling, consensus generation, phylogenetic analysis, and coverage assessment. The pipeline is designed to be modular, reproducible, and suitable for batch processing of multiple samples.

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

1. **FastQ Quality Control (`fastp_batch.sh`)**  
   Performs quality trimming, adapter removal, and filtering of raw FASTQ reads.

2. **Host Read Removal (`host_removal_batch.sh`)**  
   Removes human or other host reads to retain only non-host sequences.

3. **Mapping (`mapping_batch.sh`)**  
   Maps non-host reads to the Marburg reference genome using **Minimap2**, producing sorted BAM files.

4. **Variant Calling (`variant_calling_batch.sh`)**  
   Calls variants from mapped BAM files using **samtools mpileup** and **iVar**, producing TSV variant files.

5. **Consensus Generation (`consensus_batch.sh`)**  
   Generates consensus sequences from BAM files using iVar, with customizable thresholds for base calling.

6. **Coverage Analysis (`coverage_batch.sh`)**  
   Computes per-base depth, genome coverage statistics, and the proportion of ambiguous bases (N).

7. **Multiple Sequence Alignment (`msa_batch.sh`)**  
   Aligns consensus sequences using **MAFFT** for downstream phylogenetic analysis.

8. **Phylogenetic Analysis (`iqtree_batch.sh`)**  
   Builds phylogenetic trees from multiple sequence alignments using **IQ-TREE**.

9. **MultiQC Reporting (`multiqc_batch.sh`)**  
   Aggregates QC reports across all samples for easy review.

10. **Pipeline Launcher (`run_marg_full_pipeline.sh`)**  
    Executes the entire workflow in sequence, handling all intermediate directories and logging.

11. **Utility Script (`copy_marv_fasta.sh`)**  
    Copies and prepares reference and consensus FASTA files for downstream analyses.

---

## Software Requirements

The pipeline uses the following software tools:

- [Git](https://git-scm.com/)
- [Conda](https://docs.conda.io/en/latest/)
- [samtools](http://www.htslib.org/)
- [iVar](https://andersen-lab.github.io/ivar/html/)
- [minimap2](https://github.com/lh3/minimap2)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [IQ-TREE](http://www.iqtree.org/)
- [MultiQC](https://multiqc.info/)
- Standard UNIX utilities: `awk`, `grep`, `tr`, `wc`

**Conda environments** are used to manage dependencies:

- `mapping_env` → for Minimap2 mapping
- `ivar_env` → for variant calling and consensus generation

---

## Installation

1. Clone this repository:
```bash
git clone https://github.com/betselotz/Marburg-Virus-VSP-Illumina-Analysis.git
cd Marburg-Virus-VSP-Illumina-Analysis
```


git clone https://github.com/betselotz/Marburg-Virus-VSP-Illumina-Analysis.git
cd Marburg-Virus-VSP-Illumina-Analysis
