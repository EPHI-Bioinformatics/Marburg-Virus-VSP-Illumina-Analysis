# MARV-GEN
## Marburg Virus VSP Illumina Analysis Pipeline
This repository contains an analysis pipeline for processing Illumina sequencing data of Marburg virus (MARV) samples. The workflow includes raw data quality control, host read removal, mapping, BAM QC, variant calling, consensus generation, coverage assessment, clade assignment, phylogenetic analysis, MultiQC reporting, and downstream genomic analyses including lineage-defining SNP identification, diversity metrics, selective pressure estimation, and publication-quality visualizations.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Steps](#pipeline-steps)
- [Software Requirements](#software-requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Lineage Defining SNPs Analysis](#SNPs-analysis)
- [Visualization](#Visualization)
- [Logging](#logging)
- [Authors](#authors)

---

## Overview

This pipeline automates the analysis of Illumina sequencing reads for Marburg virus. It is optimized for:

- High-throughput batch processing
- Flexible use of computational resources (multi-threaded)
- Comprehensive logging for reproducibility
- Integration with downstream analyses (phylogenetics, SNPs, diversity, dN/dS, genome visualization)
- Aggregated QC reporting using MultiQC

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
| 10   | `msa_batch.sh`              | Combines reference and consensus sequences, aligns with MAFFT, trims with trimAl, generates CSV metadata. |
| 11   | `iqtree_batch.sh`           | Builds phylogenetic trees from MSA using IQ-TREE, restores original leaf names, and creates summary.      |
| 11   | `treetime_batch.sh`         | Generates time-resolved phylogenies using TreeTime from IQ-TREE output.                                   |
| 11   | `multiqc_batch.sh`          | Aggregates QC reports from Fastp, Qualimap, and Nextclade; filters Nextclade to selected samples.         |
| 12   | `run_MARV-GEN_full_pipeline.sh` | Launches the entire workflow in sequence, handling intermediate directories and logging.              |

 
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
- [Nextclade](https://clades.nextstrain.org/) – clade asignment
- [IQ-TREE](http://www.iqtree.org/) – Phylogenetic analysis
- [treetime](https://github.com/neherlab/treetime) – time-scaled phylogenies
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
 I. FastQ Quality Control (fastp)
```bash
conda create -n fastp_env -c bioconda -y fastp
```
 II. Host Read Removal (Bowtie2)
```bash
conda create -n host_env -c bioconda -y bowtie2
```
 III. Mapping (Minimap2 + Samtools)
```bash
conda create -n mapping_env -c bioconda -y minimap2 samtools
```
IV. BAM QC (Qualimap)
```bash
conda create -n qualimap_env -c bioconda -y qualimap
```
V. Variant Calling & Consensus (iVar + Samtools)
```bash
conda create -n ivar_env -c bioconda -y samtools ivar
```
VI. Multiple Sequence Alignment (MAFFT)
```bash
conda create -n mafft_env -c bioconda -y mafft trimal
```
VII. Clade Identification (Nextclade)
```bash
conda create -n nextclade_env -c bioconda -y nextclade
```
VIII. Phylogenetic Analysis (IQ-TREE)
```bash
conda create -n iqtree_env -c bioconda -y iqtree
```
IX. Time-Resolved Phylogeny (TreeTime)
```bash
conda create -n treetime_env -c bioconda -y treetime
```
X. MultiQC Reporting
```bash
conda create -n multiqc_env -c bioconda -y multiqc
```



3. Prepare input directories:
```bash
raw_reads/                                 # Raw FASTQ files (paired-end or single-end)
reference_genomes/MARV_downloads/          # Downloaded MARV genomes (unfiltered)
reference_genomes/MARV_compare/            # Filtered MARV reference genomes ≥18,000 bp
database/nextclade_marburg_dataset/        # Nextclade reference dataset (tree, reference.fasta, genome_annotation.gff3, etc.)
metadata/                                  # Sample metadata for TreeTime, IQ-TREE, and other analyses

```

## Usage
1. Place raw FASTQ files in a designated directory (e.g., raw_reads/).
2. Prepare the reference genomes in reference_genomes/MARV_downloads/ and filtered genomes in reference_genomes/MARV_compare/.
3. Prepare metadata for TreeTime/IQ-TREE in metadata/all_seq_metadata.csv.
4.  Run the full pipeline:
```bash
bash scripts/run_MARV-GEN_full_pipeline.sh
```
5. Alternatively, run individual steps as needed:
```bash
bash scripts/fastp_batch.sh
```
#### Step 1: FastQ Quality Control
```bash
bash scripts/fastp_batch.sh
```
#### Step 2: Host read removal
```bash
bash scripts/host_removal_batch.sh
```
#### Step 3: Mapping to reference genome
```bash
bash scripts/mapping_batch.sh
```
#### Step 4: BAM QC
```bash
bash scripts/qualimap_batch.sh
```
#### Step 5: Variant calling
```bash
bash scripts/variant_calling_batch.sh
```
#### Step 6: Consensus generation
```bash
bash scripts/consensus_batch.sh
```
#### Step 7: Coverage calculation
```bash
bash scripts/coverage_batch.sh
```
#### Step 8: Nextclade clade assignment
```bash
bash scripts/nextclade_batch.sh
```
#### Step 9: MultiQC reporting (Fastp, Qualimap, Nextclade)
```bash
bash scripts/multiqc_batch.sh
```
#### Step 10: Multiple sequence alignment (MAFFT)
bash scripts/msa_batch.sh
```
#### Step 11: Phylogenetic tree construction (IQ-TREE)
```bash
bash scripts/iqtree_batch.sh
```
#### Step 12: Time-resolved phylogeny (TreeTime)
```bash
bash scripts/treetime_batch.sh
```


## Directory Structure
```bash
MARV-GEN-VSP-Illumina-Analysis/
├── raw_reads/                                 # Raw FASTQ files (paired or single-end)
│   ├── MARV_X_1.fastq.gz
│   └── MARV_X_2.fastq.gz
├── reference_genomes/
│   ├── Marburg_reference.fasta                 # NC_001608.4.fasta
│   ├── Marburg_reference.gb                    # NC_001608.4.gb
│   ├── EF446131.1.fasta                        # outgroup.fasta
│   ├── MARV_downloads/                         # all Downloaded MARV genomes in fasta from NCBI
│   └── MARV_compare/                           # Filtered reference genomes (≥18,000 bp)
├── database/
│   └── nextclade_marburg_dataset/             # Nextclade reference dataset
│       ├── marburg_tree.json
│       ├── pathogen.json
│       ├── reference.fasta
│       ├── examples.fasta
│       ├── genome_annotation.gff3
│       ├── CHANGELOG.md
│       └── README.md
├── metadata/
│   └── all_seq_metadata.csv                    # Sample metadata for TreeTime/IQ-TREE
├── results/
│   ├── 01_fastp/                              # Fastp QC output
│   ├── 02_clean_reads/                         # Host-cleaned reads
│   ├── 03_nonhuman_reads/                      # Filtered non-human reads
│   ├── 04_mapped_bam/                          # BAM files from Minimap2
│   ├── 05_mapping_qc/                          # BAM QC (Qualimap)
│   ├── 06_variants/                            # Variant calling outputs
│   ├── 07_consensus/                           # Consensus sequences (iVar)
│   ├── 08_coverage/                            # Coverage statistics
│   ├── 09_nextclade/                           # Nextclade clade assignment results
│   ├── 10_msa/                                 # Multiple sequence alignment outputs (MAFFT)
│   ├── 11_phylogeny/                            # IQ-TREE phylogenetic trees
│   ├── 12_treetime/                             # Time-resolved phylogeny (TreeTime)
│       ├── visualization /
│       └── MARV.A.1/
│              └── visualization/
│   └── 13_multiqc/
│       ├── fastp/                              # MultiQC report for Fastp
│       ├── nextclade/                           # MultiQC report for Nextclade
│       ├── qualimap/                            # MultiQC report for Qualimap
│       └── combined/                            # MultiQC combined report
├── scripts/                                    # All pipeline scripts
├── logs/                                       # Pipeline logs
└── README.md
             
```
## SNPs analysis

The MARV-GEN repository includes scripts for advanced genomic analyses of the Marburg virus sequences, separate from the core Illumina processing pipeline.

### Lineage-Defining SNPs

**Script:** `marburg_lineage_snps.py`

- Identifies Ethiopian MARV lineage-defining SNPs.  
- Maps SNPs to coding sequences (CDS) and reports amino acid changes.  
- Outputs CSV/TSV tables for SNP positions and amino acid changes.

**Outputs:**

- `ethiopian_lineage_defining_snps_aa.csv`  
- `ethiopian_snps_protein.tsv`  

---

### Diversity Metrics & Selective Pressure

**Script:** `marburg_dn_ds_analysis.py`

- Calculates nucleotide diversity (π) and Watterson’s θ for Ethiopian sequences.  
- Computes observed nonsynonymous (N) and synonymous (S) substitutions.  
- Estimates dN/dS (ω) per gene to detect selective pressures.

**Outputs:**

- `selective_pressure_analysis.csv`  
- Console summary of genetic distance and diversity metrics.  

---

## Visualization

These scripts generate publication-quality visualizations to complement the analysis results.

### Genome Map of Lineage-Defining SNPs

**Script:** `marburg_lineage_snps.py` (genome map portion)

- Plots Ethiopian lineage-defining SNPs on the Marburg virus genome.  
- Highlights non-synonymous (red) vs. synonymous/non-coding (green) SNPs.  
- Shows coding sequences (CDS) as colored blocks with gene labels.  
- Produces high-resolution PNG figure for manuscripts.

**Output:**

- `ethiopian_snp_genome_map_aa_pro.png`  

---

### Time-Scaled Phylogenetic Tree

**Script:** `treetime_vis_2_MARV.A.1.py`

- Prunes MARV.A.1 tree to Ethiopian sequences.  
- Highlights Ethiopian sequences with larger markers.  
- Annotates bootstrap values and adds a temporal axis.  
- Produces PDF, PNG, SVG, and EPS high-resolution figures.

**Outputs:**

- `MARV.A1_Bootstrap_Labeled_NoTitle_Tree.{pdf,png,svg,eps}`  


## Logging
Each script writes per-sample logs in logs/.
A summary table is generated for mapping, variant calling, consensus, and coverage statistics.
Logs capture runtime, errors, and pipeline decisions.

## Authors

- **Betselot Zerihun Ayano** – [GitHub @betselotz](https://github.com/betselotz)  
- **Melak Getu Bire** – [GitHub @MelakG13](https://github.com/MelakG13)

## License

This repository is open for academic and research use.


