#!/usr/bin/env python3
import os
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---------------------------
# CONFIGURATION
# ---------------------------
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
ALIGNMENT_FILE = os.path.join(SCRIPT_DIR, "../results/10_msa/marburg_aligned.fasta")
GENBANK_FILE = os.path.join(SCRIPT_DIR, "../reference_genomes/Marburg_reference.gb")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/snp_difference")
CSV_FILE = os.path.join(OUTPUT_DIR, "ethiopian_lineage_defining_snps_aa.csv")
TSV_FILE = os.path.join(OUTPUT_DIR, "ethiopian_snps_protein.tsv")
NS_CSV_FILE = os.path.join(OUTPUT_DIR, "ethiopian_ns_snps.csv")
GENOME_MAP_FILE = os.path.join(OUTPUT_DIR, "ethiopian_snp_genome_map_aa.png")
ETHIOPIAN_PREFIX = "ET"
HIGHLIGHT_GENES = ["GP", "NP"]  # optional: genes to highlight in the genome map

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Colors
COLOR_NS = '#d62728'  # Non-synonymous
COLOR_S = '#2ca02c'   # Synonymous
CDS_COLOR = 'skyblue'
CDS_EDGE_COLOR = 'steelblue'

# ---------------------------
# READ ALIGNMENT & REFERENCE
# ---------------------------
alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
seq_ids = [rec.id for rec in alignment]
seqs = [str(rec.seq) for rec in alignment]

# Create DataFrame: each column = sequence, each row = nucleotide position
df = pd.DataFrame({seq_ids[i]: list(seqs[i]) for i in range(len(seqs))})

ethiopian_ids = [id for id in seq_ids if ETHIOPIAN_PREFIX in id]
reference_ids = [id for id in seq_ids if id not in ethiopian_ids]

# ---------------------------
# IDENTIFY LINEAGE-DEFINING SNPs
# ---------------------------
lineage_def_positions = []

for pos, row in df.iterrows():
    et_bases = set(row[ethiopian_ids])
    ref_bases = set(row[reference_ids])
    # Only one base among Ethiopian sequences and different from all reference sequences
    if len(et_bases) == 1 and et_bases.isdisjoint(ref_bases) and "-" not in et_bases:
        lineage_def_positions.append(pos)

print(f"Total lineage-defining SNPs identified: {len(lineage_def_positions)}")

# ---------------------------
# MAP SNPs TO CODONS AND AA CHANGES
# ---------------------------
record = SeqIO.read(GENBANK_FILE, "genbank")
gene_names, aa_changes, aa_colors = [], [], []

for pos in lineage_def_positions:
    gene_name, aa_change, color = "-", "-", COLOR_S
    for feature in record.features:
        if feature.type == "CDS" and feature.location.start <= pos <= feature.location.end:
            gene_name = feature.qualifiers.get("gene", feature.qualifiers.get("locus_tag", ["-"]))[0]
            cds_seq = feature.extract(record.seq)
            codon_start = (pos - feature.location.start) // 3 * 3
            codon_seq = cds_seq[codon_start:codon_start + 3]

            # Skip incomplete codons
            if len(codon_seq) == 3:
                aa_ref = str(codon_seq.translate())
                codon_list = list(codon_seq)
                codon_pos_in_codon = (pos - feature.location.start) % 3
                codon_list[codon_pos_in_codon] = df.iloc[pos][ethiopian_ids[0]].upper()
                aa_mut = str(Seq("".join(codon_list)).translate())
                aa_change = f"{aa_ref}{codon_start//3 + 1}{aa_mut}"
                # Only mark as non-synonymous if AA actually changes
                if aa_ref != aa_mut:
                    color = COLOR_NS
            else:
                print(f"Warning: Incomplete codon at position {pos+1}, skipping AA change")
            break
    gene_names.append(gene_name)
    aa_changes.append(aa_change)
    aa_colors.append(color)

# ---------------------------
# SAVE RESULTS
# ---------------------------
tsv_df = pd.DataFrame({
    "Position": [p+1 for p in lineage_def_positions],
    "Gene": gene_names,
    "Amino_Acid_Change": aa_changes,
    "Color": aa_colors
})

tsv_df.to_csv(CSV_FILE, index=False)
tsv_df.to_csv(TSV_FILE, sep='\t', index=False)
# Non-synonymous only
tsv_df[tsv_df['Color'] == COLOR_NS].to_csv(NS_CSV_FILE, index=False)

print(f"Non-synonymous SNPs: {sum(tsv_df['Color'] == COLOR_NS)}")
print(f"Synonymous SNPs: {sum(tsv_df['Color'] == COLOR_S)}")

# ---------------------------
# PLOT GENOME MAP WITH SNPs
# ---------------------------
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

plt.figure(figsize=(16, 6))
ax = plt.gca()
genome_length = len(record.seq)

# Genome backbone
ax.hlines(1, 0, genome_length, color='black', linewidth=1.5, zorder=1)

# CDS blocks with gene names
label_heights = [1.18, 1.35] 
cds_count = 0
for feature in record.features:
    if feature.type == "CDS":
        start, end = int(feature.location.start), int(feature.location.end)
        gene_name = feature.qualifiers.get("gene", feature.qualifiers.get("locus_tag", ["-"]))[0]

        # Always light blue for all genes
        facecolor = CDS_COLOR  

        ax.add_patch(plt.Rectangle((start, 0.92), end - start, 0.16, 
                                    facecolor=facecolor, alpha=0.6, edgecolor=CDS_EDGE_COLOR, lw=1, zorder=2))
        h = label_heights[cds_count % 2]
        ax.text((start+end)/2, h, gene_name, ha='center', va='bottom', fontsize=11, fontweight='bold', style='italic')
        cds_count += 1


# Plot SNPs
for idx, row in tsv_df.iterrows():
    pos = row['Position']
    if row['Color'] == COLOR_NS:
        ax.vlines(pos, 0.80, 1.20, color=COLOR_NS, linewidth=2, alpha=1.0, zorder=4)
        ax.text(pos, 0.75, row['Amino_Acid_Change'], color=COLOR_NS, 
                ha='center', va='top', fontsize=9, fontweight='bold', rotation=45)
    else:
        ax.vlines(pos, 0.88, 1.12, color=COLOR_S, linewidth=0.8, alpha=0.5, zorder=3)

# Legend
legend_elements = [
    Patch(facecolor=COLOR_NS, label='Non-synonymous (Red Lines)'),
    Patch(facecolor=COLOR_S, label='Synonymous (Green Lines)'),
    Patch(facecolor=CDS_COLOR, edgecolor=CDS_EDGE_COLOR, label='Coding Sequence (CDS)')
]
ax.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=3, frameon=False)

# Final polishing
ax.set_ylim(0.5, 1.6)
ax.set_xlim(0, genome_length)
for spine in ['top', 'right', 'left', 'bottom']: ax.spines[spine].set_visible(False)
ax.get_yaxis().set_visible(False)
ax.set_xlabel("Genome Position (nt)", fontsize=11)
ax.set_title("Distribution of Lineage-Defining SNPs", fontsize=14, pad=40, fontweight='bold')

plt.tight_layout()
plt.savefig(GENOME_MAP_FILE, dpi=600, bbox_inches='tight')
print(f"Files saved in {OUTPUT_DIR}")

