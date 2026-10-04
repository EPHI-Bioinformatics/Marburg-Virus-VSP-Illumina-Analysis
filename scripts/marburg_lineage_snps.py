#!/usr/bin/env python3
"""
Marburg Virus Ethiopian Clade Lineage-Defining SNP Analysis
Version: 4.2 (strict consensus + codon-level annotation + clean 3-row legend)

Identify lineage-defining mutations on the branch leading to the Ethiopian
Marburg virus outbreak using a strict consensus comparison.
"""

import os
import sys
import re
import csv
import logging
from collections import Counter
from math import log

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# ---------------------------------------------------------------------------
# CONFIGURATION – paths for your machine
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

GB_PATH      = "/media/betselotz/Expansion/MARV-Gen/reference_genomes/EF446132.gb"
ALN_PATH     = "/media/betselotz/Expansion/MARV-Gen/results/10_msa/MARV_aligned.fasta"
LINEAGE_FILE = "/media/betselotz/Expansion/MARV-Gen/reference_genomes/lineage_assignments.csv"  # optional
OUTPUT_DIR   = "/media/betselotz/Expansion/MARV-Gen/results/snp_difference"

REF_HEADER_MATCH = "EF446132"
GROUP_KEYWORD    = "ET_MARV"
A1_LINEAGE_LABEL = "MARV.A.1"

VALID_BASES = {"A", "C", "G", "T"}

os.makedirs(OUTPUT_DIR, exist_ok=True)

PRIMARY_CSV      = os.path.join(OUTPUT_DIR, "lineage_defining_mutations.csv")
PRIMARY_XLSX     = os.path.join(OUTPUT_DIR, "lineage_defining_mutations.xlsx")
PER_SAMPLE_CSV   = os.path.join(OUTPUT_DIR, "ethiopian_snps_per_sample.csv")
NS_CSV           = os.path.join(OUTPUT_DIR, "ethiopian_ns_snps.csv")
EXTENDED_TABLE3  = os.path.join(OUTPUT_DIR, "extended_table_3_ns_snps.csv")
GENE_SUMMARY_CSV = os.path.join(OUTPUT_DIR, "gene_summary.csv")
SELECTIVE_CSV    = os.path.join(OUTPUT_DIR, "selective_pressure_analysis.csv")
SUPPORT_CSV      = os.path.join(OUTPUT_DIR, "snp_support_summary.csv")
SUMMARY_TXT      = os.path.join(OUTPUT_DIR, "summary_statistics.txt")
GENOME_MAP_BASE  = os.path.join(OUTPUT_DIR, "ethiopian_snp_genome_map_aa")

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# ---------------------------------------------------------------------------
# 1. Parse FASTA
# ---------------------------------------------------------------------------
def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = []
            else:
                seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}

seqs = read_fasta(ALN_PATH)
aln_len = len(next(iter(seqs.values())))
assert all(len(v) == aln_len for v in seqs.values()), "Unequal sequence lengths"

ref_name = next((k for k in seqs if REF_HEADER_MATCH in k), None)
if ref_name is None:
    print("WARNING: No sequence containing 'EF446132' found. Headers:")
    for h in seqs:
        print(" ", h)
    sys.exit("ERROR: Set REF_HEADER_MATCH correctly.")
ref_aln = seqs[ref_name]
print(f"Reference sequence: {ref_name}")

# ---------------------------------------------------------------------------
# 2. Groups
# ---------------------------------------------------------------------------
eth_names = [k for k in seqs if GROUP_KEYWORD in k]
non_eth_all = [k for k in seqs if GROUP_KEYWORD not in k]

if os.path.exists(LINEAGE_FILE):
    try:
        lineage_df = pd.read_csv(LINEAGE_FILE, dtype=str)
        lineage_df.columns = [c.strip().lower() for c in lineage_df.columns]
        if {"sequence_id", "lineage"}.issubset(lineage_df.columns):
            lineage_map = dict(zip(lineage_df["sequence_id"], lineage_df["lineage"]))
            non_eth_names = [i for i in non_eth_all if lineage_map.get(i) == A1_LINEAGE_LABEL]
            if not non_eth_names:
                logging.warning(f"No sequences matched lineage '{A1_LINEAGE_LABEL}'. Using all non-Ethiopian.")
                non_eth_names = non_eth_all
        else:
            non_eth_names = non_eth_all
    except Exception:
        non_eth_names = non_eth_all
else:
    non_eth_names = non_eth_all

print(f"Ethiopia group (n={len(eth_names)}):")
for n in eth_names:
    print(f"  {n}")
print(f"Non-Ethiopia background (n={len(non_eth_names)}):")
for n in non_eth_names:
    print(f"  {n}")

if len(eth_names) == 0 or len(non_eth_names) == 0:
    sys.exit("ERROR: One of the two groups is empty.")

# ---------------------------------------------------------------------------
# 3. Coordinate mapping
# ---------------------------------------------------------------------------
ref_to_aln = {}
pos = 0
for col, ch in enumerate(ref_aln, start=1):
    if ch != "-":
        pos += 1
        ref_to_aln[pos] = col
aln_to_ref = {v: k for k, v in ref_to_aln.items()}
ref_len = pos
print(f"Reference ungapped length: {ref_len} bp")

# ---------------------------------------------------------------------------
# 4. CDS features
# ---------------------------------------------------------------------------
def parse_genbank_cds(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    lines = text.split("\n")
    cds_list = []
    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\r")
        m = re.match(r"^\s{5}CDS\s+(\S+)", line)
        if m:
            loc_str = m.group(1)
            loc_match = re.search(r"(\d+)\.\.(\d+)", loc_str)
            if not loc_match:
                i += 1
                continue
            start, end = int(loc_match.group(1)), int(loc_match.group(2))
            gene = product = None
            codon_start = 1
            j = i + 1
            qual_buffer = ""
            while j < len(lines):
                qline = lines[j].rstrip("\r")
                if re.match(r"^\s{5}\S", qline):
                    break
                qual_buffer += qline.strip() + " "
                j += 1
            gm = re.search(r'/gene="([^"]+)"', qual_buffer)
            pm = re.search(r'/product="([^"]+)"', qual_buffer)
            cm = re.search(r"/codon_start=(\d+)", qual_buffer)
            if gm: gene = gm.group(1)
            if pm: product = pm.group(1)
            if cm: codon_start = int(cm.group(1))
            cds_list.append({
                "gene": gene or "unknown",
                "product": product or "",
                "start": start,
                "end": end,
                "codon_start": codon_start,
            })
            i = j
        else:
            i += 1
    return cds_list

cds_features = parse_genbank_cds(GB_PATH)
print(f"\nParsed {len(cds_features)} CDS features:")
for c in cds_features:
    print(f"  {c['gene']:6s} {c['start']:6d}-{c['end']:6d}  {c['product']}")

def genome_pos_to_cds(genome_pos):
    for c in cds_features:
        if c["start"] <= genome_pos <= c["end"]:
            return c
    return None

# ---------------------------------------------------------------------------
# 5. Strict consensus
# ---------------------------------------------------------------------------
def strict_consensus_base(names, col):
    bases = set()
    for n in names:
        b = seqs[n][col - 1]
        if b not in VALID_BASES:
            return None
        bases.add(b)
    if len(bases) == 1:
        return bases.pop()
    return None

eth_consensus = [strict_consensus_base(eth_names, col) for col in range(1, aln_len + 1)]
non_eth_consensus = [strict_consensus_base(non_eth_names, col) for col in range(1, aln_len + 1)]

# ---------------------------------------------------------------------------
# 6. Differing columns
# ---------------------------------------------------------------------------
diff_columns = []
for col in range(1, aln_len + 1):
    a = non_eth_consensus[col - 1]
    b = eth_consensus[col - 1]
    if a is not None and b is not None and a != b and col in aln_to_ref:
        diff_columns.append(col)

print(f"\nTotal lineage-defining differences (mapped to reference): {len(diff_columns)}")

# ---------------------------------------------------------------------------
# 7. Annotate
# ---------------------------------------------------------------------------
def get_codon_from_consensus(consensus_list, ref_positions, fallback_seq=None):
    codon = ""
    for rp in ref_positions:
        col = ref_to_aln.get(rp)
        if col is None:
            return None
        base = consensus_list[col - 1]
        if base is None:
            if fallback_seq is not None:
                base = fallback_seq[col - 1]
                if base not in VALID_BASES:
                    return None
            else:
                return None
        codon += base
    return codon

noncoding_cols = []
cds_codon_keys = {}

for col in diff_columns:
    genome_pos = aln_to_ref[col]
    cds = genome_pos_to_cds(genome_pos)
    if cds is None:
        noncoding_cols.append(col)
        continue
    nt_in_cds = genome_pos - cds["start"] + 1
    codon_num = (nt_in_cds - 1) // 3 + 1
    cds_codon_keys[(cds["gene"], codon_num)] = cds

results = []

for col in noncoding_cols:
    genome_pos = aln_to_ref[col]
    a_base = non_eth_consensus[col - 1]
    b_base = eth_consensus[col - 1]
    results.append({
        "CDS_gene": "non-CDS (intergenic / UTR)",
        "CDS_product": None,
        "nt_position_in_CDS": None,
        "codon_position_(aa)": None,
        "ref_genome_position_(nt)": str(genome_pos),
        "alignment_position_(col)": str(col),
        "codon_MARV.A.1": a_base,
        "codon_Ethiopia.2025": b_base,
        "aa_MARV.A.1": None,
        "aa_Ethiopia.2025": None,
        "aa_mutation": "n/a (non-coding)",
        "effect": "non-coding",
        "Ancestral_A1_Base": a_base,
        "Ethiopian_Base": b_base,
        "ET_Support": sum(1 for n in eth_names if seqs[n][col-1] == b_base),
        "ET_N": len(eth_names),
        "A1_Support": sum(1 for n in non_eth_names if seqs[n][col-1] == a_base),
        "A1_N": len(non_eth_names),
    })

for (gene, codon_num), cds in cds_codon_keys.items():
    codon_first_nt_in_cds = (codon_num - 1) * 3 + 1
    codon_first_genome_pos = cds["start"] + codon_first_nt_in_cds - 1
    ref_positions = [codon_first_genome_pos,
                     codon_first_genome_pos + 1,
                     codon_first_genome_pos + 2]
    nt_in_cds = codon_first_nt_in_cds

    codon_a = get_codon_from_consensus(non_eth_consensus, ref_positions, fallback_seq=ref_aln)
    codon_b = get_codon_from_consensus(eth_consensus, ref_positions,
                                       fallback_seq=seqs[eth_names[0]])

    aa_a = CODON_TABLE.get(codon_a, "?") if codon_a else None
    aa_b = CODON_TABLE.get(codon_b, "?") if codon_b else None

    if aa_a and aa_b:
        aa_mut = f"{aa_a}{codon_num}{aa_b}"
        effect = "synonymous" if aa_a == aa_b else "non-synonymous"
    else:
        aa_mut = "n/a (incomplete codon consensus)"
        effect = "ambiguous"

    first_diff_col = None
    for rp in ref_positions:
        c = ref_to_aln.get(rp)
        if c and non_eth_consensus[c-1] != eth_consensus[c-1]:
            first_diff_col = c
            break
    if first_diff_col is None:
        first_diff_col = ref_to_aln[ref_positions[0]]

    a_base = non_eth_consensus[first_diff_col - 1]
    b_base = eth_consensus[first_diff_col - 1]

    results.append({
        "CDS_gene": gene,
        "CDS_product": cds["product"],
        "nt_position_in_CDS": str(nt_in_cds),
        "codon_position_(aa)": str(codon_num),
        "ref_genome_position_(nt)": ";".join(str(p) for p in ref_positions),
        "alignment_position_(col)": ";".join(str(ref_to_aln[p]) for p in ref_positions),
        "codon_MARV.A.1": codon_a,
        "codon_Ethiopia.2025": codon_b,
        "aa_MARV.A.1": aa_a,
        "aa_Ethiopia.2025": aa_b,
        "aa_mutation": aa_mut,
        "effect": effect,
        "Ancestral_A1_Base": a_base,
        "Ethiopian_Base": b_base,
        "ET_Support": sum(1 for n in eth_names if seqs[n][first_diff_col-1] == b_base),
        "ET_N": len(eth_names),
        "A1_Support": sum(1 for n in non_eth_names if seqs[n][first_diff_col-1] == a_base),
        "A1_N": len(non_eth_names),
    })

results.sort(key=lambda r: int(str(r["ref_genome_position_(nt)"]).split(";")[0]))

# ---------------------------------------------------------------------------
# 8. Write primary table
# ---------------------------------------------------------------------------
fieldnames = [
    "CDS_gene", "CDS_product", "nt_position_in_CDS", "codon_position_(aa)",
    "ref_genome_position_(nt)", "alignment_position_(col)",
    "codon_MARV.A.1", "codon_Ethiopia.2025", "aa_MARV.A.1", "aa_Ethiopia.2025",
    "aa_mutation", "effect",
    "Ancestral_A1_Base", "Ethiopian_Base", "ET_Support", "ET_N", "A1_Support", "A1_N",
]

with open(PRIMARY_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(results)

try:
    pd.DataFrame(results).to_excel(PRIMARY_XLSX, index=False)
except Exception:
    pass

# ---------------------------------------------------------------------------
# 9. Nucleotide-level tables
# ---------------------------------------------------------------------------
per_sample_rows = []
for col in diff_columns:
    genome_pos = aln_to_ref[col]
    a_base = non_eth_consensus[col - 1]
    b_base = eth_consensus[col - 1]
    cds = genome_pos_to_cds(genome_pos)
    if cds:
        nt_in_cds = genome_pos - cds["start"] + 1
        codon_num = (nt_in_cds - 1) // 3 + 1
        codon_start_genome = cds["start"] + ((nt_in_cds - 1) // 3) * 3
        ref_pos_codon = [codon_start_genome, codon_start_genome+1, codon_start_genome+2]
        codon_a = get_codon_from_consensus(non_eth_consensus, ref_pos_codon, fallback_seq=ref_aln)
        codon_b = get_codon_from_consensus(eth_consensus, ref_pos_codon, fallback_seq=seqs[eth_names[0]])
        aa_a = CODON_TABLE.get(codon_a, "?") if codon_a else "?"
        aa_b = CODON_TABLE.get(codon_b, "?") if codon_b else "?"
        aa_change = f"{aa_a}{codon_num}{aa_b}" if aa_a != "?" and aa_b != "?" else "-"
        gene = cds["gene"]
        is_ns = aa_a != aa_b and aa_a != "?" and aa_b != "?"
        typ = "NS" if is_ns else "S"
    else:
        gene = "Intergenic"
        aa_change = "-"
        typ = "S"

    entry = {
        "Position": genome_pos,
        "Gene": gene,
        "Amino_Acid_Change": aa_change,
        "Type": typ,
        "Ancestral_A1_Base": a_base,
        "Ethiopian_Base": b_base,
        "ET_Support": sum(1 for n in eth_names if seqs[n][col-1] == b_base),
        "ET_N": len(eth_names),
        "A1_Support": sum(1 for n in non_eth_names if seqs[n][col-1] == a_base),
        "A1_N": len(non_eth_names),
    }
    for et_id in eth_names:
        entry[et_id] = f"{a_base}→{seqs[et_id][col-1]}"
    per_sample_rows.append(entry)

shared_df = pd.DataFrame(per_sample_rows)
shared_df.to_csv(PER_SAMPLE_CSV, index=False)

ns_df = shared_df[shared_df["Type"] == "NS"].copy()
ns_df.to_csv(NS_CSV, index=False)

if not ns_df.empty:
    ns_table = ns_df[["Position", "Gene", "Amino_Acid_Change",
                      "Ancestral_A1_Base", "Ethiopian_Base",
                      "ET_Support", "ET_N", "A1_Support", "A1_N"]].copy()
else:
    ns_table = pd.DataFrame(columns=["Position", "Gene", "Amino_Acid_Change",
                                     "Ancestral_A1_Base", "Ethiopian_Base",
                                     "ET_Support", "ET_N", "A1_Support", "A1_N"])
ns_table.to_csv(EXTENDED_TABLE3, index=False)

if not shared_df.empty:
    gene_summary = (shared_df.groupby("Gene")
                    .agg(Total_SNPs=("Position", "count"),
                         NS_SNPs=("Type", lambda x: (x == "NS").sum()),
                         S_SNPs=("Type", lambda x: (x == "S").sum()))
                    .reset_index())
else:
    gene_summary = pd.DataFrame(columns=["Gene", "Total_SNPs", "NS_SNPs", "S_SNPs"])
gene_summary.to_csv(GENE_SUMMARY_CSV, index=False)

support_cols = ["Position", "Gene", "Amino_Acid_Change", "Type",
                "ET_Support", "ET_N", "A1_Support", "A1_N"]
shared_df[support_cols].to_csv(SUPPORT_CSV, index=False)

# ---------------------------------------------------------------------------
# 10. Selective pressure
# ---------------------------------------------------------------------------
sel_results = []
if HAS_BIOPYTHON:
    try:
        record = SeqIO.read(GB_PATH, "genbank")
        def count_potential_sites(sequence):
            s_sites = n_sites = 0.0
            seq_str = str(sequence)
            for i in range(0, len(seq_str) - 2, 3):
                codon = seq_str[i:i+3]
                if len(codon) != 3 or any(b not in "ACGT" for b in codon):
                    continue
                orig_aa = str(Seq(codon).translate())
                for pos_in_codon in range(3):
                    orig_base = codon[pos_in_codon]
                    for base in "ACGT":
                        if base == orig_base:
                            continue
                        mut_list = list(codon)
                        mut_list[pos_in_codon] = base
                        mut_aa = str(Seq("".join(mut_list)).translate())
                        if mut_aa == orig_aa:
                            s_sites += 1/3
                        else:
                            n_sites += 1/3
            return s_sites, n_sites

        for feature in record.features:
            if feature.type != "CDS":
                continue
            gene = feature.qualifiers.get("gene", ["Unknown"])[0]
            cds_seq = feature.extract(record.seq)
            S_pot, N_pot = count_potential_sites(cds_seq)
            gene_snps = shared_df[shared_df["Gene"] == gene] if not shared_df.empty else shared_df
            Obs_N = len(gene_snps[gene_snps["Type"] == "NS"]) if not shared_df.empty else 0
            Obs_S = len(gene_snps[gene_snps["Type"] == "S"]) if not shared_df.empty else 0
            pn = Obs_N / N_pot if N_pot > 0 else 0
            ps = Obs_S / S_pot if S_pot > 0 else 0
            dn = -0.75 * log(1 - (4/3)*pn) if 0 < pn < 0.75 else pn
            ds = -0.75 * log(1 - (4/3)*ps) if 0 < ps < 0.75 else ps
            omega = dn / ds if ds > 0 else float("nan")
            sel_results.append({
                "Gene": gene,
                "Obs_NS_SNPs": Obs_N,
                "Obs_S_SNPs": Obs_S,
                "dN": round(dn, 6),
                "dS": round(ds, 6),
                "Approx_Selection_Index": round(omega, 4),
            })
    except Exception as e:
        logging.warning(f"Selective-pressure calculation skipped: {e}")

pd.DataFrame(sel_results).to_csv(SELECTIVE_CSV, index=False)

# ---------------------------------------------------------------------------
# 11. Summary
# ---------------------------------------------------------------------------
effect_counts = Counter(r["effect"] for r in results)
gene_counts = Counter(r["CDS_gene"] for r in results)
n_cds = sum(v for k, v in gene_counts.items() if k != "non-CDS (intergenic / UTR)")
n_nc = gene_counts.get("non-CDS (intergenic / UTR)", 0)

with open(SUMMARY_TXT, "w") as f:
    f.write(f"Total lineage-defining mutations (codon-level rows): {len(results)}\n")
    f.write(f"  In CDS: {n_cds}\n")
    f.write(f"  Non-coding: {n_nc}\n")
    f.write(f"Effects: {dict(effect_counts)}\n")
    f.write(f"By gene: {dict(gene_counts)}\n")
    f.write(f"\nNucleotide-level view:\n")
    f.write(f"  Total branch-defining SNPs: {len(shared_df)}\n")
    f.write(f"  Non-synonymous SNPs: {len(ns_df)}\n")
    f.write(f"Ethiopian samples: {len(eth_names)}\n")
    f.write(f"Non-Ethiopian MARV.A.1 background samples: {len(non_eth_names)}\n")
    f.write("Method: strict consensus (all sequences must agree and be unambiguous A/C/G/T)\n")

print("\n=== SUMMARY ===")
print(f"Total lineage-defining mutations (codon-level): {len(results)}")
print(f"  In CDS: {n_cds}")
print(f"  Non-coding: {n_nc}")
print(f"Effects: {dict(effect_counts)}")
print(f"By gene: {dict(gene_counts)}")
print(f"Nucleotide-level NS SNPs: {len(ns_df)}")

# ---------------------------------------------------------------------------
# 12. Genome map – clean horizontal markers + 3-row legend
# ---------------------------------------------------------------------------
def plot_snp_genome_map(df, base_filename, genome_length, cds_features):
    if df.empty:
        print("No branch-defining SNPs – skipping genome map.")
        return

    fig, ax = plt.subplots(figsize=(24, 12))

    # Genome backbone
    ax.hlines(1.0, 0, genome_length, color="black", linewidth=1.8, zorder=1)

    # CDS regions
    for c in cds_features:
        start, end = c["start"], c["end"]
        ax.add_patch(plt.Rectangle(
            (start, 0.92), end - start, 0.16,
            facecolor="skyblue", alpha=0.6,
            edgecolor="steelblue", lw=1, zorder=2
        ))
        ax.text((start + end) / 2, 1.32, c["gene"],
                ha="center", va="bottom",
                fontsize=17, fontweight="bold", style="italic",
                color="#333333", zorder=5)

    # Separate NS and S
    ns_df = df[df["Type"] == "NS"].sort_values("Position").reset_index(drop=True)
    s_df  = df[df["Type"] == "S"]

    # Synonymous – short green ticks
    for _, row in s_df.iterrows():
        pos = row["Position"]
        ax.vlines(pos, 0.90, 1.10, color="#006400", linewidth=1.4, alpha=0.9, zorder=3)

    # Non-synonymous markers
    available_markers = ["*", "o", "s", "D", "^", "v", ">", "<", "p", "X", "h", "8", "P", "d"]
    colormap = matplotlib.colormaps["tab20"].resampled(40)

    marker_map = {}
    for i, (_, row) in enumerate(ns_df.iterrows()):
        aa = row["Amino_Acid_Change"]
        if aa not in marker_map:
            idx = len(marker_map)
            marker_map[aa] = {
                "marker": available_markers[idx % len(available_markers)],
                "color": colormap(idx),
            }

    # Three clean horizontal rows for markers
    levels = [0.68, 0.48, 0.28]

    for i, (_, row) in enumerate(ns_df.iterrows()):
        pos = row["Position"]
        level = levels[i % len(levels)]
        m = marker_map[row["Amino_Acid_Change"]]

        # connector line
        ax.plot([pos, pos], [0.85, level + 0.04],
                color="gray", linestyle="--", linewidth=0.7, zorder=3)

        # red highlight on backbone
        ax.vlines(pos, 0.85, 1.15, color="#d62728", linewidth=1.6, zorder=4)

        # marker
        ax.scatter(pos, level,
                   marker=m["marker"],
                   color=m["color"],
                   s=130,
                   zorder=6,
                   edgecolors="black",
                   linewidths=0.6)

    # ---------- Legend forced into exactly 3 rows ----------
    legend_elements = [
        Patch(facecolor="skyblue", edgecolor="steelblue", label="CDS Region"),
        Line2D([0], [0], color="#006400", lw=4, label="Synonymous"),
        Line2D([0], [0], color="#d62728", lw=3, label="Non-synonymous"),
    ]
    for aa, m in marker_map.items():
        legend_elements.append(
            Line2D([0], [0],
                   marker=m["marker"], color="w",
                   markerfacecolor=m["color"],
                   markeredgecolor="black",
                   markersize=11,
                   label=f"NS: {aa}",
                   linestyle="None")
        )

    n_items = len(legend_elements)
    ncol = (n_items + 2) // 3          # force 3 rows

    ax.legend(handles=legend_elements,
              loc="upper left",
              bbox_to_anchor=(0.0, -0.13),
              ncol=ncol,
              frameon=True,
              fontsize=10,
              handletextpad=0.4,
              columnspacing=1.1,
              borderpad=0.6)

    ax.set_ylim(0.05, 1.55)
    ax.set_xlim(0, genome_length)
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel("Genome Position (nt)", fontsize=16, fontweight="bold")
    ax.set_title(
        "Mutations on the Branch Leading to the Ethiopian Outbreak",
        fontsize=18, pad=25, fontweight="bold"
    )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    for fmt in [".png", ".svg", ".pdf", ".eps"]:
        for dpi in (300, 200, 150):
            try:
                plt.savefig(base_filename + fmt, dpi=dpi, bbox_inches="tight")
                print(f"Saved: {base_filename}{fmt} (dpi={dpi})")
                break
            except RuntimeError:
                continue
    plt.close()

plot_snp_genome_map(shared_df, GENOME_MAP_BASE, ref_len, cds_features)

# ---------------------------------------------------------------------------
# 13. Ethiopian intra-group pairwise SNP differences
# ---------------------------------------------------------------------------
PAIRWISE_CSV        = os.path.join(OUTPUT_DIR, "ethiopian_pairwise_snp_differences.csv")
PAIRWISE_MATRIX_CSV = os.path.join(OUTPUT_DIR, "ethiopian_pairwise_snp_matrix.csv")

def pairwise_snp_diff(name_a, name_b):
    """Return list of (col, genome_pos, base_a, base_b) where two Ethiopian
    sequences differ, restricted to unambiguous A/C/G/T calls (gaps/N excluded)."""
    seq_a, seq_b = seqs[name_a], seqs[name_b]
    diffs = []
    for col in range(1, aln_len + 1):
        ba, bb = seq_a[col - 1], seq_b[col - 1]
        if ba in VALID_BASES and bb in VALID_BASES and ba != bb:
            diffs.append((col, aln_to_ref.get(col), ba, bb))
    return diffs

pairwise_rows = []
matrix = {a: {b: 0 for b in eth_names} for a in eth_names}

for i, a in enumerate(eth_names):
    for b in eth_names[i + 1:]:
        diffs = pairwise_snp_diff(a, b)
        matrix[a][b] = matrix[b][a] = len(diffs)
        for col, genome_pos, ba, bb in diffs:
            cds = genome_pos_to_cds(genome_pos) if genome_pos else None
            pairwise_rows.append({
                "Sample_A": a,
                "Sample_B": b,
                "ref_genome_position_(nt)": genome_pos,
                "alignment_position_(col)": col,
                "Gene": cds["gene"] if cds else "non-CDS (intergenic / UTR)",
                "Base_A": ba,
                "Base_B": bb,
            })

pairwise_df = pd.DataFrame(pairwise_rows).sort_values(
    ["Sample_A", "Sample_B", "ref_genome_position_(nt)"]
) if pairwise_rows else pd.DataFrame(columns=[
    "Sample_A", "Sample_B", "ref_genome_position_(nt)",
    "alignment_position_(col)", "Gene", "Base_A", "Base_B"
])
pairwise_df.to_csv(PAIRWISE_CSV, index=False)

matrix_df = pd.DataFrame(matrix).reindex(index=eth_names, columns=eth_names)
matrix_df.to_csv(PAIRWISE_MATRIX_CSV)

print("\n=== ETHIOPIAN INTRA-GROUP PAIRWISE SNP DIFFERENCES ===")
for i, a in enumerate(eth_names):
    for b in eth_names[i + 1:]:
        print(f"  {a} vs {b}: {matrix[a][b]} SNP(s)")
print(f"Pairwise SNP list   : {PAIRWISE_CSV}")
print(f"Pairwise SNP matrix : {PAIRWISE_MATRIX_CSV}")

print(f"\n=== ANALYSIS COMPLETE ===")
print(f"Primary codon-level table : {PRIMARY_CSV}")
print(f"Nucleotide-level tables   : {PER_SAMPLE_CSV}, {NS_CSV}")
print(f"Genome map                : {GENOME_MAP_BASE}.*")
print(f"All outputs written to    : {OUTPUT_DIR}")
