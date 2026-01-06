#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from math import log

# --- CONFIGURATION ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ALIGNMENT_FILE = os.path.join(SCRIPT_DIR, "../results/10_msa/marburg_aligned.fasta")
GENBANK_FILE = os.path.join(SCRIPT_DIR, "../reference_genomes/Marburg_reference.gb")
SNP_CSV = os.path.join(SCRIPT_DIR, "../results/snp_difference/ethiopian_lineage_defining_snps_aa.csv")
OUTPUT_SELECTIVE = os.path.join(SCRIPT_DIR, "../results/snp_difference/selective_pressure_analysis.csv")

# --- FUNCTIONS ---
def calculate_diversity(alignment, prefix="ET"):
    """Calculates Pi, Theta, and Segregating Sites for the prefix group."""
    ethiopian_seqs = [rec for rec in alignment if prefix in rec.id]
    n = len(ethiopian_seqs)
    L = alignment.get_alignment_length()
    if n < 2: 
        return 0, 0, 0
    
    matrix = np.array([list(str(rec.seq).upper()) for rec in ethiopian_seqs])
    
    # 1. Nucleotide Diversity (Pi)
    total_diffs, pairs = 0, 0
    for i in range(n):
        for j in range(i + 1, n):
            mask = (matrix[i] != '-') & (matrix[j] != '-') & (matrix[i] != 'N') & (matrix[j] != 'N')
            diffs = np.sum(matrix[i][mask] != matrix[j][mask])
            valid_sites = np.sum(mask)
            if valid_sites > 0:
                total_diffs += diffs / valid_sites
            pairs += 1
    pi = total_diffs / pairs if pairs > 0 else 0

    # 2. Watterson's Theta
    s_sites = 0
    for col in range(L):
        column_data = matrix[:, col]
        unique_bases = set(column_data)
        unique_bases.discard('-'); unique_bases.discard('N')
        if len(unique_bases) > 1: 
            s_sites += 1
    
    a_n = sum(1.0/i for i in range(1, n))
    theta_w = (s_sites / L) / a_n if a_n > 0 else 0
    return pi, theta_w, s_sites

def count_potential_sites(sequence):
    """Calculates potential S and N sites for dN/dS."""
    s_sites, n_sites = 0, 0
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) < 3: 
            continue
        if any(base not in "ACGT" for base in codon.upper()):
            continue
        for pos in range(3):
            orig_aa = codon.translate()
            for base in "ACGT":
                if base == codon[pos]: 
                    continue
                mut_codon = list(codon)
                mut_codon[pos] = base
                if Seq("".join(mut_codon)).translate() == orig_aa:
                    s_sites += 1/3
                else: 
                    n_sites += 1/3
    return s_sites, n_sites

def calculate_dn_ds(alignment, prefix="ET"):
    """Calculates observed N/S substitutions and dN/dS for each gene."""
    obs_N, obs_S = 0, 0
    seqs = [rec for rec in alignment if prefix in rec.id]
    if len(seqs) < 2: 
        return obs_N, obs_S, 0, 0, float('nan')
    
    L = alignment.get_alignment_length()
    seq1, seq2 = str(seqs[0].seq), str(seqs[1].seq)
    
    for k in range(0, L, 3):
        codon1, codon2 = seq1[k:k+3], seq2[k:k+3]
        if len(codon1) < 3 or len(codon2) < 3:
            continue
        # Skip invalid codons
        if any(base not in "ACGT" for base in codon1.upper()) or any(base not in "ACGT" for base in codon2.upper()):
            continue
        aa1, aa2 = Seq(codon1).translate(), Seq(codon2).translate()
        if aa1 != aa2:
            obs_N += 1
        elif codon1 != codon2:
            obs_S += 1

    # Rough dN/dS estimate (Nei & Gojobori method)
    pn = obs_N / max(1, (L/3))  # avoid division by zero
    ps = obs_S / max(1, (L/3))
    dn = -0.75 * log(1 - (4/3) * pn) if 0 < pn < 0.75 else pn
    ds = -0.75 * log(1 - (4/3) * ps) if 0 < ps < 0.75 else ps
    omega = dn/ds if ds > 0 else float('nan')
    
    return obs_N, obs_S, dn, ds, omega

# --- EXECUTION ---
print("Loading alignment...")
alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
gb_record = SeqIO.read(GENBANK_FILE, "genbank")
snp_df = pd.read_csv(SNP_CSV)

# 1. Diversity Metrics
print("Calculating Diversity Metrics...")
pi_val, theta_val, s_count = calculate_diversity(alignment, prefix="ET")

# 2. Genetic Distance (vs Reference)
ethiopian_ids = [rec.id for rec in alignment if "ET" in rec.id]
ref_ids = [rec.id for rec in alignment if "ET" not in rec.id]

def get_p_dist(id1, id2, align):
    s1 = next(rec.seq for rec in align if rec.id == id1)
    s2 = next(rec.seq for rec in align if rec.id == id2)
    mask = (np.array(list(s1)) != '-') & (np.array(list(s2)) != '-')
    diffs = np.sum(np.array(list(s1))[mask] != np.array(list(s2))[mask])
    return diffs / np.sum(mask)

avg_dist = np.mean([get_p_dist(et, ref, alignment) for et in ethiopian_ids for ref in ref_ids])

# 3. dN/dS Calculation per gene
results = []
for feature in gb_record.features:
    if feature.type == "CDS":
        gene = feature.qualifiers.get("gene", ["Unknown"])[0]
        cds_seq = feature.extract(gb_record.seq)
        S_pot, N_pot = count_potential_sites(cds_seq)
        gene_snps = snp_df[snp_df['Gene'] == gene]
        Obs_N = len(gene_snps[gene_snps['AA_Color'] == 'red'])
        Obs_S = len(gene_snps[gene_snps['AA_Color'] == 'green'])
        pn = Obs_N / S_pot if S_pot > 0 else 0
        ps = Obs_S / N_pot if N_pot > 0 else 0
        dn = -0.75 * log(1 - (4/3) * pn) if 0 < pn < 0.75 else pn
        ds = -0.75 * log(1 - (4/3) * ps) if 0 < ps < 0.75 else ps
        omega = dn / ds if ds > 0 else float('nan')
        results.append({"Gene": gene, "Obs_N": Obs_N, "Obs_S": Obs_S, "dN": round(dn, 6), "dS": round(ds, 6), "omega": round(omega, 4)})

# --- OUTPUT ---
print(f"\n--- FINAL RESULTS ---")
print(f"Genetic Distance (ET vs Others): {avg_dist:.4%}")
print(f"Ethiopian Nucleotide Diversity (Pi): {pi_val:.6f}")
print(f"Ethiopian Watterson's Theta: {theta_val:.6f}")
print(f"Segregating Sites in ET: {s_count}")
print("\nSelective Pressure (dN/dS):")
print(pd.DataFrame(results))
pd.DataFrame(results).to_csv(OUTPUT_SELECTIVE, index=False)

