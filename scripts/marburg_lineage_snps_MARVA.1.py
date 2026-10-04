#!/usr/bin/env python3
import os
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from itertools import combinations

# ============================
# CONFIG
# ============================
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
ALIGNMENT_FILE = os.path.join(SCRIPT_DIR, "../results/10_msa/MARV.A.1/marburg_aligned.fasta")
GENBANK_FILE = os.path.join(SCRIPT_DIR, "../reference_genomes/Marburg_reference.gb")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/core_genome")
os.makedirs(OUTPUT_DIR, exist_ok=True)

ETHIOPIAN_IDS = ["ET_MARV_23","ET_MARV_31","ET_MARV_32","ET_MARV_261","ET_MARV_262","ET_MARV_45","ET_MARV_60"]
REFERENCE_IDS = ["JN408064.1","JX458851.1","KC545387.1","KC545388.1","JX458853.1","JX458858.1","EF446132.1"]

VALID = {"A","C","G","T","a","c","g","t"}

# ============================
# LOAD DATA
# ============================
alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
df = pd.DataFrame({r.id: list(str(r.seq)) for r in alignment})
record = SeqIO.read(GENBANK_FILE, "genbank")

# ============================
# CORE GENOME MASK
# ============================
core_positions = []
for pos, row in df.iterrows():
    if all(row[s] in VALID for s in ETHIOPIAN_IDS):
        core_positions.append(pos)

core_pct = 100*len(core_positions)/len(df)
print(f"Core genome: {len(core_positions)} bp ({core_pct:.2f}%)")

# ============================
# EXPORT CORE ALIGNMENT
# ============================
core_fasta = []
for sid in df.columns:
    seq = "".join(df.loc[p,sid] for p in core_positions)
    core_fasta.append(f">{sid}\n{seq}")

with open(OUTPUT_DIR+"/core_alignment.fasta","w") as f:
    f.write("\n".join(core_fasta))

# ============================
# GLOBAL CONSENSUS (core only)
# ============================
ref_consensus = df.loc[core_positions,REFERENCE_IDS].mode(axis=1)[0]

# ============================
# CODON EFFECT
# ============================
def get_effect(feature, pos, alt):
    cds = list(feature.location)
    if pos not in cds:
        return None
    local = cds.index(pos)
    seq = feature.extract(record.seq)
    codon_start = (local//3)*3
    codon = seq[codon_start:codon_start+3]
    if len(codon)<3: return None
    aa_ref = codon.translate()
    codon = list(str(codon))
    codon[local%3] = alt
    aa_alt = Seq("".join(codon)).translate()
    return aa_ref, aa_alt

# ============================
# CORE GENE LENGTHS
# ============================
gene_core = {}
for f in record.features:
    if f.type=="CDS":
        g = f.qualifiers.get("gene",["?"])[0]
        gene_core[g]=sum(1 for p in f.location if p in core_positions)//3

# ============================
# SNP CALLING
# ============================
shared, private = [], []
dn_ds = {g:{"N":0,"S":0} for g in gene_core}

for pos in core_positions:
    row = df.loc[pos]
    ref = ref_consensus[pos]

    et = row[ETHIOPIAN_IDS]
    et_unique = set(et)

    # Shared Ethiopian
    if len(et_unique)==1 and list(et_unique)[0]!=ref:
        alt=list(et_unique)[0]
        for f in record.features:
            if f.type=="CDS" and pos in f.location:
                g=f.qualifiers.get("gene",["?"])[0]
                eff=get_effect(f,pos,alt)
                if eff:
                    a,b=eff
                    t="NS" if a!=b else "S"
                    if t=="NS": dn_ds[g]["N"]+=1
                    else: dn_ds[g]["S"]+=1
                    shared.append([pos+1,g,ref,alt,str(a)+str((pos+1))+str(b),t])

    # Private
    else:
        for s in ETHIOPIAN_IDS:
            alt=row[s]
            if alt==ref: continue
            if list(et).count(alt)==1:
                for f in record.features:
                    if f.type=="CDS" and pos in f.location:
                        g=f.qualifiers.get("gene",["?"])[0]
                        eff=get_effect(f,pos,alt)
                        if eff:
                            a,b=eff
                            t="NS" if a!=b else "S"
                            private.append([pos+1,s,g,ref,alt,str(a)+str((pos+1))+str(b),t])

pd.DataFrame(shared,columns=["Pos","Gene","Ref","Alt","AA","Type"]).to_csv(OUTPUT_DIR+"/shared_snps.csv",index=False)
pd.DataFrame(private,columns=["Pos","Sample","Gene","Ref","Alt","AA","Type"]).to_csv(OUTPUT_DIR+"/private_snps.csv",index=False)

# ============================
# dN/dS
# ============================
dnds=[]
for g in gene_core:
    N=dn_ds[g]["N"]
    S=dn_ds[g]["S"]
    L=gene_core[g]
    dnds.append([g,N,S,(N/L)/(S/L+1e-6)])

pd.DataFrame(dnds,columns=["Gene","N","S","dN_dS"]).to_csv(OUTPUT_DIR+"/gene_selection.csv",index=False)

# ============================
# PER-GENE COVERAGE
# ============================
coverage=[]
for g,L in gene_core.items():
    coverage.append([g,L*3])
pd.DataFrame(coverage,columns=["Gene","Core_bp"]).to_csv(OUTPUT_DIR+"/gene_core_coverage.csv",index=False)

# ============================
# PAIRWISE SNP DISTANCES (core)
# ============================
dist=[]
for a,b in combinations(ETHIOPIAN_IDS,2):
    d=0
    for p in core_positions:
        if df.loc[p,a]!=df.loc[p,b]:
            d+=1
    dist.append([a,b,d])

pd.DataFrame(dist,columns=["Sample1","Sample2","Core_SNPs"]).to_csv(OUTPUT_DIR+"/pairwise_core_distances.csv",index=False)

print("\nCORE-GENOME MARV ANALYSIS COMPLETE")
print(f"All results in: {OUTPUT_DIR}")

