import matplotlib
# Force headless mode BEFORE importing anything else
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
from Bio import Phylo

# 1. PATHS
TREE_FILE = "/media/betselotz/Expansion/MARV-Gen/results/11_phylogeny/MARV.A.1/marburg_ml.treefile"
METADATA_FILE = "/media/betselotz/Expansion/MARV-Gen/metadata/all_seq_metadata.csv"
OUTPUT_FILE = "/media/betselotz/Expansion/MARV-Gen/results/11_phylogeny/MARV.A.1/final_publication_tree.png"

COUNTRY_COLORS = {
    "Ethiopia": "#d62728", "Uganda": "#1f77b4", "DRC": "#2ca02c",
    "Angola": "#9467bd", "Kenya": "#ff7f0e", "South Africa": "#8c564b",
    "Sierra Leone": "#e377c2", "Rwanda": "#bcbd22", "Guinea": "#17becf",
    "Ghana": "#7f7f7f", "Netherlands": "#000000",
    "Germany": "#FFD700", "Canada": "#C0C0C0", "US": "#4B0082"
}

# 2. LOAD METADATA
metadata = {}
with open(METADATA_FILE, mode='r') as f:
    for row in csv.DictReader(f):
        metadata[row['name']] = row

# 3. DRAWING
tree = Phylo.read(TREE_FILE, "newick")

# Set up figure
fig, ax = plt.subplots(figsize=(10, 8))

# Apply colors
for clade in tree.get_terminals():
    data = metadata.get(clade.name, {})
    clade.color = COUNTRY_COLORS.get(data.get("country", ""), "#7f7f7f")

# Draw in headless mode
Phylo.draw(tree, axes=ax, do_show=False)

plt.title("Marburg Virus Phylogeny (SNP-based)")
plt.savefig(OUTPUT_FILE, dpi=300)
print(f">>> Tree saved to {OUTPUT_FILE}")
