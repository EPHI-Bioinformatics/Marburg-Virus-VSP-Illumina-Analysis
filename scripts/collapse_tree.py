import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from Bio import Phylo
import pandas as pd
import numpy as np

# -------------------------------
# Matplotlib settings
# -------------------------------
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.use('Agg')

# -------------------------------
# 1. Load Data & Directories
# -------------------------------
tree_path = "../results/12_treetime/timetree.nexus"
ml_tree_path = "../results/11_phylogeny/marburg_ml.treefile"
metadata_path = "../metadata/all_seq_metadata.csv"
outdir = "../results/12_treetime/visualization_eth_nld"
os.makedirs(outdir, exist_ok=True)

tree = Phylo.read(tree_path, "nexus")
ml_tree = Phylo.read(ml_tree_path, "newick")

df = pd.read_csv(metadata_path)
df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
metadata = df.set_index("name").to_dict("index")

# -------------------------------
# 2. Filter, Prune, and Root
# -------------------------------
# Define targets
TARGET_COUNTRIES = ["Ethiopia", "Netherlands"]
target_taxa = [name for name, info in metadata.items() if info.get("country") in TARGET_COUNTRIES]

# Prune tips that are NOT in our target list
all_tips = [t.name for t in tree.get_terminals()]
for tip in all_tips:
    if tip not in target_taxa:
        tree.prune(tip)

# Dynamic Outgroup Selection: Find the latest sample among the targets
target_df = df[df['country'].isin(TARGET_COUNTRIES)].copy()
target_df['year'] = pd.to_numeric(target_df['year'], errors='coerce')
latest_sample_name = target_df.loc[target_df['year'].idxmax(), 'name']

# Re-root the tree
tree.root_with_outgroup({"name": latest_sample_name})

# -------------------------------
# 3. Coordinate Calculation
# -------------------------------
def calculate_coords(tr):
    x = tr.depths()
    terms = tr.get_terminals()
    y = {tip: i for i, tip in enumerate(terms, 1)}
    for node in tr.get_nonterminals(order='postorder'):
        y[node] = sum(y[child] for child in node.clades) / len(node.clades)
    return x, y, terms

x_coords, y_coords, terminals = calculate_coords(tree)

# -------------------------------
# 4. Styling & Figure Setup
# -------------------------------
COUNTRY_COLORS = {"Ethiopia": "#d62728", "Netherlands": "#000000"}
fig = plt.figure(figsize=(16, 16))
ax = fig.add_subplot(1, 1, 1)

# Draw the tree
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: "")

# Annotate tips
for leaf in terminals:
    country = metadata.get(leaf.name, {}).get("country", "Unknown")
    color = COUNTRY_COLORS.get(country, "#7f7f7f")
    ax.scatter(x_coords[leaf], y_coords[leaf], color=color, s=150, zorder=10, edgecolors="black")

# -------------------------------
# 5. Formatting
# -------------------------------
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("Time (Years)", fontsize=14, fontweight="bold")
plt.title(f"Phylogeny Rooted at {latest_sample_name}", fontsize=16)

# Save
plt.savefig(f"{outdir}/Ethiopia_Netherlands_Tree.pdf", bbox_inches="tight")
plt.savefig(f"{outdir}/Ethiopia_Netherlands_Tree.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"Success: Tree pruned to {len(target_taxa)} sequences and rooted at {latest_sample_name}.")
