import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo
import pandas as pd

# -------------------------------
# Matplotlib font settings
# -------------------------------
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.use('Agg')

# -------------------------------
# 1. Load Data & Directories
# -------------------------------
tree_path = "../results/12_treetime/MARV.A.1/timetree.nexus"
ml_tree_path = "../results/11_phylogeny/MARV.A.1/marburg_ml.treefile"
metadata_path = "../metadata/all_seq_metadata.csv"
outdir = "../results/12_treetime/MARV.A.1/visualization"
os.makedirs(outdir, exist_ok=True)

tree = Phylo.read(tree_path, "nexus")
ml_tree = Phylo.read(ml_tree_path, "newick")

try:
    df = pd.read_csv(metadata_path)
    df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    metadata = df.set_index("name").to_dict("index")
except FileNotFoundError:
    metadata = {}

# -------------------------------
# 2. Define MARV.A.1 Clade & Prune Tree
# -------------------------------
MARV_A1_SAMPLES = [
    "JN408064.1","JX458851.1","KC545387.1","KC545388.1","JX458853.1",
    "JX458858.1","ET_MARV_23","ET_MARV_31","ET_MARV_32","ET_MARV_261",
    "ET_MARV_262","ET_MARV_45","ET_MARV_60", "EF446132.1"
]

all_tips = [t.name for t in tree.get_terminals()]
for tip in all_tips:
    if tip not in MARV_A1_SAMPLES:
        tree.prune(tip)

# -------------------------------
# 3. Coordinate & Time Calculations
# -------------------------------
x_coords = tree.depths()
terminals = tree.get_terminals()
y_coords = {tip: i for i, tip in enumerate(terminals, 1)}
for node in tree.get_nonterminals(order='postorder'):
    y_coords[node] = sum(y_coords[child] for child in node.clades) / len(node.clades)

max_x = max(x_coords.values())

# -------------------------------
# 4. Styling Maps
# -------------------------------
COUNTRY_COLORS = {
    "Ethiopia": "#d62728", "Uganda": "#1f77b4", "DRC": "#2ca02c",
    "Angola": "#9467bd", "Kenya": "#ff7f0e", "South Africa": "#8c564b",
    "Sierra Leone": "#e377c2", "Rwanda": "#bcbd22", "Guinea": "#17becf",
    "Ghana": "#7f7f7f", "Netherlands": "#000000",
    "Germany": "#FFD700", "Canada": "#C0C0C0", "US": "#4B0082"
}

def get_node_style(strain_name):
    row = metadata.get(strain_name, {})
    source = str(row.get("source", "")).lower()
    country = str(row.get("country", ""))
    host_marker = "o" if "human" in source else ("D" if "bat" in source else ">")
    color = COUNTRY_COLORS.get(country, "#7f7f7f")
    return host_marker, color

# -------------------------------
# 5. Figure Setup
# -------------------------------
fig = plt.figure(figsize=(26, 14)) 
ax = fig.add_subplot(1, 1, 1)
# Adjusted bottom margin since x-axis labels are gone
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.05, right=0.80)

Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: "")
ax.set_title("")

# -------------------------------
# 6. Tip Annotations & Labels
# -------------------------------
extension_length = 2.0
label_offset = 2.0 

for leaf in terminals:
    host_marker, color = get_node_style(leaf.name)
    is_ethiopia = metadata.get(leaf.name, {}).get("country") == "Ethiopia"
    circle_size = 800 if is_ethiopia else 180 
    
    ax.plot([x_coords[leaf], x_coords[leaf] + extension_length],
            [y_coords[leaf], y_coords[leaf]],
            color='black', lw=1.0, zorder=1)
    
    ax.scatter(x_coords[leaf] + extension_length, y_coords[leaf],
               marker=host_marker, color=color, s=circle_size,
               edgecolor="black", zorder=10)
    
    display_name = "" if is_ethiopia else leaf.name
    if display_name:
        ax.text(x_coords[leaf] + extension_length + label_offset, y_coords[leaf], 
                display_name, va='center', fontsize=16, fontweight='bold')


# -------------------------------
# 8. Axis Formatting (X-Axis Removed)
# -------------------------------
# Maintain the same scale/view but hide all x-axis elements
ax.set_xlim(-5, max_x + extension_length + 20) 

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# Remove the outer frame (spines) for a clean look
for spine in ax.spines.values():
    spine.set_visible(False)

# -------------------------------
# 9. Export
# -------------------------------
output_filename = f"{outdir}/MARV.A1_Clean_Phylogeny_NoAxis"

for fmt in ["pdf", "png", "svg", "eps"]:
    plt.savefig(f"{output_filename}.{fmt}", format=fmt, 
                dpi=600 if fmt=="png" else None, 
                bbox_inches="tight", pad_inches=0.1)

print(f"✅ Success! Tree exported without x-axis to: {outdir}")
plt.close()
