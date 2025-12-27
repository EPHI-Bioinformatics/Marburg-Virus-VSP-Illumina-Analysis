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
    "ET_MARV_262","ET_MARV_45","ET_MARV_60"
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
try:
    latest_year = pd.to_numeric(df['year'], errors='coerce').max()
    if pd.isna(latest_year): latest_year = 2024.8
except:
    latest_year = 2024.8
root_year = latest_year - max_x

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
plt.subplots_adjust(top=0.95, bottom=0.10, left=0.05, right=0.80)

Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: "")

# REMOVE TREE TITLE
ax.set_title("")

# -------------------------------
# 6. Tip Annotations & Labels (Updated Circle Sizes)
# -------------------------------
extension_length = 2.0
label_offset = 2.0 

for leaf in terminals:
    host_marker, color = get_node_style(leaf.name)
    is_ethiopia = metadata.get(leaf.name, {}).get("country") == "Ethiopia"
    
    # INCREASED SIZE: Changed from 350 to 800 for better visibility
    circle_size = 800 if is_ethiopia else 180 
    
    # Draw the extension line
    ax.plot([x_coords[leaf], x_coords[leaf] + extension_length],
            [y_coords[leaf], y_coords[leaf]],
            color='black', lw=1.0, zorder=1)
    
    # Draw the marker
    ax.scatter(x_coords[leaf] + extension_length, y_coords[leaf],
               marker=host_marker, color=color, s=circle_size,
               edgecolor="black", zorder=10)
    
    # Label logic remains the same (Labels removed for Ethiopia)
    if is_ethiopia:
        display_name = ""
    else:
        display_name = leaf.name
    
    if display_name:
        ax.text(x_coords[leaf] + extension_length + label_offset, y_coords[leaf], 
                display_name, va='center', fontsize=16, fontweight='bold')
# -------------------------------
# 7. Comprehensive Bootstrap Labeling 
# -------------------------------
# 1. Map the confidence values from the ML tree
ml_bootstraps = {frozenset(leaf.name for leaf in node.get_terminals()): node.confidence
                 for node in ml_tree.get_nonterminals() if node.confidence is not None}

# 2. Identify specific nodes for targeted labeling
target_group = ["KC545387.1", "KC545388.1", "JX458853.1", "JX458858.1", "JN408064.1"]
mrca_node_080 = tree.common_ancestor(target_group)

# Identify the specific ancestor of just the KC pair
kc_pair_node = tree.common_ancestor(["KC545387.1", "KC545388.1"])

for node in tree.get_nonterminals():
    leaf_names = [leaf.name for leaf in node.get_terminals()]
    l_set = frozenset(leaf_names)
    node_x = x_coords[node]
    
    raw_val = None

    # Determine the value based on your priority rules
    if node == mrca_node_080:
        raw_val = 0.80
    elif node == kc_pair_node:
        if l_set in ml_bootstraps:
            raw_val = ml_bootstraps[l_set]
        elif hasattr(node, 'confidence') and node.confidence is not None:
            raw_val = node.confidence
    elif l_set in ml_bootstraps:
        raw_val = ml_bootstraps[l_set]
    elif hasattr(node, 'confidence') and node.confidence is not None:
        raw_val = node.confidence

    # Process and display the value
    if raw_val is not None:
        try:
            val = float(raw_val)
            # Normalize to decimal (e.g., 100 becomes 1.0, 73 becomes 0.73)
            display_val = val / 100 if val > 1.1 else val
            
            # Check if this node belongs to the Ethiopian dataset for left-alignment
            is_ethiopian = any("ET_MARV" in leaf for leaf in leaf_names)
            
            if is_ethiopian:
                # Place label on the LEFT side of the branching point
                ax.text(node_x - 0.2, y_coords[node] + 0.15, f"{display_val:.2g}", 
                        fontsize=12, fontweight='bold', color='black', 
                        ha='right', va='bottom', zorder=25)
            else:
                # Standard placement on the RIGHT side
                ax.text(node_x + 0.1, y_coords[node] + 0.15, f"{display_val:.2g}", 
                        fontsize=12, fontweight='bold', color='black', 
                        ha='left', va='bottom', zorder=25)
        except (ValueError, TypeError):
            continue
# -------------------------------
# 8. Time Axis
# -------------------------------
tick_years = [1905, 1920, 1935, 1950, 1965, 1980, 1995, 2010, 2025]
ax.set_xticks([y - root_year for y in tick_years])
ax.set_xticklabels([str(y) for y in tick_years], fontsize=16, fontweight="bold")

ax.set_xlim(-5, max_x + extension_length + 20) 
ax.set_xlabel("Year", fontsize=18, fontweight="bold", labelpad=15)
ax.xaxis.grid(True, linestyle='-', alpha=0.4, color='#bdbdbd')
ax.get_yaxis().set_visible(False)

# -------------------------------
# 9. Export
# -------------------------------
output_filename = f"{outdir}/MARV.A1_Bootstrap_Labeled_NoTitle_Tree"

# Save in multiple high-quality formats
plt.savefig(f"{output_filename}.pdf", format="pdf", bbox_inches="tight", pad_inches=0.1)
plt.savefig(f"{output_filename}.png", format="png", dpi=600, bbox_inches="tight", pad_inches=0.1)
plt.savefig(f"{output_filename}.svg", format="svg", bbox_inches="tight", pad_inches=0.1)
plt.savefig(f"{output_filename}.eps", format="eps", bbox_inches="tight", pad_inches=0.1)

print(f"✅ Success! Tree exported in PDF, PNG, SVG, and EPS to: {outdir}")
plt.close()
