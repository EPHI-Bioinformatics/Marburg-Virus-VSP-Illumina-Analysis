import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from Bio import Phylo
import pandas as pd

# -------------------------------
# Matplotlib font settings
# -------------------------------
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
# Using Agg backend for high-quality file generation without a GUI
matplotlib.use('Agg')

# -------------------------------
# 1. Load Data & Directories
# -------------------------------
tree_path = "../results/12_treetime/timetree.nexus"
ml_tree_path = "../results/11_phylogeny/marburg_ml.treefile"
metadata_path = "../metadata/all_seq_metadata.csv"
outdir = "../results/12_treetime/visualization"
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
# 2. Clade Definitions & Pruning
# -------------------------------
SPECIFIC_CLADES = {
    "MARV.A.1": ["ET_MARV_262", "ET_MARV_261", "ET_MARV_23", "ET_MARV_32", "ET_MARV_31", "ET_MARV_45", "ET_MARV_60", "JN408064.1", "JX458851.1", "JX458853.1", "JX458858.1", "KC545387.1", "KC545388.1"],
    "MARV.B": ["JX458825.1", "JX458826.1", "JX458829.1", "JX458828.1", "JX458827.1", "JX458833.1", "JX458832.1", "JX458831.1", "DQ447650.1", "JX458830.1", "JX458849.1", "JX458834.1", "JX458835.1", "JX458847.1", "JX458850.1", "DQ447651.1", "JX458845.1", "JX458846.1", "JX458844.1", "JX458836.1", "JX458837.1", "JX458841.1", "JX458843.1", "JX458848.1", "JX458842.1", "JX458839.1", "JX458840.1", "JX458838.1"],
    "MARV.A.3": ["MN258361.1", "MN187404.1", "MN187406.1", "MN187403.1", "MN187405.1", "OK665848.1", "OL702894.1", "OQ672470.1", "OQ672471.1", "DQ447659.1", "DQ447660.1", "KY047763.1", "DQ447654.1", "DQ447656.1", "DQ447657.1", "DQ447653.1", "DQ447655.1", "DQ447658.1", "KY425629.1", "KU978782.1", "MT586762.1", "KR867677.1", "KR063674.1"],
    "MARV.B.2": ["JX458856.1", "JX458852.1", "FJ750959.1", "FJ750958.1", "JX45884.1", "MH638314.1", "MH638315.1", "FJ750957.1", "JX458855.1", "KP985768.1", "PQ552727", "PQ552726", "PQ552735", "PQ552738", "PQ552731", "PQ552736", "PQ552737", "PQ552730", "PQ552725", "PQ552742", "PQ552739", "PQ552733", "PQ552728", "PQ552732", "PQ552729", "PQ552741", "PQ552734", "PQ552740"],
    "RAVV.2": ["MT321489.1", "DQ447652.1", "FJ750956.1", "JX458857.1", "FJ750954.1", "FJ750953.1", "FJ750955.1"],
    "RAVV.1": ["KU179482.1", "DQ447649.1", "NC_024781.1"]
}
all_clade_taxa = [taxa for sublist in SPECIFIC_CLADES.values() for taxa in sublist]
all_tips = [t.name for t in tree.get_terminals()]
taxa_to_keep = [t for t in all_tips if t in all_clade_taxa]
for tip in all_tips:
    if tip not in taxa_to_keep:
        tree.prune(tip)

# -------------------------------
# 2.5 Manual Tip Shuffling
# -------------------------------
# Find the common ancestor of the two tips you want to swap
mrca = tree.common_ancestor("JN408064.1", "JX458851.1")

# Rotate by reversing the order of the child clades
mrca.clades.reverse()

# -------------------------------
# 3. Coordinate & Time Calculation
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

CLADE_STYLES = {
    "MARV.A.1": {"color": "#3C5488"}, 
    "MARV.A.3": {"color": "#00A087"},
    "MARV.B": {"color": "#E64B35"}, 
    "MARV.B.2": {"color": "#F39B7F"},
    "RAVV.1": {"color": "#1f77b4"}, 
    "RAVV.2": {"color": "#654321"},
    "default": {"color": "#000000"}
}

def get_node_style(strain_name):
    row = metadata.get(strain_name, {})
    source = str(row.get("source", "")).lower()
    country = str(row.get("country", ""))
    if "human" in source: host_marker = "o"
    elif "bat" in source: host_marker = "D"
    elif "other" in source: host_marker = ">"
    else: host_marker = "o"
    color = COUNTRY_COLORS.get(country, "#7f7f7f")
    return host_marker, color

# -------------------------------
# 5. Figure Setup
# -------------------------------
fig = plt.figure(figsize=(26, 26))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(top=0.98, bottom=0.05, left=0.03, right=0.94)
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: "")

# -------------------------------
# 6. Tip Annotations
# -------------------------------
extension_length = 2.0
for leaf in terminals:
    host_marker, color = get_node_style(leaf.name)
    is_ethiopia = metadata.get(leaf.name, {}).get("country") == "Ethiopia"
    circle_size = 450 if is_ethiopia else 250 
    
    ax.plot([x_coords[leaf], x_coords[leaf] + extension_length], [y_coords[leaf], y_coords[leaf]],
            color='black', lw=1.0, zorder=1)
    ax.scatter(x_coords[leaf] + extension_length, y_coords[leaf], marker=host_marker, color=color,
               s=circle_size, edgecolor="black", zorder=10)

# 7. Clade Bars 
# -------------------------------
bar_x_start = max_x + extension_length + 2.5
for clade in ["MARV.A.1", "MARV.A.3", "MARV.B", "MARV.B.2", "RAVV.1", "RAVV.2"]:
    t_names = SPECIFIC_CLADES.get(clade, [])
    c_y = [y_coords[t] for t in terminals if t.name in t_names]
    if c_y:
        y_min, y_max = min(c_y), max(c_y)
        # This line fetches the color defined in CLADE_STYLES above
        col = CLADE_STYLES.get(clade, CLADE_STYLES["default"])["color"] 
        
        ax.add_patch(Rectangle((bar_x_start, y_min - 0.4), 2.5, (y_max - y_min) + 0.8, 
                               color=col, alpha=0.8, zorder=5))
        ax.text(bar_x_start + 3.5, (y_min + y_max) / 2, clade, 
                va='center', fontsize=18, fontweight='bold', color=col)

# -------------------------------
# 8. Bootstraps & Text Labels (Modified)
# -------------------------------
ml_bootstraps = {}
for node in ml_tree.get_nonterminals():
    if node.confidence is not None:
        ml_bootstraps[frozenset(leaf.name for leaf in node.get_terminals())] = node.confidence

ancestral_cutoff = max_x * 0.75 

for node in tree.get_nonterminals():
    leaf_names = [leaf.name for leaf in node.get_terminals()]
    l_set = frozenset(leaf_names)
    num_leaves = len(leaf_names)
    node_x = x_coords[node]
    
    if l_set in ml_bootstraps:
        val = ml_bootstraps[l_set]
        scaled_conf = val/100 if val > 1.1 else val
        
        # --- SHAPES REMOVED ---
        # The ax.scatter call that drew the Diamonds/Squares/Circles has been deleted.

        show_label = False
        if node_x < ancestral_cutoff and num_leaves >= 3:
            show_label = True
        if num_leaves > 15:
            show_label = True
            
        if show_label:
            label_str = f"{scaled_conf:g}" if scaled_conf >= 1.0 else f"{scaled_conf:.2f}"
            # Placement adjusted slightly to center on the node since the shape is gone
            ax.text(node_x + 0.2, y_coords[node] + 0.1, label_str, 
                    fontsize=14, color="#222222", weight="bold",
                    ha="left", va="bottom", zorder=25)
# -------------------------------
# 9. Time Axis & Y-Axis Styling
# -------------------------------
# X-Axis Styling (Year)
tick_years = [1905, 1920, 1935, 1950, 1965, 1980, 1995, 2010, 2025]
ax.set_xticks([y - root_year for y in tick_years])
ax.set_xticklabels([str(y) for y in tick_years], fontsize=18, fontweight="bold")
ax.set_xlabel("Year", fontsize=22, fontweight="bold", labelpad=15)

# Y-Axis Styling (Taxa/Labels)
# Setting the Y-label and formatting tick labels (20, 40, 60, etc.)
ax.set_ylabel("Taxa Count", fontsize=22, fontweight="bold", labelpad=15)
for label in ax.get_yticklabels():
    label.set_fontsize(18)
    label.set_fontweight("bold")

# Make the actual axis lines (spines) bolder
ax.spines['left'].set_linewidth(2.0)
ax.spines['bottom'].set_linewidth(2.0)

left_margin_shift = -5 
ax.set_xlim(left_margin_shift, bar_x_start + 8)
ax.xaxis.grid(True, linestyle='-', alpha=0.4, color='#bdbdbd')
ax.set_title("")

# -------------------------------
# 10. Legends (Modified)
# -------------------------------
leg_x, leg_y = 0.02, 0.02
s_el = [Line2D([0], [0], color="none", label="Host Source (tips shapes):"),
        Line2D([0], [0], marker="o", color="none", label="Human", mfc="none", mec="black", ms=14),
        Line2D([0], [0], marker="D", color="none", label="Bat", mfc="none", mec="black", ms=12)]
l1 = ax.legend(handles=s_el, loc="lower left", bbox_to_anchor=(leg_x, leg_y), bbox_transform=ax.transAxes, fontsize=20, frameon=False); ax.add_artist(l1)

excluded_countries = {"Germany", "Canada", "US"}
filtered_countries = [c for c in COUNTRY_COLORS.keys() if c not in excluded_countries]

c_el = [Line2D([0], [0], color="none", label="Geographic Origin (tips colour):")] + \
        [Line2D([0], [0], marker="s", color="none", label=c, mfc=COUNTRY_COLORS[c], ms=12) for c in filtered_countries]

l2 = ax.legend(handles=c_el, loc="lower left", bbox_to_anchor=(leg_x + 0.20, leg_y), bbox_transform=ax.transAxes, fontsize=20, frameon=False, ncol=1); ax.add_artist(l2)

# --- BOOTSTRAP LEGEND REMOVED ---
# The b_el list and the third ax.legend call have been removed.

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# -------------------------------
# 11. Final High-Quality Export
# -------------------------------
output_filename = f"{outdir}/MARV_Final_FullFill_Tree_Big_Pruned"

# Save as PDF (Vector format - best for publication)
plt.savefig(f"{output_filename}.pdf", format="pdf", bbox_inches="tight", pad_inches=0.1)

# Save as PNG (600 DPI - ultra high quality for presentations/web)
plt.savefig(f"{output_filename}.png", format="png", dpi=600, bbox_inches="tight", pad_inches=0.1)

# Save as SVG (Vector format - best for further editing in Illustrator/Inkscape)
plt.savefig(f"{output_filename}.svg", format="svg", bbox_inches="tight", pad_inches=0.1)

# Save as EPS (Vector format - often required for certain journal submissions)
plt.savefig(f"{output_filename}.eps", format="eps", bbox_inches="tight", pad_inches=0.1)

print(f"Tree exported successfully in PDF, PNG (600 DPI), SVG, and EPS formats to: {outdir}")
plt.close()
