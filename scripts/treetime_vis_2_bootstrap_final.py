import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from Bio import Phylo
import pandas as pd

# -------------------------------
# Matplotlib font settings for publication
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
    "MARV.A.1": ["ET_MARV_262", "ET_MARV_261", "ET_MARV_23", "ET_MARV_32", "ET_MARV_31", "ET_MARV_45", "ET_MARV_60", "JN408064.1", "JX458851.1", "JX458853.1", "JX458858.1", "KC545387.1", "KC545388.1", "EF446132.1"],
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
# 2.5 Targeted Node Rotation
# -------------------------------
mrca_kc = tree.common_ancestor("KC545387.1", "KC545388.1")
mrca_jx = tree.common_ancestor("JX458853.1", "JX458858.1")
parent = tree.common_ancestor(mrca_kc, mrca_jx)

if mrca_kc in parent.clades and mrca_jx in parent.clades:
    other_siblings = [c for c in parent.clades if c not in [mrca_kc, mrca_jx]]
    parent.clades = [mrca_kc, mrca_jx] + other_siblings

mrca_kc.clades.sort(key=lambda x: x.name if x.name else "")
mrca_jx.clades.sort(key=lambda x: x.name if x.name else "")

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
# 3.5 Time Calculation
# -------------------------------
max_x = max(x_coords.values())
try:
    latest_year = pd.to_numeric(df['year'], errors='coerce').max()
    if pd.isna(latest_year): latest_year = 2025.11
except:
    latest_year = 2025.11
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
    host_marker = "o" if "human" in source else ("D" if "bat" in source else ">")
    color = COUNTRY_COLORS.get(country, "#7f7f7f")
    return host_marker, color

# -------------------------------
# 5. Figure Setup (Matplotlib)
# -------------------------------
fig = plt.figure(figsize=(26, 26))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(top=0.98, bottom=0.05, left=0.03, right=0.94)

Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: "")
ax.set_title("") 

# -------------------------------
# 6. Tip Annotations
# -------------------------------
extension_length = 2.0
for leaf in terminals:
    host_marker, color = get_node_style(leaf.name)
    is_ethiopia = metadata.get(leaf.name, {}).get("country") == "Ethiopia"
    circle_size = 450 if is_ethiopia else 250 
    
    ax.plot([x_coords[leaf], x_coords[leaf] + extension_length], [y_coords[leaf], y_coords[leaf]],
            color='black', lw=1.2, zorder=1)
    ax.scatter(x_coords[leaf] + extension_length, y_coords[leaf], marker=host_marker, color=color,
               s=circle_size, edgecolor="black", zorder=10)

# -------------------------------
# 7. Clade Bars
# -------------------------------
bar_x_start = max_x + extension_length + 2.5
for clade in ["MARV.A.1", "MARV.A.3", "MARV.B", "MARV.B.2", "RAVV.1", "RAVV.2"]:
    t_names = SPECIFIC_CLADES.get(clade, [])
    c_y = [y_coords[t] for t in terminals if t.name in t_names]
    if c_y:
        y_min, y_max = min(c_y), max(c_y)
        col = CLADE_STYLES.get(clade, CLADE_STYLES["default"])["color"] 
        ax.add_patch(Rectangle((bar_x_start, y_min - 0.4), 2.5, (y_max - y_min) + 0.8, 
                               color=col, alpha=0.8, zorder=5))
        ax.text(bar_x_start + 3.5, (y_min + y_max) / 2, clade, 
                va='center', fontsize=20, fontweight='bold', color=col)

# -------------------------------
# 8. Bootstraps
# -------------------------------
ml_bootstraps = {frozenset(leaf.name for leaf in node.get_terminals()): node.confidence 
                 for node in ml_tree.get_nonterminals() if node.confidence is not None}

ancestral_cutoff = max_x * 0.75 
for node in tree.get_nonterminals():
    leaf_names = [leaf.name for leaf in node.get_terminals()]
    l_set = frozenset(leaf_names)
    if l_set in ml_bootstraps:
        val = ml_bootstraps[l_set]
        scaled_conf = val/100 if val > 1.1 else val
        if (x_coords[node] < ancestral_cutoff and len(leaf_names) >= 3) or len(leaf_names) > 15:
            label_str = f"{scaled_conf:g}" if scaled_conf >= 1.0 else f"{scaled_conf:.2f}"
            ax.text(x_coords[node] + 0.2, y_coords[node] + 0.1, label_str, 
                    fontsize=14, color="#222222", weight="bold", ha="left", va="bottom", zorder=25)

# -------------------------------
# 9. Styling (REVISED TO REMOVE Y-AXIS LINE)
# -------------------------------
tick_years = [1905, 1920, 1935, 1950, 1965, 1980, 1995, 2010, 2025]
ax.set_xticks([y - root_year for y in tick_years])
ax.set_xticklabels([str(y) for y in tick_years], fontsize=18, fontweight="bold")

# Hide y-axis elements
ax.get_yaxis().set_visible(False)
ax.spines['left'].set_visible(False) # <--- REMOVES THE VERTICAL TAXA LINE

ax.set_xlabel("Year", fontsize=24, fontweight="bold", labelpad=15)
ax.spines['bottom'].set_linewidth(2.0)
ax.set_xlim(-5, bar_x_start + 10)
ax.xaxis.grid(True, linestyle='--', alpha=0.4, color='#bdbdbd')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# -------------------------------
# 10. Legends
# -------------------------------
leg_x, leg_y = 0.02, 0.02
s_el = [Line2D([0], [0], color="none", label="Host Source (tips shapes):"),
        Line2D([0], [0], marker="o", color="none", label="Human", mfc="none", mec="black", ms=14),
        Line2D([0], [0], marker="D", color="none", label="Bat", mfc="none", mec="black", ms=12)]
l1 = ax.legend(handles=s_el, loc="lower left", bbox_to_anchor=(leg_x, leg_y), bbox_transform=ax.transAxes, fontsize=20, frameon=False); ax.add_artist(l1)

filtered_countries = [c for c in COUNTRY_COLORS.keys() if c not in {"Canada", "US"}]
c_el = [Line2D([0], [0], color="none", label="Geographic Origin (tips colour):")] + \
        [Line2D([0], [0], marker="s", color="none", label=c, mfc=COUNTRY_COLORS[c], ms=12) for c in filtered_countries]
l2 = ax.legend(handles=c_el, loc="lower left", bbox_to_anchor=(leg_x + 0.20, leg_y), bbox_transform=ax.transAxes, fontsize=20, frameon=False, ncol=1); ax.add_artist(l2)

# -------------------------------
# 11. Final High-Quality Static Export
# -------------------------------
output_base = f"{outdir}/MARV_Final_Ordered_FullTree"
plt.savefig(f"{output_base}.png", format="png", dpi=600, bbox_inches="tight", pad_inches=0.1)
plt.savefig(f"{output_base}.pdf", format="pdf", bbox_inches="tight", pad_inches=0.1)
plt.savefig(f"{output_base}.svg", format="svg", bbox_inches="tight", pad_inches=0.1)
plt.savefig(f"{output_base}.eps", format="eps", bbox_inches="tight", pad_inches=0.1)
plt.close()

# -------------------------------
# 12. Interactive HTML Export (Plotly)
# -------------------------------
import plotly.graph_objects as go

def generate_interactive_tree():
    fig_html = go.Figure()

    # 1. Draw Branches
    for node in tree.get_nonterminals():
        for child in node.clades:
            fig_html.add_trace(go.Scatter(x=[x_coords[node], x_coords[child]], y=[y_coords[child], y_coords[child]], mode='lines', line=dict(color='black', width=1.5), hoverinfo='none', showlegend=False))
            fig_html.add_trace(go.Scatter(x=[x_coords[node], x_coords[node]], y=[y_coords[node], y_coords[child]], mode='lines', line=dict(color='black', width=1.5), hoverinfo='none', showlegend=False))

    # 2. Draw Tips
    for leaf in terminals:
        host_marker, color = get_node_style(leaf.name)
        plotly_marker = "circle" if host_marker == "o" else "diamond" if host_marker == "D" else "triangle-right"
        is_ethiopia = metadata.get(leaf.name, {}).get("country") == "Ethiopia"
        fig_html.add_trace(go.Scatter(x=[x_coords[leaf] + extension_length], y=[y_coords[leaf]], mode='markers', marker=dict(symbol=plotly_marker, size=15 if is_ethiopia else 10, color=color, line=dict(width=1, color='black')), name=leaf.name, hovertemplate=f"<b>Accession:</b> {leaf.name}<br><b>Country:</b> {metadata.get(leaf.name, {}).get('country', 'Unknown')}<br><b>Source:</b> {metadata.get(leaf.name, {}).get('source', 'Unknown')}<extra></extra>"))

    # 3. Draw Bootstraps
    for node in tree.get_nonterminals():
        leaf_names = [leaf.name for leaf in node.get_terminals()]
        l_set = frozenset(leaf_names)
        if l_set in ml_bootstraps:
            val = ml_bootstraps[l_set]
            scaled_conf = val/100 if val > 1.1 else val
            fig_html.add_trace(go.Scatter(x=[x_coords[node]], y=[y_coords[node]], mode='markers', marker=dict(size=8, color='rgba(0,0,0,0)'), hovertemplate=f"Bootstrap: {scaled_conf:.2f}<extra></extra>", showlegend=False))

    # 4. Clade Rectangles
    for clade, t_names in SPECIFIC_CLADES.items():
        c_y = [y_coords[t] for t in terminals if t.name in t_names]
        if c_y:
            y_min, y_max = min(c_y), max(c_y)
            col = CLADE_STYLES.get(clade, CLADE_STYLES["default"])["color"]
            fig_html.add_shape(type="rect", x0=bar_x_start, x1=bar_x_start + 2.5, y0=y_min - 0.4, y1=y_max + 0.4, fillcolor=col, opacity=0.8, line_width=0, layer="below")
            fig_html.add_annotation(x=bar_x_start + 3.5, y=(y_min + y_max) / 2, text=f"<b>{clade}</b>", showarrow=False, xanchor="left", font=dict(color=col, size=16))

    # 5. Layout Styling (Removed Y-Axis Line)
    fig_html.update_layout(
        template="simple_white", width=1200, height=1200,
        xaxis=dict(title="Year", tickmode='array', tickvals=[y - root_year for y in tick_years], ticktext=[str(y) for y in tick_years], range=[-5, bar_x_start + 15], showgrid=True, gridcolor='lightgrey', gridwidth=1, mirror=True, linecolor='black', linewidth=2),
        yaxis=dict(title="", showticklabels=False, showline=False, autorange="reversed"), # <--- REMOVED LINE
        showlegend=False, margin=dict(l=50, r=50, t=50, b=50)
    )
    fig_html.write_html(f"{output_base}.html")

generate_interactive_tree()
