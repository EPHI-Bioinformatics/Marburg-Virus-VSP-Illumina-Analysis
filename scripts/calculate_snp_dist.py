import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import seaborn as sns  # Added for heatmap
from Bio import AlignIO
from collections import defaultdict

# ---------------------------------------------------------
# 1. CALCULATION FUNCTION (exclude any non-ACGT for clustering)
# ---------------------------------------------------------
def calculate_snp_distance(alignment, min_frac=0.7):
    n = len(alignment)
    ids = [rec.id for rec in alignment]
    genome_len = alignment.get_alignment_length()
    full_dist = np.full((n, n), np.nan)
    callable_sites = np.zeros((n, n), dtype=int)
    
    # Convert alignment to uppercase array for fast masking
    seq_array = np.array([list(str(rec.seq).upper()) for rec in alignment])
    
    # Mask positions where any sample is not A/C/G/T
    valid_mask = np.all(np.isin(seq_array, ["A", "C", "G", "T"]), axis=0)
    
    for i in range(n):
        s1 = seq_array[i]
        callable_sites[i, i] = np.sum(valid_mask)
        for j in range(i + 1, n):
            s2 = seq_array[j]
            # Only consider positions valid in the entire alignment
            pair_mask = valid_mask
            n_callable = np.sum(pair_mask)
            callable_sites[i, j] = callable_sites[j, i] = n_callable
            
            if n_callable > 0 and (n_callable / genome_len) >= min_frac:
                snps = np.sum(s1[pair_mask] != s2[pair_mask])
                full_dist[i, j] = full_dist[j, i] = snps

    return (
        pd.DataFrame(full_dist, index=ids, columns=ids),
        pd.DataFrame(callable_sites, index=ids, columns=ids)
    )

# ---------------------------------------------------------
# 2. EXECUTION & CSV EXPORT
# ---------------------------------------------------------
MSA_FILE = "../results/10_msa/Ethiopian_only/ethiopian_only_msa.fasta"
OUT_DIR = "../results/transmission_dynamics"
os.makedirs(OUT_DIR, exist_ok=True)

MIN_CALLABLE_FRAC = 0.5  # Treat samples <50% callable as low-quality

alignment = AlignIO.read(MSA_FILE, "fasta")
snp_df, callable_df = calculate_snp_distance(alignment, min_frac=0.7)

# Compute callable fraction per sample
sample_callable_fraction = callable_df.sum(axis=1) / callable_df.shape[1]
sample_callable_fraction.to_csv(f"{OUT_DIR}/sample_callable_fraction.csv", header=["Callable_Fraction"])

# Identify low-quality samples (<50% callable)
low_quality_samples = sample_callable_fraction.index[sample_callable_fraction < MIN_CALLABLE_FRAC].tolist()
if low_quality_samples:
    print(f"⚠️ Low-quality samples detected (<{MIN_CALLABLE_FRAC*100}% callable): {low_quality_samples}")

# Export standard CSVs
snp_df.to_csv(f"{OUT_DIR}/full_pairwise_snp_distances.csv")
snp_df.to_csv(f"{OUT_DIR}/pairwise_snp_distances.csv")
callable_df.to_csv(f"{OUT_DIR}/pairwise_callable_sites.csv")

# ---------------------------------------------------------
# 3. BUILD MINIMUM SPANNING TREE
# ---------------------------------------------------------
G_full = nx.Graph()
G_full.add_nodes_from(snp_df.index)

for i in range(len(snp_df)):
    for j in range(i + 1, len(snp_df)):
        d = snp_df.iloc[i, j]
        if not np.isnan(d):
            G_full.add_edge(snp_df.index[i], snp_df.index[j], weight=int(d))

# Minimum spanning tree
G_mst = nx.minimum_spanning_tree(G_full, weight="weight")

# Export MST edges
mst_edges = []
for u, v, d in G_mst.edges(data=True):
    mst_edges.append([u, v, d["weight"]])

pd.DataFrame(
    mst_edges,
    columns=["Sample_1", "Sample_2", "SNP_Distance"]
).to_csv(f"{OUT_DIR}/mst_edges.csv", index=False)

# ---------------------------------------------------------
# 4. ASSIGN GROUP IDs BASED ON MST CONNECTIVITY
# ---------------------------------------------------------
clusters = list(nx.connected_components(G_mst))
cluster_nodes_map = {}
for idx, cluster in enumerate(clusters):
    for node in cluster:
        cluster_nodes_map[node] = idx

# Assign Group IDs
node_group_ids = []
current_singleton_id = len(clusters)
for node in snp_df.index:
    if node in cluster_nodes_map:
        node_group_ids.append(cluster_nodes_map[node])
    else:
        # Low-quality samples become singletons
        node_group_ids.append(current_singleton_id)
        current_singleton_id += 1

# Color mapping: low-quality samples in red
cmap = matplotlib.colormaps['tab20']
plot_node_colors = [
    'red' if node in low_quality_samples else cmap(gid % 20)
    for node, gid in zip(snp_df.index, node_group_ids)
]

# ---------------------------------------------------------
# 5. transmission_clusters.csv
# ---------------------------------------------------------
cluster_export = []
for i, node in enumerate(snp_df.index):
    label = f"Cluster_{node_group_ids[i]}" if node in cluster_nodes_map else "Singleton"
    cluster_export.append([node, label])

pd.DataFrame(cluster_export, columns=["Sample_ID", "Cluster_ID"]).to_csv(
    f"{OUT_DIR}/transmission_clusters.csv", index=False
)

# ---------------------------------------------------------
# 6. cluster_size_summary.csv
# ---------------------------------------------------------
summary = defaultdict(int)
for c in clusters:
    summary[len(c)] += 1

pd.DataFrame(
    list(summary.items()),
    columns=["Cluster_Size", "Number_of_Clusters"]
).to_csv(f"{OUT_DIR}/cluster_size_summary.csv", index=False)

# --- Section 7: VISUALIZATION (Full Network) ---

pos = nx.circular_layout(G_full)  # Circular layout for the full network

plt.figure(figsize=(14, 14))

# Draw all edges
nx.draw_networkx_edges(
    G_full, pos, 
    width=1.5, 
    edge_color="gray", 
    alpha=0.3, 
    style='dashed'
)

# Draw nodes
nx.draw_networkx_nodes(
    G_full, pos,
    node_size=8000,
    node_color=plot_node_colors,
    edgecolors="black",
    linewidths=2
)

# Draw labels
nx.draw_networkx_labels(G_full, pos, font_size=10, font_weight="bold")

# Draw SNP labels for EVERY edge
edge_labels = nx.get_edge_attributes(G_full, "weight")
nx.draw_networkx_edge_labels(
    G_full, pos,
    edge_labels=edge_labels,
    font_size=9,
    label_pos=0.3, 
    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
)

plt.title("Full Pairwise SNP Distances: Marburg Outbreak Connectivity", fontsize=20)
plt.axis("off")
plt.savefig(f"{OUT_DIR}/marburg_complete_network.png", dpi=300, bbox_inches="tight")

# ---------------------------------------------------------
# 7.5 VISUALIZATION: CLUSTERED HEATMAP
# ---------------------------------------------------------
heatmap_df = snp_df.copy()
np.fill_diagonal(heatmap_df.values, 0) # Ensure diagonal is 0 for clustering

# If NaNs exist (e.g. from low-quality samples), fill with max+5 to separate them
if heatmap_df.isna().any().any():
    fill_val = heatmap_df.max().max() + 5
    heatmap_df = heatmap_df.fillna(fill_val)

g = sns.clustermap(
    heatmap_df,
    annot=True, 
    fmt=".0f", 
    cmap="YlOrRd", 
    linewidths=.5,
    figsize=(12, 10),
    cbar_kws={'label': 'SNP Distance'}
)
g.fig.suptitle('Clustered Heatmap of Marburg Virus Pairwise SNP Distances', fontsize=18, y=1.02)
plt.savefig(f"{OUT_DIR}/marburg_snp_heatmap.png", dpi=300, bbox_inches="tight")

# ---------------------------------------------------------
# 8. GENERATE SNP DIFFERENCE TABLE (strict ACGT mask)
# ---------------------------------------------------------
snp_diff_records = []

seq_array = np.array([list(str(rec.seq).upper()) for rec in alignment])
ids = [rec.id for rec in alignment]

valid_mask = np.all(np.isin(seq_array, ["A", "C", "G", "T"]), axis=0)
alignment_len = alignment.get_alignment_length()

for i in range(len(alignment)):
    s1_id = ids[i]
    s1_seq = seq_array[i]
    
    for j in range(i + 1, len(alignment)):
        s2_id = ids[j]
        s2_seq = seq_array[j]
        
        diff_positions = np.where((s1_seq != s2_seq) & valid_mask)[0]
        
        for pos in diff_positions:
            snp_diff_records.append([
                s1_id,
                s2_id,
                pos + 1,        # 1-based genome position
                s1_seq[pos],
                s2_seq[pos]
            ])

# Export SNP differences
snp_diff_df = pd.DataFrame(
    snp_diff_records,
    columns=["Sample_1", "Sample_2", "Position", "Sample_1_Base", "Sample_2_Base"]
)

snp_diff_df.to_csv(f"{OUT_DIR}/pairwise_snp_differences.csv", index=False)
print(f"✅ SNP differences and Heatmap saved to {OUT_DIR}")
