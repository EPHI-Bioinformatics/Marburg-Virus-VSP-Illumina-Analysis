import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from Bio import AlignIO
from collections import defaultdict

# ---------------------------------------------------------
# 1. CALCULATION FUNCTION (exclude N for clustering)
# ---------------------------------------------------------
def calculate_snp_distance(alignment, min_frac=0.7):
    n = len(alignment)
    ids = [rec.id for rec in alignment]
    genome_len = alignment.get_alignment_length()
    full_dist = np.full((n, n), np.nan)
    callable_sites = np.zeros((n, n), dtype=int)
    
    for i in range(n):
        s1 = np.array(list(str(alignment[i].seq).upper()))
        callable_sites[i, i] = np.sum((s1 != "-") & (s1 != "N"))
        for j in range(i + 1, n):
            s2 = np.array(list(str(alignment[j].seq).upper()))
            mask = ((s1 != "-") & (s2 != "-") & (s1 != "N") & (s2 != "N"))
            n_callable = np.sum(mask)
            callable_sites[i, j] = callable_sites[j, i] = n_callable
            
            if n_callable > 0 and (n_callable / genome_len) >= min_frac:
                snps = np.sum(s1[mask] != s2[mask])
                full_dist[i, j] = full_dist[j, i] = snps
            else:
                full_dist[i, j] = full_dist[j, i] = np.nan

    return (pd.DataFrame(full_dist, index=ids, columns=ids), 
            pd.DataFrame(callable_sites, index=ids, columns=ids))

# ---------------------------------------------------------
# 2. EXECUTION & CSV EXPORT
# ---------------------------------------------------------
MSA_FILE = "../results/10_msa/Ethiopian_only/ethiopian_only_msa.fasta"
OUT_DIR = "../results/transmission_dynamics"
os.makedirs(OUT_DIR, exist_ok=True)

SNP_THRESHOLD = 2
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
# 3. CLUSTER DETECTION
# ---------------------------------------------------------
G_cluster = nx.Graph()
G_cluster.add_nodes_from(snp_df.index)

for i in range(len(snp_df)):
    for j in range(i + 1, len(snp_df)):
        if not np.isnan(snp_df.iloc[i,j]) and snp_df.iloc[i,j] <= SNP_THRESHOLD:
            G_cluster.add_edge(snp_df.index[i], snp_df.index[j])

clusters = list(nx.connected_components(G_cluster))
cluster_nodes_map = {}
for idx, cluster in enumerate(clusters):
    for node in cluster:
        cluster_nodes_map[node] = idx

# Assign Group IDs to every node
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
# 4. transmission_clusters.csv
# ---------------------------------------------------------
cluster_export = []
for i, node in enumerate(snp_df.index):
    label = f"Cluster_{node_group_ids[i]}" if node in cluster_nodes_map else "Singleton"
    cluster_export.append([node, label])
pd.DataFrame(cluster_export, columns=["Sample_ID", "Cluster_ID"]).to_csv(
    f"{OUT_DIR}/transmission_clusters.csv", index=False
)

# ---------------------------------------------------------
# 5. cluster_size_summary.csv
# ---------------------------------------------------------
summary = defaultdict(int)
for c in clusters:
    summary[len(c)] += 1
pd.DataFrame(
    list(summary.items()),
    columns=["Cluster_Size", "Number_of_Clusters"]
).to_csv(f"{OUT_DIR}/cluster_size_summary.csv", index=False)

# ---------------------------------------------------------
# 6. VISUALIZATION
# ---------------------------------------------------------
G = nx.Graph()
G.add_nodes_from(snp_df.index)
for i in range(len(snp_df)):
    for j in range(i + 1, len(snp_df)):
        d = snp_df.iloc[i, j]
        if not np.isnan(d):
            G.add_edge(
                snp_df.index[i],
                snp_df.index[j],
                snp=int(d),
                length=1 + d
            )

pos = nx.kamada_kawai_layout(G, weight='length')

# Shift cluster 0 nodes slightly to the right
cluster_0_nodes = [node for node, gid in zip(snp_df.index, node_group_ids) if gid == 0]
for node in pos:
    if node in cluster_0_nodes:
        pos[node][0] += 0.5

plt.figure(figsize=(18, 18))

# Draw edges
nx.draw_networkx_edges(G, pos, width=1.2, edge_color="gray", alpha=0.2)
triangle_edges = [(u, v) for u, v in G.edges() if u in cluster_0_nodes and v in cluster_0_nodes]
nx.draw_networkx_edges(G, pos, edgelist=triangle_edges, width=4.5, edge_color="black")

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_size=9000, node_color=plot_node_colors, edgecolors="black", linewidths=2.5)
nx.draw_networkx_labels(G, pos, font_size=10, font_weight="bold")

# Draw SNP edge labels
edge_labels = nx.get_edge_attributes(G, 'snp')
nx.draw_networkx_edge_labels(
    G, pos, edge_labels=edge_labels,
    font_color='black', font_size=15, font_weight='bold',
    label_pos=0.5, rotate=False,
    bbox=dict(facecolor='white', alpha=0.9, edgecolor='none', boxstyle='round,pad=0.3')
)

plt.title("Marburg Virus Transmission Cluster Analysis", fontsize=26, pad=30)
plt.axis("off")
plt.savefig(f"{OUT_DIR}/marburg_triangular_clusters.png", dpi=300, bbox_inches="tight")

print(f"✅ Success! Generated plot and all CSVs (including low-quality singletons) in {OUT_DIR}")

