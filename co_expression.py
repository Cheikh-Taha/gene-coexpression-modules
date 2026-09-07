import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import pearsonr
import community as community_louvain

# =====================================================
# 1. LOAD DATA
# =====================================================

print("Loading data...")

df = pd.read_csv(
    "GSE60424_norm_counts_TPM_GRCh38.p13_NCBI.tsv",
    sep="\t"
)

# GeneID as index
df = df.set_index("GeneID")

print("Original shape:", df.shape)

# =====================================================
# 2. REMOVE LOW EXPRESSED GENES
# =====================================================

df = df[df.mean(axis=1) > 1]

print("After filtering:", df.shape)

# =====================================================
# 3. KEEP MOST VARIABLE GENES
# =====================================================

N_GENES = 1000

variance = df.var(axis=1)

top_genes = variance.nlargest(N_GENES).index

expr = df.loc[top_genes]

print("Selected genes:", expr.shape)

# =====================================================
# 4. CORRELATION MATRIX
# =====================================================

print("Computing correlation matrix...")

corr_matrix = expr.T.corr(method="pearson")

print("Correlation matrix shape:", corr_matrix.shape)

# =====================================================
# 5. BUILD GRAPH
# =====================================================

print("Building graph...")

G = nx.Graph()

threshold = 0.85

genes = corr_matrix.index.tolist()

for i in range(len(genes)):
    for j in range(i + 1, len(genes)):

        corr = corr_matrix.iloc[i, j]

        if corr >= threshold:
          G.add_edge(
            str(genes[i]),
            str(genes[j]), 
            weight=float(corr)
    )

print("\nGRAPH STATISTICS")
print("----------------")
print("Nodes:", G.number_of_nodes())
print("Edges:", G.number_of_edges())

# =====================================================
# 6. CONNECTED COMPONENTS (DFS)
# =====================================================

print("\nFinding connected components...")

components = list(nx.connected_components(G))

sizes = [len(c) for c in components]

print("Number of connected components:", len(components))
print("Largest component size:", max(sizes))

# Save components

component_results = []

for idx, comp in enumerate(components):

    for gene in comp:

        component_results.append(
            [gene, idx]
        )

component_df = pd.DataFrame(
    component_results,
    columns=["Gene", "Component"]
)

component_df.to_csv(
    "connected_components.csv",
    index=False
)

# =====================================================
# 7. LOUVAIN COMMUNITY DETECTION
# =====================================================

print("\nRunning Louvain...")

partition = community_louvain.best_partition(G)

n_modules = len(set(partition.values()))

print("Number of modules:", n_modules)

module_df = pd.DataFrame(
    partition.items(),
    columns=["Gene", "Module"]
)

module_df.to_csv(
    "louvain_modules.csv",
    index=False
)

# =====================================================
# 8. HUB GENES
# =====================================================

print("\nComputing centrality...")

degree_centrality = nx.degree_centrality(G)

hub_df = pd.DataFrame(
    degree_centrality.items(),
    columns=["Gene", "DegreeCentrality"]
)

hub_df = hub_df.sort_values(
    by="DegreeCentrality",
    ascending=False
)

hub_df.to_csv(
    "hub_genes.csv",
    index=False
)

print("\nTOP 20 HUB GENES")
print("----------------")

print(hub_df.head(20))

# =====================================================
# 9. NETWORK SUMMARY
# =====================================================

density = nx.density(G)

print("\nNETWORK SUMMARY")
print("----------------")
print("Nodes:", G.number_of_nodes())
print("Edges:", G.number_of_edges())
print("Density:", density)
print("Connected Components:", len(components))
print("Modules:", n_modules)

# =====================================================
# 10. SAVE EDGE LIST
# =====================================================

edge_list = []

for u, v, d in G.edges(data=True):

    edge_list.append([
        u,
        v,
        d["weight"]
    ])

edges_df = pd.DataFrame(
    edge_list,
    columns=[
        "Gene1",
        "Gene2",
        "Correlation"
    ]
)

edges_df.to_csv(
    "coexpression_network.csv",
    index=False
)

modules = pd.read_csv("louvain_modules.csv")

module_dict = dict(
    zip(modules["Gene"], modules["Module"])
)

for node in G.nodes():
    if node in module_dict:
        G.nodes[node]["module"] = module_dict[node]

nx.write_gexf(
    G,
    "gene_network_modules.gexf"
)

print("\nFiles generated:")
print("---------------")
print("coexpression_network.csv")
print("connected_components.csv")
print("louvain_modules.csv")
print("hub_genes.csv")

print("\nAnalysis completed successfully.")

