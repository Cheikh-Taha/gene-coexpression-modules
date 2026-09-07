import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering


df = pd.read_csv(
    "GSE60424_norm_counts_TPM_GRCh38.p13_NCBI.tsv",
    sep="\t"
)

df = df.set_index("GeneID")

print(df.shape)

df = df[df.mean(axis=1) > 1]

variance = df.var(axis=1)

top_genes = variance.nlargest(1000).index

expr = df.loc[top_genes]

corr_matrix = expr.T.corr()

similarity = (corr_matrix + 1) / 2

model = SpectralClustering(
    n_clusters=10,
    affinity='precomputed',
    assign_labels='kmeans',
    random_state=42
)

labels = model.fit_predict(similarity)

result = pd.DataFrame({
    "GeneID": expr.index,
    "Cluster": labels
})

print(result.head())

result.to_csv("spectral_clusters.csv", index=False)

print(result["Cluster"].value_counts())

import matplotlib.pyplot as plt

result["Cluster"].value_counts().plot(kind="bar")

plt.title("Spectral Clustering Modules")
plt.xlabel("Cluster")
plt.ylabel("Number of genes")
plt.show()