# gene-coexpression-modules
Detection of gene co-expression modules using graph theory and unsupervised machine learning on the GSE60424 RNA-seq dataset.
# Gene Co-Expression Module Detection

##  Overview

This project focuses on the detection of functional modules in a human gene co-expression network using RNA-seq expression data.

The project combines **bioinformatics, graph theory, and artificial intelligence** to identify groups of genes with similar expression patterns and characterize highly connected genes within the network.

The analysis is based on the **GSE60424** dataset from the Gene Expression Omnibus (GEO).

---

##  Objectives

The main objectives are:

- Preprocess and filter RNA-seq gene expression data
- Identify highly variable genes
- Calculate gene-gene Pearson correlations
- Construct a weighted gene co-expression network
- Analyze network connectivity using graph theory
- Identify connected components using DFS
- Detect communities/modules using the Louvain algorithm
- Identify hub genes using degree centrality
- Apply Spectral Clustering as an unsupervised machine learning approach
- Visualize the resulting gene network and modules

---

##  Dataset

**Dataset:** GSE60424

**Data type:** RNA-seq normalized expression data (TPM)

The dataset contains expression profiles from human immune-related samples.

After preprocessing, the project selects the **1000 most variable genes** from the original dataset for network construction.

Final expression matrix:


1000 genes × 134 samples

## Methodology

The project follows this pipeline:
```text
RNA-seq Expression Data
          ↓
Data Preprocessing
          ↓
Gene Filtering
          ↓
Variance Analysis
          ↓
Top Variable Genes
          ↓
Pearson Correlation
          ↓
Co-Expression Network
          ↓
Graph Analysis
     ↙          ↘
   DFS        Louvain
     ↓          ↓
Connected    Modules
Components
              ↓
         Hub Genes
              ↓
     Spectral Clustering
              ↓
       AI-based Modules


