
# Single-Cell RNA-seq Immune Microenvironment Analysis

## Project Overview

This project performs single-cell RNA sequencing (scRNA-seq) analysis of immune cells using the PBMC3k dataset from Scanpy.

The workflow includes:

- Quality control (QC)
- Normalization
- Highly variable gene selection
- PCA dimensionality reduction
- UMAP visualization
- Leiden clustering
- Immune cell type annotation
- Interactive Streamlit visualization

This project demonstrates a complete computational pipeline for immune microenvironment profiling using single-cell transcriptomics.

---

# Dataset

- PBMC3k dataset
- Source: Scanpy built-in dataset
- ~2700 peripheral blood mononuclear cells (PBMCs)

---

# Methods

## Preprocessing
- Cell filtering
- Gene filtering
- Library size normalization
- Log transformation

## Feature Engineering
- Highly variable gene (HVG) detection

## Dimensionality Reduction
- Principal Component Analysis (PCA)
- UMAP embedding

## Clustering
- Leiden community detection algorithm

## Immune Annotation
Marker genes used:

| Marker | Cell Type |
|---|---|
| CD3D | T cells |
| MS4A1 | B cells |
| NKG7 | NK cells |
| LYZ | Monocytes |

---

# Results

## UMAP Clustering

![Leiden UMAP](figures/umap_leiden.png)

---

## Immune Cell Type Annotation

![Cell Types](figures/umap_celltypes.png)

---

# Streamlit App

Interactive visualization app built using Streamlit.

Features:
- UMAP visualization
- Cluster exploration
- Immune cell inspection

---

# Technologies Used

- Python
- Scanpy
- AnnData
- Streamlit
- NumPy
- Pandas
- Matplotlib

---

# Future Improvements

- Differential expression analysis
- Cell-cell communication analysis
- Pathway enrichment
- Integration of multiple datasets
- Spatial transcriptomics extension

---

# 👤 Author

Divya Mishra

Bioinformatics | Computational Biology | Genomics | Translational Data Science
