
# 🧬 scRNA-seq Immune Microenvironment Analysis (PBMC)

An end-to-end single-cell RNA-seq pipeline to study immune cell heterogeneity and intercellular communication using PBMC data.

---

## 📊 Project Overview

This project analyzes immune microenvironment structure using scRNA-seq data and integrates:

- Cell clustering (Leiden)
- Differential gene expression analysis
- Ligand–receptor communication (LIANA)
- GO and KEGG pathway enrichment
- Interactive Streamlit dashboard

---

## 🔬 Workflow

### 1. Data Acquisition
- PBMC scRNA-seq dataset loaded using Scanpy
- Quality control filtering applied

### 2. Preprocessing
- Normalization (log1p)
- Highly variable gene selection
- PCA dimensionality reduction
- Neighborhood graph construction

### 3. Clustering & Visualization
- UMAP embedding
- Leiden clustering for immune cell populations

### 4. Differential Expression
- Rank-based gene markers per cluster
- Identification of immune-specific signatures

### 5. Cell–Cell Communication
- Ligand–receptor inference using LIANA
- Immune signaling network reconstruction

### 6. Pathway Enrichment
- GO Biological Process enrichment
- KEGG pathway analysis
- Identification of immune functional programs

---

##  Key Findings

- Distinct immune cell populations identified (T cells, B cells, monocytes, dendritic cells)
- Strong immune activation signatures observed in T-cell clusters
- Key signaling pathways include:
  - T cell receptor signaling
  - Cytokine–cytokine receptor interaction
  - Antigen processing and presentation

---

##  Tech Stack

- Python
- Scanpy
- AnnData
- LIANA
- GSEApy
- Pandas
- Streamlit
- GitHub

---

##  Streamlit App

Interactive dashboard includes:

- UMAP visualization of immune clusters
- Gene expression explorer
- Cell type composition
- Ligand–receptor interaction viewer
- GO/KEGG pathway enrichment explorer

---

##  Repository Structure

##  Repository Structure


scRNA-immune-microenvironment-analysis/
│
├── app.py
├── pbmc_analysis.py
├── DE_analysis.py
│
├── data/
│ └── adata_processed.h5ad
│
├── results/
│ ├── go_enrichment_results.csv
│ └── kegg_enrichment_results.csv
│
├── figures/
│ ├── umap_leiden.png
│ ├── umap_celltypes.png
│ └── pathways/
│ ├── go_barplot.png
│ └── kegg_barplot.png
│
└── README.md


---

## How to Run

```bash
pip install -r requirements.txt
streamlit run app.py
 Biological Significance

This project demonstrates how single-cell transcriptomics can be used to:

Map immune cell diversity
Infer functional cell states
Reconstruct intercellular signaling networks
Identify immune pathway activation patterns
👨‍💻 Author

Dr. Divya Mishra, Ph.D.
Bioinformatics & Computational Immunology Project
Built using Scanpy + LIANA + Streamlit
