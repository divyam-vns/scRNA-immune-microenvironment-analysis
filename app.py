
import streamlit as st
import scanpy as sc
import pandas as pd
import os

st.set_page_config(page_title="Immune scRNA-seq Explorer", layout="wide")

st.title("🧬 scRNA-seq Immune Microenvironment Dashboard")

# -----------------------------
# SAFE DATA LOADING
# -----------------------------
@st.cache_data
def load_data():
    if os.path.exists("data/adata_processed.h5ad"):
        adata = sc.read_h5ad("data/adata_processed.h5ad")
    else:
        st.warning("adata file not found, loading demo PBMC dataset")
        adata = sc.datasets.pbmc3k()
    return adata

adata = load_data()

# -----------------------------
# SIDEBAR NAVIGATION
# -----------------------------
view = st.sidebar.selectbox(
    "Select Analysis View",
    [
        "Overview",
        "UMAP",
        "Gene Expression",
        "Cell Types",
        "Pathway Enrichment"
    ]
)

# -----------------------------
# OVERVIEW
# -----------------------------
if view == "Overview":
    st.subheader("Dataset Overview")
    st.write(adata)

# -----------------------------
# UMAP
# -----------------------------
elif view == "UMAP":
    st.subheader("UMAP Clusters")

    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color="leiden", show=False)
        st.pyplot()
    else:
        st.warning("UMAP not available")

# -----------------------------
# GENE EXPRESSION
# -----------------------------
elif view == "Gene Expression":
    st.subheader("Gene Expression Viewer")

    gene = st.text_input("Enter Gene Name", "CD3D")

    if gene in adata.var_names:
        sc.pl.umap(adata, color=gene, show=False)
        st.pyplot()
    else:
        st.error(f"{gene} not found in dataset")

# -----------------------------
# CELL TYPES
# -----------------------------
elif view == "Cell Types":
    st.subheader("Cluster Composition")

    if "leiden" in adata.obs:
        st.bar_chart(adata.obs["leiden"].value_counts())
    else:
        st.warning("No clustering found")

# -----------------------------
# PATHWAYS
# -----------------------------
elif view == "Pathway Enrichment":
    st.subheader("GO / KEGG Pathways")

    if os.path.exists("results/go_enrichment_results.csv"):
        go = pd.read_csv("results/go_enrichment_results.csv")
        st.write("GO Top Results")
        st.dataframe(go.head(10))
    else:
        st.warning("GO results not found")

    if os.path.exists("results/kegg_enrichment_results.csv"):
        kegg = pd.read_csv("results/kegg_enrichment_results.csv")
        st.write("KEGG Top Results")
        st.dataframe(kegg.head(10))
    else:
        st.warning("KEGG results not found")
