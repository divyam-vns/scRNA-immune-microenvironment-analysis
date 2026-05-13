
import streamlit as st
import scanpy as sc
import pandas as pd
import os

st.set_page_config(page_title="scRNA Immune Dashboard", layout="wide")

st.title("🧬 scRNA-seq Immune Microenvironment Dashboard")

# -----------------------------
# SAFE DATA LOADING + FULL PIPELINE
# -----------------------------
@st.cache_data
def load_data():
    if os.path.exists("data/adata_processed.h5ad"):
        adata = sc.read_h5ad("data/adata_processed.h5ad")
        return adata

    # fallback dataset
    adata = sc.datasets.pbmc3k()

    # FULL SAFE PREPROCESSING PIPELINE
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()

    sc.pp.scale(adata)
    sc.tl.pca(adata)

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    return adata


adata = load_data()

# -----------------------------
# SIDEBAR NAVIGATION
# -----------------------------
view = st.sidebar.selectbox(
    "Select View",
    ["Overview", "UMAP", "Gene Expression", "Cell Composition", "Pathways"]
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
    st.subheader("UMAP Clustering")

    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color="leiden", show=False)
        st.pyplot()
    else:
        st.error("UMAP not available")

# -----------------------------
# GENE EXPRESSION (SAFE)
# -----------------------------
elif view == "Gene Expression":
    st.subheader("Gene Expression Viewer")

    gene = st.text_input("Enter gene", "CD3D")

    if gene in adata.var_names:
        sc.pl.umap(adata, color=gene, show=False)
        st.pyplot()
    else:
        st.warning(f"{gene} not found in dataset")

# -----------------------------
# CELL COMPOSITION
# -----------------------------
elif view == "Cell Composition":
    st.subheader("Cluster Composition")

    if "leiden" in adata.obs:
        st.bar_chart(adata.obs["leiden"].value_counts())
    else:
        st.warning("No clustering available")

# -----------------------------
# PATHWAYS (SAFE CSV HANDLING)
# -----------------------------
elif view == "Pathways":
    st.subheader("GO / KEGG Enrichment")

    go_path = "results/go_enrichment_results.csv"
    kegg_path = "results/kegg_enrichment_results.csv"

    if os.path.exists(go_path):
        go = pd.read_csv(go_path)
        st.write("GO Enrichment")
        st.dataframe(go.head(10))
    else:
        st.warning("GO results missing")

    if os.path.exists(kegg_path):
        kegg = pd.read_csv(kegg_path)
        st.write("KEGG Enrichment")
        st.dataframe(kegg.head(10))
    else:
        st.warning("KEGG results missing")
