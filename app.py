
import streamlit as st
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt

st.set_page_config(page_title="Immune scRNA-seq Explorer", layout="wide")

st.title("🧬 Immune Microenvironment scRNA-seq Explorer")

# -----------------------------
# DATA LOADING (ROBUST)
# -----------------------------
@st.cache_data
def load_data():
    if os.path.exists("data/adata_processed.h5ad"):
        adata = sc.read_h5ad("data/adata_processed.h5ad")
    else:
        adata = sc.datasets.pbmc3k()

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        adata = adata[:, adata.var.highly_variable].copy()

        sc.pp.scale(adata)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.5)

        sc.tl.rank_genes_groups(adata, groupby="leiden")

    return adata


adata = load_data()

# -----------------------------
# GET TOP MARKER GENES (SAFE)
# -----------------------------
def get_top_markers(adata, n=10):
    try:
        markers = adata.uns["rank_genes_groups"]["names"]
        df = pd.DataFrame(markers)
        return list(df.iloc[:n, 0].values)
    except:
        return ["CD3D", "MS4A1", "LYZ", "NKG7"]

top_genes = get_top_markers(adata)

# -----------------------------
# TABS (NEW UX)
# -----------------------------
tab1, tab2, tab3, tab4 = st.tabs([
    "📊 Overview",
    "🧬 UMAP",
    "🔬 Gene Expression",
    "🧪 Pathways"
])

# =============================
# TAB 1: OVERVIEW (IMPROVED)
# =============================
with tab1:
    st.subheader("Dataset Overview")

    col1, col2, col3 = st.columns(3)

    col1.metric("Cells", adata.n_obs)
    col2.metric("Genes", adata.n_vars)
    col3.metric("Clusters", len(adata.obs["leiden"].unique()) if "leiden" in adata.obs else "NA")

    st.write("### Cluster Distribution")
    if "leiden" in adata.obs:
        st.bar_chart(adata.obs["leiden"].value_counts())

    st.write("### Top Marker Genes")
    st.write(top_genes)

# =============================
# TAB 2: UMAP
# =============================
with tab2:
    st.subheader("UMAP Clusters")

    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color="leiden", show=False)
        st.pyplot()
    else:
        st.warning("UMAP not available")

# =============================
# TAB 3: GENE EXPRESSION (FIXED)
# =============================
with tab3:
    st.subheader("Gene Expression Viewer")

    gene = st.text_input("Enter gene name", "CD3D")

    if gene in adata.var_names:
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=gene, ax=ax, show=False)
        st.pyplot(fig)
    else:
        st.error(f"{gene} not found")

    st.write("### Suggested Marker Genes")
    st.write(top_genes)

# =============================
# TAB 4: PATHWAYS
# =============================
with tab4:
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
