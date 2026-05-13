
import streamlit as st
import scanpy as sc
import pandas as pd

st.set_page_config(page_title="Immune Microenvironment Explorer", layout="wide")

st.title("🧬 Immune Microenvironment Explorer (scRNA-seq + LIANA)")

# -------------------------
# LOAD DATA (SAFE OPTION)
# -------------------------
@st.cache_data
def load_data():
    try:
        adata = sc.read("data/adata_processed.h5ad")
    except:
        adata = sc.datasets.pbmc3k()

        # minimal preprocessing (safe pipeline)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata)

    adata.obs["celltype"] = adata.obs["leiden"]
    return adata

adata = load_data()

# -------------------------
# LOAD LIANA RESULTS
# -------------------------
@st.cache_data
def load_liana():
    return pd.read_csv("results/liana_cellcell_interactions.csv")

liana_df = load_liana()

# -------------------------
# SIDEBAR NAVIGATION
# -------------------------
view = st.sidebar.selectbox(
    "Select View",
    ["UMAP Clusters", "Gene Expression", "Cell Counts", "Cell Communication"]
)

# -------------------------
# UMAP VIEW
# -------------------------
if view == "UMAP Clusters":
    st.subheader("Immune Cell Clusters (UMAP)")

    fig = sc.pl.umap(adata, color="celltype", return_fig=True, show=False)
    st.pyplot(fig)

# -------------------------
# GENE EXPRESSION
# -------------------------
elif view == "Gene Expression":
    st.subheader("Gene Expression Viewer")

    gene = st.text_input("Enter gene (e.g., CD3D, MS4A1, NKG7, LYZ)")

    if gene:
        if gene in adata.var_names:
            fig = sc.pl.umap(adata, color=gene, return_fig=True, show=False)
            st.pyplot(fig)
        else:
            st.warning("Gene not found in dataset")

# -------------------------
# CELL COUNTS
# -------------------------
elif view == "Cell Counts":
    st.subheader("Cluster Distribution")
    st.write(adata.obs["celltype"].value_counts())

# -------------------------
# CELL COMMUNICATION (LIANA)
# -------------------------
elif view == "Cell Communication":
    st.subheader("🧬 Immune Signaling Network (LIANA)")

    st.dataframe(liana_df.head(20))

    sender = st.selectbox("Sender cluster", sorted(liana_df["source"].unique()))
    receiver = st.selectbox("Receiver cluster", sorted(liana_df["target"].unique()))

    filtered = liana_df[
        (liana_df["source"] == sender) &
        (liana_df["target"] == receiver)
    ].sort_values("lrscore", ascending=False)

    st.write("Top ligand–receptor interactions")
    st.dataframe(filtered.head(15))
