import streamlit as st
import pandas as pd
import scanpy as sc

st.title("🧬 Immune Microenvironment Explorer + Cell Communication")

# -------------------------
# LOAD DATA
# -------------------------
adata = sc.datasets.pbmc3k()

adata.obs["celltype"] = sc.tl.leiden(adata, resolution=0.5, copy=True)

liana_df = pd.read_csv("results/liana_cellcell_interactions.csv")

# -------------------------
# SIDEBAR
# -------------------------
view = st.sidebar.selectbox(
    "Select View",
    ["UMAP Clusters", "Gene Expression", "Cell Counts", "Cell Communication"]
)

# -------------------------
# UMAP
# -------------------------
if view == "UMAP Clusters":
    st.subheader("Immune Cell Clusters")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    fig = sc.pl.umap(adata, color="leiden", return_fig=True, show=False)
    st.pyplot(fig)

# -------------------------
# GENE EXPRESSION
# -------------------------
elif view == "Gene Expression":
    gene = st.text_input("Enter gene (CD3D, MS4A1, LYZ, NKG7)")

    if gene:
        if gene in adata.var_names:
            fig = sc.pl.umap(adata, color=gene, return_fig=True, show=False)
            st.pyplot(fig)
        else:
            st.warning("Gene not found")

# -------------------------
# CELL COUNTS
# -------------------------
elif view == "Cell Counts":
    st.write(adata.obs["celltype"].value_counts())

# -------------------------
# CELL-CELL COMMUNICATION (NEW)
# -------------------------
elif view == "Cell Communication":
    st.subheader("Ligand-Receptor Interactions (LIANA)")

    st.dataframe(liana_df.head(20))

    sender = st.selectbox("Sender cluster", sorted(liana_df["source"].unique()))
    receiver = st.selectbox("Receiver cluster", sorted(liana_df["target"].unique()))

    filtered = liana_df[
        (liana_df["source"] == sender) &
        (liana_df["target"] == receiver)
    ].sort_values("lrscore", ascending=False)

    st.write("Top interactions")
    st.dataframe(filtered.head(10))
