import streamlit as st
import scanpy as sc

st.title("Immune Microenvironment Explorer")

adata = sc.datasets.pbmc3k()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVGs (DO NOT subset full data)
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

# PCA uses HVGs only (temporary view)
adata_hvg = adata[:, adata.var.highly_variable]

sc.pp.scale(adata_hvg)
sc.tl.pca(adata_hvg)
sc.pp.neighbors(adata_hvg)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.5)

# store embeddings back into full object
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
adata.obs["leiden"] = adata_hvg.obs["leiden"]

# sidebar
view = st.sidebar.selectbox(
    "Select View",
    ["UMAP Clusters", "Gene Expression", "Cell Counts"]
)

# UMAP
if view == "UMAP Clusters":
    st.subheader("Leiden Clusters")
    fig = sc.pl.umap(adata, color="leiden", return_fig=True, show=False)
    st.pyplot(fig)

# Gene expression (NOW WORKS FOR ALL GENES)
elif view == "Gene Expression":
    gene = st.text_input("Enter gene (CD3D, MS4A1, LYZ, NKG7)")

    if gene:
        if gene in adata.var_names:
            fig = sc.pl.umap(adata, color=gene, return_fig=True, show=False)
            st.pyplot(fig)
        else:
            st.warning("Gene not found")

# counts
elif view == "Cell Counts":
    st.write(adata.obs["leiden"].value_counts())
