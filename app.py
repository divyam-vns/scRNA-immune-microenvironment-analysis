import streamlit as st
import scanpy as sc

st.title("🧬 Immune Microenvironment Explorer")

st.write("Loading PBMC3k dataset...")

adata = sc.datasets.pbmc3k()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

adata = adata[:, adata.var.highly_variable]

sc.pp.scale(adata)
sc.tl.pca(adata)

sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.5)

st.subheader("UMAP Clusters")

fig = sc.pl.umap(
    adata,
    color="leiden",
    return_fig=True,
    show=False
)

st.pyplot(fig)
