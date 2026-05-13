import streamlit as st
import scanpy as sc

st.title("🧬 Immune Microenvironment Explorer (PBMC scRNA-seq)")

adata = sc.datasets.pbmc3k()

# preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# sidebar navigation
view = st.sidebar.selectbox(
    "Select View",
    ["UMAP Clusters", "Gene Expression", "Cell Counts"]
)

# -----------------------
# VIEW 1: UMAP
# -----------------------
if view == "UMAP Clusters":
    st.subheader("Leiden Clusters")

    fig = sc.pl.umap(
        adata,
        color="leiden",
        return_fig=True,
        show=False
    )
    st.pyplot(fig)

# -----------------------
# VIEW 2: Gene Expression
# -----------------------
elif view == "Gene Expression":
    gene = st.text_input("Enter gene name (e.g., CD3D, MS4A1, LYZ, NKG7)")

    if gene:
        if gene in adata.var_names:
            fig = sc.pl.umap(
                adata,
                color=gene,
                return_fig=True,
                show=False
            )
            st.pyplot(fig)
        else:
            st.warning("Gene not found in dataset")

# -----------------------
# VIEW 3: Cell Counts
# -----------------------
elif view == "Cell Counts":
    st.write(adata.obs["leiden"].value_counts())
