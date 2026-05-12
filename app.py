import streamlit as st
import scanpy as sc

st.title("🧬 Immune Microenvironment Explorer (PBMC scRNA-seq)")

# Load processed dataset
adata = sc.read("data/adata_processed.h5ad")

st.write("Dataset loaded:", adata)

# Sidebar selection
cell_types = adata.obs["cell_type"].unique().tolist()
selected = st.sidebar.selectbox("Select Cell Type", cell_types)

# Filter
subset = adata[adata.obs["cell_type"] == selected]

st.subheader(f"UMAP: {selected}")

sc.pl.umap(subset, color="cell_type", show=False)
st.pyplot()
