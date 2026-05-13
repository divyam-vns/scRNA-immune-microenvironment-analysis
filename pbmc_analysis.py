
import scanpy as sc

# Load PBMC dataset
adata = sc.datasets.pbmc3k()

# Quality control
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

adata = adata[:, adata.var.highly_variable]

# Scaling + PCA
sc.pp.scale(adata)
sc.tl.pca(adata)

# Neighbors + UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# Leiden clustering
sc.tl.leiden(adata, resolution=0.5)

# Save processed dataset
adata.write("data/adata_processed.h5ad")

print("Pipeline completed successfully")
