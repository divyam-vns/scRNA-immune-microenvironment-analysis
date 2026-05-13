
import scanpy as sc

def run_de_analysis(adata):
    """
    Differential expression analysis using Scanpy
    """

    # Rank genes per cluster
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon"
    )

    # Save DE plot
    sc.pl.rank_genes_groups(
        adata,
        n_genes=20,
        sharey=False,
        show=False,
        save="_leiden_de_markers.png"
    )

    return adata
