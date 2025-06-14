import scanpy as sc

def load_and_preprocess(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    return adata

def apply_pca_umap_clustering(adata):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return adata


def compute_top_marker_genes(adata, n_genes=5):
    if 'leiden' not in adata.obs:
        raise ValueError("Leiden clusters not found.")

    sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
    result = sc.get.rank_genes_groups_df(adata, group=None)

    # Get top N genes per cluster
    top_markers = (
        result.groupby("group")
        .head(n_genes)
        .reset_index(drop=True)
        .rename(columns={"group": "Cluster", "names": "Gene", "logfoldchanges": "LogFC", "pvals_adj": "Adjusted_pval"})
    )
    return top_markers
