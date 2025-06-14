import scanpy as sc
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd

def plot_umap(adata):
    label_column = get_best_label_column(adata)

    # Ensure UMAP is already computed
    if 'X_umap' not in adata.obsm:
        raise ValueError("UMAP has not been computed yet.")

    # Get UMAP coordinates
    umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    umap_df[label_column] = adata.obs[label_column].values

    # Optional: add cell names to hover
    umap_df['Cell'] = umap_df.index

    # Create interactive Plotly scatter plot
    fig = px.scatter(
        umap_df,
        x='UMAP1',
        y='UMAP2',
        color=label_column,
        hover_data=['Cell', label_column],
        title=f"Interactive UMAP: colored by {label_column}",
        height=600
    )

    st.plotly_chart(fig, use_container_width=True)



def plot_gene_heatmap(adata, gene):
    if gene in adata.var_names:
        sc.pl.matrixplot(
            adata,
            var_names=[gene],
            groupby="leiden",
            cmap="viridis",
            standard_scale="var",  # normalize per gene
            use_raw=False,
            show=False
        )
        st.pyplot(plt.gcf())
    else:
        st.warning("Gene not found.")



def get_best_label_column(adata):
    # Priority list of known label columns
    preferred_labels = ['cell_type', 'bulk_labels', 'annotation', 'labels', 'ident']

    for col in preferred_labels:
        if col in adata.obs.columns:
            return col  # Found a real label column

    # If no label column exists, generate human-friendly labels from 'leiden'
    if 'leiden' in adata.obs.columns:
        adata.obs['cluster_label'] = adata.obs['leiden'].apply(lambda x: f"Cluster {x}")
        return 'cluster_label'

    return None  # fallback: should not happen
