import scanpy as sc
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd
import plotly.express as px
import pandas as pd
import streamlit as st

def plot_umap(adata, gene_for_hover=None):
    label_column = get_best_label_column(adata)

    if 'X_umap' not in adata.obsm:
        raise ValueError("UMAP coordinates not found. Make sure UMAP is computed.")

    # Create DataFrame with UMAP + cluster labels
    df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    df[label_column] = adata.obs[label_column].astype(str)
    df['Cell'] = adata.obs_names

    # Add gene expression values (optional)
    if gene_for_hover and gene_for_hover in adata.var_names:
        df[f'{gene_for_hover}_expr'] = adata[:, gene_for_hover].X.toarray().flatten()
    elif gene_for_hover:
        st.warning(f"Gene `{gene_for_hover}` not found in dataset.")

    # Define hover columns
    hover_cols = ['Cell', label_column]
    if gene_for_hover and f'{gene_for_hover}_expr' in df.columns:
        hover_cols.append(f'{gene_for_hover}_expr')

    # Plot with Plotly
    fig = px.scatter(
        df,
        x='UMAP1',
        y='UMAP2',
        color=label_column,
        hover_data=hover_cols,
        title=f"UMAP colored by: {label_column}",
        height=600
    )

    st.plotly_chart(fig, use_container_width=True)


def plot_gene_heatmap(adata, gene):
    if gene in adata.var_names:
        # Clear previous plot to avoid overlap
        plt.clf()

        # Create heatmap with smaller figure size
        sc.pl.matrixplot(
            adata,
            var_names=[gene],
            groupby="leiden",
            cmap="viridis",
            use_raw=False,
            standard_scale="var",
            dendrogram=False,
            figsize=(6, 2),  # ⬅️ smaller width & height to fit screen
            show=False
        )
        st.pyplot(plt.gcf())
    else:
        st.warning("Gene not found in dataset.")

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
