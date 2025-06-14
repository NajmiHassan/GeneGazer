import streamlit as st
from visualizer import plot_umap, plot_gene_heatmap, get_best_label_column
from sc_processing import compute_top_marker_genes

def render_dataset_selector():
    all_datasets = st.session_state.get('all_datasets', [])
    if not all_datasets:
        st.warning("No datasets uploaded yet. Please upload one in the 'Load Data' tab.")
        return None

    dataset_names = [ds['label'] for ds in all_datasets]
    selected_label = st.selectbox("📁 Select a dataset to view:", dataset_names)

    selected_dataset = next(ds for ds in all_datasets if ds['label'] == selected_label)
    st.session_state['adata'] = selected_dataset['adata']
    return selected_dataset['adata']

def render_umap_section(adata):
    st.subheader("🔍 Interactive UMAP")

    # Sorted full list of genes
    gene_options = sorted(list(adata.var_names))

    hover_gene = st.selectbox(
        "🔎 Type or select a gene to highlight on hover (optional):",
        options=[""] + gene_options  # empty for "None"
    )

    plot_umap(adata, hover_gene if hover_gene else None)


def render_gene_heatmap_section(adata):
    st.subheader("🎯 Gene Heatmap (by Cluster)")

    # Sorted full list of genes
    gene_options = sorted(list(adata.var_names))

    # Selectbox with search enabled
    gene = st.selectbox(
        "🔎 Type or select a gene to plot:",
        options=gene_options
    )

    plot_gene_heatmap(adata, gene)


def render_marker_gene_table(adata):
    with st.expander("📌 Show Top Marker Genes per Cluster (Auto-Detected)"):
        try:
            top_genes = compute_top_marker_genes(adata)
            st.dataframe(top_genes)
        except Exception as e:
            st.warning(f"Could not compute markers: {str(e)}")

def render_visualization():
    st.title("📊 Visualizations")

    adata = render_dataset_selector()
    if adata is None:
        return

    label_column = get_best_label_column(adata)
    st.caption(f"UMAP is colored by: **{label_column}**")

    render_umap_section(adata)
    render_marker_gene_table(adata)
    render_gene_heatmap_section(adata)
