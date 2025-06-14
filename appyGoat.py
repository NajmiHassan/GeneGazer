import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile

st.set_page_config(page_title="RNA-seq Viewer", layout="wide")

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

st.sidebar.title("Navigation")
menu = st.sidebar.radio("Select a section:", ["📘 Instructions", "📁 Load Data", "📊 Visualize", "🤖 AI Assistant"])

if menu == "📘 Instructions":
    st.title("Single-Cell RNA-seq Viewer")
    st.markdown("""
    This tool lets you explore single-cell RNA-seq data with ease:

    1. Go to **Load Data** and upload a `.h5ad` file or use a sample dataset.
    2. Go to **Visualize** to generate interactive plots:
       - Clustering (UMAP)
       - Gene heatmap
    3. Use the **AI Assistant** to ask questions.

    **Accepted format**:
    `.h5ad` files exported from Scanpy or Seurat.
    """)
    st.image("https://i.imgur.com/zVfGZkP.png", caption="Example of UMAP visualization")

elif menu == "📁 Load Data":
    st.title("Load RNA-seq Data")

    uploaded_file = st.file_uploader("Upload a `.h5ad` file", type=["h5ad"])
    if uploaded_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            tmp.write(uploaded_file.read())
            temp_path = tmp.name
        try:
            adata = sc.read_h5ad(temp_path)
            adata = load_and_preprocess(adata)
            adata = apply_pca_umap_clustering(adata)
            st.session_state['adata'] = adata
            st.success("File successfully loaded and processed!")
        except Exception as e:
            st.error(f"Error during processing: {str(e)}")

    st.markdown("### Or choose one of the public datasets:")
    col1, col2 = st.columns(2)

    with col1:
        if st.button("Use PBMC 3k (Scanpy)"):
            try:
                adata = sc.datasets.pbmc3k()
                adata = load_and_preprocess(adata)
                adata = apply_pca_umap_clustering(adata)
                st.session_state['adata'] = adata
                st.success("PBMC 3k successfully loaded and processed!")
            except Exception as e:
                st.error(f"Error loading PBMC 3k: {str(e)}")

    with col2:
        if st.button("Use Mouse Brain Cortex (Scanpy)"):
            try:
                adata = sc.datasets.mouse_brain_cortex_1k()
                adata = load_and_preprocess(adata)
                adata = apply_pca_umap_clustering(adata)
                st.session_state['adata'] = adata
                st.success("Mouse brain cortex successfully loaded and processed!")
            except Exception as e:
                st.error(f"Error loading mouse brain cortex: {str(e)}")

elif menu == "📊 Visualize":
    st.title("Visualizations")

    if 'adata' not in st.session_state:
        st.warning("Please load a file from the 'Load Data' tab.")
    else:
        adata = st.session_state['adata']
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("UMAP Visualization")
            sc.pl.umap(adata, color='leiden', show=False)
            st.pyplot(plt.gcf())

        with col2:
            st.subheader("Gene Heatmap")
            gene = st.text_input("Enter a gene name:", "IL7R")
            if gene in adata.var_names:
                sc.pl.matrixplot(adata, var_names=[gene], groupby="leiden", cmap="viridis", use_raw=False, show=False)
                st.pyplot(plt.gcf())
            else:
                st.error("Gene not found.")

elif menu == "🤖 AI Assistant":
    st.title("AI Assistant: Ask about RNA-seq")

    with st.chat_message("assistant"):
        st.markdown("Hi! Ask anything about RNA-seq, genes, or how to use this tool!")

    question = st.chat_input("Type your question here...")
    if question:
        with st.chat_message("user"):
            st.markdown(question)

        with st.chat_message("assistant"):
            if "gene" in question.lower():
                st.markdown("A gene is a DNA sequence that contains instructions to produce proteins.")
            elif "umap" in question.lower():
                st.markdown("UMAP is a technique that reduces complex data into 2D for easier visualization.")
            elif "h5ad" in question.lower():
                st.markdown("`.h5ad` is a file format used by Scanpy to store gene expression data.")
            else:
                st.markdown("I'm a simple assistant, but I can help with basic RNA-seq questions.")

st.markdown("---")
st.info("Let's shine in this hackathon, team Goat!")

