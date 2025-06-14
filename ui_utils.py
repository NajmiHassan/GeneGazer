import streamlit as st
import tempfile
from file_handler import detect_file_type, load_file_dynamic, get_csv_download_link
from sc_processing import load_and_preprocess, apply_pca_umap_clustering
from visualization_ui import render_visualization
import scanpy as sc
import streamlit as st

def render_instructions():
    st.title("Single-Cell RNA-seq Viewer")
    st.markdown("""
This tool lets you explore single-cell RNA-seq data with ease:

1. Go to **Load Data** and upload any of the following:
   - `.h5ad` file (Scanpy)
   - `.csv` gene expression matrix (cells in columns, genes in rows)
   - `.mtx` + `genes.tsv` + `barcodes.tsv` (10x Genomics)
2. Go to **Visualize** to generate:
   - UMAP plots
   - Gene heatmaps
3. Use the **AI Assistant** to ask about RNA-seq.

Supported formats: `.h5ad`, `.csv`, `.mtx + genes.tsv + barcodes.tsv`
    """)
    st.markdown(get_csv_download_link(), unsafe_allow_html=True)
    st.image("https://i.imgur.com/zVfGZkP.png", caption="Example of UMAP visualization")


def render_load_data():
    st.title("Load RNA-seq Data")
    uploaded_file = st.file_uploader("Upload a `.h5ad`, `.csv`, or `.mtx` file", type=["h5ad", "csv", "mtx"])

    extra = {}
    if uploaded_file and uploaded_file.name.endswith("mtx"):
        st.info("📎 Please upload `genes.tsv` and `barcodes.tsv` to complete `.mtx` input.")
        extra["genes"] = st.file_uploader("Upload gene list (.mtx_rows or genes.tsv)", type=["tsv", "txt", "mtx_rows"], key="genes_file")
        extra["barcodes"] =st.file_uploader("Upload barcode list (.mtx_cols or barcodes.tsv)", type=["tsv", "txt", "mtx_cols"], key="barcodes_file")

    if uploaded_file:
        file_type = detect_file_type(uploaded_file.name)

        if file_type not in ["h5ad", "csv", "mtx"]:
            st.error("Unsupported file format.")
            return

        with tempfile.NamedTemporaryFile(delete=False, suffix=f".{file_type}") as tmp:
            tmp.write(uploaded_file.read())
            file_path = tmp.name

        if file_type == "mtx" and (not extra["genes"] or not extra["barcodes"]):
            st.warning("Please upload both `genes.tsv` and `barcodes.tsv`.")
            return

        try:
            adata = load_file_dynamic(file_path, file_type, extra if file_type == "mtx" else None)
            adata = load_and_preprocess(adata)
            adata = apply_pca_umap_clustering(adata)

            st.session_state['adata'] = adata
            st.session_state.setdefault('all_datasets', []).append({
                "label": uploaded_file.name,
                "adata": adata
            })

            st.success(f"{file_type.upper()} file loaded and processed successfully!")

        except Exception as e:
            st.error(f"Error: {str(e)}")

    col1, col2 = st.columns(2)

    with col1:
        if st.button("Use PBMC 3k (Scanpy)"):
            try:
                adata = sc.datasets.pbmc3k()
                adata = load_and_preprocess(adata)
                adata = apply_pca_umap_clustering(adata)
                st.session_state['adata'] = adata
                st.session_state.setdefault('all_datasets', []).append({
                    "label": "PBMC_3k (Scanpy)",
                    "adata": adata
                })
                st.success("PBMC 3k loaded! Now Move to the Visualize tab.")
            except Exception as e:
                st.error(str(e))

    with col2:
        if st.button("Use Paul15 (Scanpy)"):
            try:
                adata = sc.datasets.paul15()
                adata = load_and_preprocess(adata)
                adata = apply_pca_umap_clustering(adata)
                st.session_state['adata'] = adata
                st.session_state.setdefault('all_datasets', []).append({
                    "label": "Paul15",
                    "adata": adata
                })
                st.success("Paul15 loaded! Now Move to the Visualize tab.")
            except Exception as e:
                st.error(str(e))

def render_visualizationTab():
    render_visualization()



def render_ai_assistant():
    st.title("AI Assistant: Ask about RNA-seq")

    with st.chat_message("assistant"):
        st.markdown("Hi! Ask anything about RNA-seq or this tool.")

    question = st.chat_input("Type your question here...")
    if question:
        with st.chat_message("user"):
            st.markdown(question)

        with st.chat_message("assistant"):
            if "gene" in question.lower():
                st.markdown("A gene is a segment of DNA that encodes a protein.")
            elif "umap" in question.lower():
                st.markdown("UMAP is a dimensionality reduction technique to visualize clusters.")
            elif "h5ad" in question.lower():
                st.markdown("`.h5ad` is a format used by Scanpy to store annotated data matrices.")
            else:
                st.markdown("I'm a basic assistant. Try asking about genes, UMAP, or formats.")
