import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile
from file_handler import detect_file_type, load_file_dynamic


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
menu = st.sidebar.radio("Select a section:", ["📘 Instructions", "📁 Load Data", "📊 Visualize", "🗂️ All Datasets", "🤖 AI Assistant"])


if menu == "📘 Instructions":
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

    st.image("https://i.imgur.com/zVfGZkP.png", caption="Example of UMAP visualization")

elif menu == "📁 Load Data":
    st.title("Load RNA-seq Data")

    st.markdown("### Upload your data file")

    uploaded_file = st.file_uploader("Upload a `.h5ad`, `.csv`, or `.mtx` file", type=["h5ad", "csv", "mtx"])

    extra_genes = st.file_uploader("Upload genes.tsv (for .mtx)", type=["tsv"]) if uploaded_file and uploaded_file.name.endswith("mtx") else None
    extra_barcodes = st.file_uploader("Upload barcodes.tsv (for .mtx)", type=["tsv"]) if uploaded_file and uploaded_file.name.endswith("mtx") else None

    if uploaded_file:
        from file_handler import detect_file_type, load_file_dynamic
        import tempfile

        file_type = detect_file_type(uploaded_file.name)

        with tempfile.NamedTemporaryFile(delete=False, suffix=f".{file_type}") as tmp:
            tmp.write(uploaded_file.read())
            file_path = tmp.name

        extra = {"genes": extra_genes, "barcodes": extra_barcodes} if file_type == "mtx" else None

        try:
            adata = load_file_dynamic(file_path, file_type, extra)
            adata = load_and_preprocess(adata)
            adata = apply_pca_umap_clustering(adata)
            # Save current dataset as active
            st.session_state['adata'] = adata

            if 'all_datasets' not in st.session_state:
             st.session_state['all_datasets'] = []

            dataset_label = uploaded_file.name if uploaded_file else "PBMC_3k" if "pbmc3k" in file_path else "Mouse_Brain"
            st.session_state['all_datasets'].append({
               "label": dataset_label,
              "adata": adata
            })

            st.success(f"{file_type.upper()} file loaded and processed successfully!")
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

elif menu == "🗂️ All Datasets":
    st.title("Previously Uploaded Datasets")

    if 'all_datasets' not in st.session_state or len(st.session_state['all_datasets']) == 0:
        st.warning("No datasets uploaded yet.")
    else:
        dataset_names = [ds['label'] for ds in st.session_state['all_datasets']]
        selected_label = st.selectbox("Choose a dataset to view:", dataset_names)

        selected_dataset = next(ds for ds in st.session_state['all_datasets'] if ds['label'] == selected_label)
        adata = selected_dataset['adata']

        st.subheader(f"UMAP for: {selected_label}")
        sc.pl.umap(adata, color='leiden', show=False)
        st.pyplot(plt.gcf())

        gene = st.text_input("Enter a gene name to plot:", "IL7R")
        if gene in adata.var_names:
            sc.pl.matrixplot(adata, var_names=[gene], groupby="leiden", cmap="viridis", use_raw=False, show=False)
            st.pyplot(plt.gcf())
        else:
            st.warning("Gene not found in this dataset.")



st.markdown("---")
st.info("Let's shine in this hackathon, team Goat!")

