import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import tempfile
import os
import requests

st.set_page_config(page_title="🔬 RNA-seq Viewer", layout="wide")

# ======= SIDEBAR ========
st.sidebar.title("🔍 Navigation")
menu = st.sidebar.radio("Choose a section:", ["📘 Instructions", "📁 Load Data", "📊 Visualize", "🤖 AI Agent"])
st.sidebar.markdown("---")
st.sidebar.info("Developed at TraeIA Hackathon 🧠")

# ======= INSTRUCTIONS ========
if menu == "📘 Instructions":
    st.title("🧬 Single-Cell RNA-seq Viewer")
    st.markdown("""
    This tool allows you to easily explore single-cell RNA-seq data:

    1. Go to **\"Load Data\"** and upload your `.h5ad` file or use our example dataset.
    2. Go to **\"Visualize\"** to generate interactive plots:
        - Clustering plot (UMAP)
        - Gene expression heatmap
    3. Use the **\"AI Agent\"** tab to ask questions.

    #### Accepted format:
    - `.h5ad` files exported from **Scanpy** or **Seurat**
    """)
    st.image("https://i.imgur.com/zVfGZkP.png", caption="Example of UMAP visualization")

# ======= LOAD DATA ========
elif menu == "📁 Load Data":
    st.title("📤 Load RNA-seq Data")

    uploaded_file = st.file_uploader("Upload a `.h5ad` file", type=["h5ad"])
    if uploaded_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            tmp.write(uploaded_file.read())
            temp_path = tmp.name
        st.session_state['adata'] = sc.read_h5ad(temp_path)
        st.success("File loaded successfully!")

    st.markdown("### Or choose one of the public datasets:")
    col1, col2 = st.columns(2)

    with col1:
        if st.button("🔬 Use PBMC 3k"):
            url = "https://figshare.com/ndownloader/files/28125275"
            try:
                with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
                    response = requests.get(url)
                    tmp.write(response.content)
                    temp_path = tmp.name
                st.session_state['adata'] = sc.read_h5ad(temp_path)
                st.success("PBMC 3k loaded successfully!")
            except Exception as e:
                st.error(f"Error: {str(e)}")

    with col2:
        if st.button("🧠 Use Mouse Brain (subsampled)"):
            url = "https://github.com/scverse/scvi-tools/raw/main/tests/data/mouse_brain_subsampled.h5ad"
            try:
                with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
                    response = requests.get(url)
                    tmp.write(response.content)
                    temp_path = tmp.name
                st.session_state['adata'] = sc.read_h5ad(temp_path)
                st.success("Mouse Brain (subsampled) loaded successfully!")
            except Exception as e:
                st.error(f"Error: {str(e)}")

# ======= VISUALIZE ========
elif menu == "📊 Visualize":
    st.title("📈 Visualizations")

    if 'adata' not in st.session_state:
        st.warning("⚠️ Please upload a file in the 'Load Data' tab first.")
    else:
        adata = st.session_state['adata']

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("🔵 UMAP")
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata)
            sc.pl.umap(adata, color='leiden', show=False)
            st.pyplot(plt.gcf())

        with col2:
            st.subheader("🔥 Heatmap")
            gene = st.text_input("Enter a gene name to visualize:", "IL7R")
            if gene in adata.var_names:
                sc.pl.matrixplot(adata, var_names=[gene], groupby="leiden", cmap="viridis", use_raw=False, show=False)
                st.pyplot(plt.gcf())
            else:
                st.error("Gene not found.")

# ======= SIMPLE AI AGENT ========
elif menu == "🤖 AI Agent":
    st.title("🤖 AI Agent: Help and Explanations")

    with st.chat_message("assistant"):
        st.markdown("Hi! Ask me anything about RNA-seq, genes, or how to use this tool!")

    user_input = st.chat_input("Type your question here...")
    if user_input:
        with st.chat_message("user"):
            st.markdown(user_input)

        with st.chat_message("assistant"):
            if "gene" in user_input.lower():
                st.markdown("🧬 A gene is a DNA sequence that contains instructions for making proteins.")
            elif "umap" in user_input.lower():
                st.markdown("🔵 UMAP is a dimensionality reduction technique that helps visualize complex data in 2D.")
            elif "h5ad" in user_input.lower():
                st.markdown("📁 `.h5ad` file is a format used by Scanpy to store gene expression data.")
            else:
                st.markdown("Sorry, I'm still a simple agent 😅 — but I can help with basic questions!")

# Footer
st.markdown("---")
st.info("Just say the word. Let’s shine at this hackathon, Goat! 🚀🧠")
