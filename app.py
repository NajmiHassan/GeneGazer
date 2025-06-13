import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import sklearn
from sc_processing import load_and_preprocess, run_pca_umap_cluster

st.set_page_config(page_title="Single-Cell RNA-Seq Visualization Tool", layout="wide")
st.title("Single-Cell RNA-Seq Visualization Tool")

st.sidebar.header("Upload & Settings")

uploaded_file = st.sidebar.file_uploader("Upload your .h5ad file", type=["h5ad"])

if uploaded_file is not None:
    st.success("File uploaded successfully!")
    adata = load_and_preprocess(uploaded_file)
    adata = run_pca_umap_cluster(adata)
    st.subheader("UMAP Clustering Plot")
    fig, ax = plt.subplots(figsize=(7,5))
    sc.pl.umap(adata, color=["leiden"], ax=ax, show=False)
    st.pyplot(fig)
else:
    st.info("Please upload a .h5ad file to begin.")
