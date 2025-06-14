import scanpy as sc
import streamlit as st
import matplotlib.pyplot as plt

def plot_umap(adata):
    sc.pl.umap(adata, color='leiden', show=False)
    st.pyplot(plt.gcf())

def plot_gene_heatmap(adata, gene):
    if gene in adata.var_names:
        sc.pl.matrixplot(adata, var_names=[gene], groupby="leiden", cmap="viridis", use_raw=False, show=False)
        st.pyplot(plt.gcf())
    else:
        st.warning("Gene not found in dataset.")
