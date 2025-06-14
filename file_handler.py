# file_handler.py

import scanpy as sc
import pandas as pd
import os
import tempfile
from scipy import io
import anndata as ad

def detect_file_type(file_name):
    if file_name.endswith(".h5ad"):
        return "h5ad"
    elif file_name.endswith(".csv"):
        return "csv"
    elif file_name.endswith(".mtx"):
        return "mtx"
    else:
        return "unknown"

def load_file_dynamic(file, file_type, extra_files=None):
    if file_type == "h5ad":
        return sc.read_h5ad(file)

    elif file_type == "csv":
        df = pd.read_csv(file, index_col=0).T  # transpose to cells x genes
        adata = sc.AnnData(df)
        return adata

    elif file_type == "mtx" and extra_files:
        matrix = io.mmread(file).T.tocsr()
        genes = pd.read_csv(extra_files['genes'], header=None, sep="\t")
        barcodes = pd.read_csv(extra_files['barcodes'], header=None, sep="\t")
        adata = ad.AnnData(X=matrix)
        adata.var_names = genes[1]
        adata.obs_names = barcodes[0]
        return adata

    else:
        raise ValueError("Unsupported or incomplete file format.")
