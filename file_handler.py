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
    # --------------------------
    # .h5ad format (AnnData)
    # --------------------------
    if file_type == "h5ad":
        return sc.read_h5ad(file)

    # --------------------------
    # .csv format (expression matrix)
    # --------------------------
    elif file_type == "csv":
        df = pd.read_csv(file, index_col=0).T  # cells × genes
        try:
            df = df.apply(pd.to_numeric, errors="coerce")
            df = df.dropna(axis=1, how="any")
        except Exception as e:
            raise ValueError(f"Error parsing CSV file: {str(e)}")

        if df.shape[1] == 0:
            raise ValueError("CSV file has no valid numeric gene expression columns.")

        adata = sc.AnnData(df)
        return adata

    # --------------------------
    # .mtx + metadata files
    # --------------------------
    elif file_type == "mtx":
        if not extra_files or not extra_files.get("genes") or not extra_files.get("barcodes"):
            raise ValueError("Missing genes file or barcodes file for .mtx format.")

        try:
            # Save genes metadata (e.g., .tsv or .mtx_rows)
            with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as gfile:
                gfile.write(extra_files['genes'].read())
                genes_path = gfile.name

            # Save barcodes metadata (e.g., .tsv or .mtx_cols)
            with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as bfile:
                bfile.write(extra_files['barcodes'].read())
                barcodes_path = bfile.name

            # Load matrix
            matrix = io.mmread(file).T.tocsr()

            # Load genes (rows)
            try:
                genes = pd.read_csv(genes_path, header=None, sep="\t")
            except:
                genes = pd.read_csv(genes_path, header=None)

            # Load barcodes (columns)
            try:
                barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")
            except:
                barcodes = pd.read_csv(barcodes_path, header=None)

            # Validate
            if matrix.shape[0] != barcodes.shape[0]:
                raise ValueError(f"Matrix rows ({matrix.shape[0]}) != barcodes ({barcodes.shape[0]})")
            if matrix.shape[1] != genes.shape[0]:
                raise ValueError(f"Matrix columns ({matrix.shape[1]}) != genes ({genes.shape[0]})")

            adata = ad.AnnData(X=matrix)
            adata.obs_names = barcodes[0]
            adata.var_names = genes[1] if genes.shape[1] > 1 else genes[0]
            return adata

        except Exception as e:
            raise ValueError(f"Error loading .mtx format: {str(e)}")

    # --------------------------
    # Unsupported file
    # --------------------------
    else:
        raise ValueError("Unsupported or incomplete file format. Please upload a .h5ad, .csv, or .mtx with metadata.")
