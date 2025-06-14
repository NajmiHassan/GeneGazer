
# 🧬 Single-Cell RNA-seq Visualization Tool

This is a user-friendly web app built with **Streamlit** and **Scanpy** for exploring and visualizing single-cell RNA sequencing data. It supports multiple data formats, generates clustering and gene expression visualizations, and includes an AI assistant for beginners.

---

## 🚀 Features

- ✅ Upload support for:
  - `.h5ad` (AnnData format from Scanpy)
  - `.csv` gene expression matrix (cells in columns or rows)
  - `.mtx` + `genes.tsv` + `barcodes.tsv` (10x Genomics format)
- 📈 Visualize UMAP plots for clusters
- 🎯 Gene expression heatmaps by cluster
- 🧠 Basic AI assistant for RNA-seq help
- 🗂️ View and revisit previously uploaded datasets
- 📥 Downloadable `.csv` template
- 🔒 Ensures only valid file formats are uploaded

---

## 📁 Folder Structure

```
.
├── app.py                # Main Streamlit app
├── file_handler.py       # Detect and load different file formats
├── sc_processing.py      # Preprocessing and clustering pipeline
├── README.md             # This file
├── requirements.txt      # Required Python packages
└── TestingData/          # Testing files to test different file types. 
└── Date/                 # Public Datasets for testing.  
```

---

## 📦 Installation

### 1. Clone this repository

```bash
git clone https://github.com/NajmiHassan/GeneGazer.git
cd GeneGazer
```

### 2. (Optional but recommended) Create a virtual environment

```bash
python -m venv venv
source venv/bin/activate    # On Windows: venv\Scripts\activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

If `requirements.txt` is not available, install manually:

```bash
pip install streamlit scanpy matplotlib pandas seaborn anndata scipy scikit-learn
```

---

## 🧪 Running the App

```bash
streamlit run app.py
```

It will open the application in your browser at `http://localhost:8501`.

---

## 📤 File Upload Instructions

### 🔹 Supported Formats

| Format         | Description |
|----------------|-------------|
| `.h5ad`        | Scanpy AnnData object |
| `.csv`         | Gene x Cell matrix (or transposed) |
| `.mtx` bundle  | Must upload `.mtx`, `genes.tsv` or `genes.mtx_rows`, and `barcodes.tsv` or `barcodes.mtx_cols` |

> ⚠️ For `.mtx` upload, **you must upload all three files together** from 10x Genomics outputs.

---

## 📥 Downloadable Template

You can download a test `.csv` file directly from the app (under **Load Data** section) to fill your data and then upload. Don't try to upload template because it doesn't have enough data. 

Example structure:

|        | Cell1 | Cell2 | Cell3 |
|--------|-------|-------|-------|
| GeneA  | 0     | 1     | 2     |
| GeneB  | 2     | 0     | 4     |
| GeneC  | 3     | 5     | 1     |

---

## 📊 Available Visualizations

- **UMAP Clustering Plot** (Leiden clustering)
- **Gene Heatmap** (Matrix plot grouped by clusters)
- **Multiple Datasets Tab** to switch between uploaded datasets

---

## 🧠 AI Assistant

A built-in chatbot can answer basic RNA-seq questions like:

- “What is UMAP?”
- “What is `.h5ad`?”
- “What is gene expression?”

---

## 🛠️ Troubleshooting

| Problem | Solution |
|--------|----------|
| `.csv` error | Ensure numeric values, and transpose if needed |
| `.mtx` error | Upload `mtx`, `genes.tsv`, and `barcodes.tsv` together |
| File not uploading | Check file extensions and sizes |

---

## 📚 Public Dataset Sources

- [ArrayExpress (EMBL-EBI)](https://www.ebi.ac.uk/biostudies/arrayexpress/)
- [GEO Datasets (NCBI)](https://www.ncbi.nlm.nih.gov/geo/)
- [10x Genomics](https://www.10xgenomics.com/resources/datasets)

---

## 📄 License

This project is licensed under the **MIT License**. You are free to use, modify, and distribute it.

---

## 🙌 Acknowledgements

- [Scanpy](https://scanpy.readthedocs.io/en/stable/)
- [Streamlit](https://streamlit.io/)
- [10x Genomics](https://www.10xgenomics.com/)
- Team Goat for contributions 🚀

---

## ✉️ Contact

For questions or support, please open an issue or contact the developer.
