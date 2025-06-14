
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

## 📥 Download Testing Data

- **MTX File Type and Supported Files:** [Download](https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-10371/download/zip?fileType=normalised&accessKey=)
- **CSV File Type:** [Download](https://storage.googleapis.com/kaggle-data-sets/1584326/2606779/compressed/brain_counts.csv.zip?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=gcp-kaggle-com%40kaggle-161607.iam.gserviceaccount.com%2F20250614%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20250614T080124Z&X-Goog-Expires=259200&X-Goog-SignedHeaders=host&X-Goog-Signature=de63917d3fe6236b05a2c7d194fed8be2136daeec7cb2e557bf6804cfb9a2f5388cc13d1f6fc4088c28a1bcff274018380bd23d838c60860494d90d4ffe921d2531d228a983e54d50ba1c675f8186763c9678a8164d6b11504ed089d41d2e4a586744c6fd8edb7af5d0974800d00cbde1909d1a37ab7ebb7cae54a6c64089dc8ff1b07c515a87514e9a7bf7c962c968fb604262d7a4138d88188df8a81463eea78a6eb1fb3d6c77bdf3ea6071e47d2e80e74e4e1fd3269a084d19ad684b608efe6fc7f534bcb30d8d917dde66eeb0424987b81187a9a6273a1a8fab088fb08714a876dfb617db6527da90cdef7210cdb4ef6346a15e0085c4605c1ca620e853b)
- **H5ad File Type:** [Download](https://storage.googleapis.com/kaggle-data-sets/1584326/2606779/compressed/glioblastoma_normalized.h5ad.zip?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=gcp-kaggle-com%40kaggle-161607.iam.gserviceaccount.com%2F20250614%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20250614T062652Z&X-Goog-Expires=259200&X-Goog-SignedHeaders=host&X-Goog-Signature=8ce6b9f4169cd46ab3e9bd0a9bd6bb846cac6a65fbb3f3d27acb1f17c0ffcd2a2e7f6449776ce2b1f73c0a50f63911197bebd291d6d2030a53d876f21aa41949ffef915c97044cc30d6f42cfcbf2d916401a6425c94ee0f42e77f1c60949449198a6c711d14fea6df317c10f283a0d152bcb0072c3a72205ba47eadb6424c6ad4a1cf1ef663827e0030f6c7c3fb0f1185ee1c348d7bb8be297ca756809a0426bbc5571f74e5f2d9047e88bf916f5fb0095102fe47d47c3b953dd347d968d4d3b502a51d5c3494ef32ce5dd43825f72fdfd0d2b896627d5e3d7378a8aefbf6eea07f967de872686da0610d03b399375db157ac3fb295d21139d6147060f6a2e58)

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
