# 🧬 GeneGazer

**GeneGazer** is an interactive web-based tool for visualizing and exploring single-cell RNA-seq datasets. It supports a wide range of input formats, provides cluster-based UMAP visualization, gene expression heatmaps, marker gene discovery, and integrates a powerful AI assistant powered by Gemini for bioinformatics guidance.

---

## 🚀 Features

- Upload and process `.h5ad`, `.csv`, or `.mtx` (with `genes.tsv` + `barcodes.tsv`) formats
- Preprocessing: normalization, log-transformation, HVG selection, scaling
- Interactive UMAP plot with cluster annotations and gene-level hover
- Cluster-wise gene heatmap with expression summary
- Automatically suggested top marker genes per cluster
- AI assistant powered by Gemini 1.5 Flash to answer RNA-seq related questions

---

## 🗂️ Project Structure

```
GeneGazer/
├── app.py
├── ui_utils.py
├── visualization_ui.py
├── visualizer.py
├── sc_processing.py
├── file_handler.py
├── ai_assistant.py
└── .streamlit/
    └── secrets.toml
```

---

## 📁 Module Overview

### `app.py`

Main entry point of the app. Uses `st.tabs()` to display:
- Instructions
- Data uploader
- Visualization panel
- Gemini AI Assistant

Sets up the project title, layout, and overall routing.

### `ui_utils.py`

Contains UI rendering functions for:
- Instructions page
- Load data tab

It also imports `render_visualization()` from `visualization_ui.py` to keep modularity.

### `visualization_ui.py`

Handles all visual interaction with processed data:
- `render_dataset_selector`: Lets user switch between uploaded datasets
- `render_umap_section`: Interactive UMAP plot with cluster labels and optional gene hover
- `render_gene_heatmap_section`: Select and visualize gene expression by cluster
- `render_marker_gene_table`: Automatically suggests top marker genes using differential expression
- `render_visualization`: Orchestrates all the above sections in the Visualize tab

### `visualizer.py`

Contains plotting logic using:
- `plot_umap`: UMAP with color-coded clusters using Plotly
- `plot_gene_heatmap`: Matrix plot to show gene expression across clusters
- `get_best_label_column`: Smartly chooses the best available column to label clusters (`cell_type`, `leiden`, etc.)

### `sc_processing.py`

Handles Scanpy-based data processing:
- `load_and_preprocess`: Filtering, normalization, HVG selection, scaling
- `apply_pca_umap_clustering`: Dimensionality reduction and clustering
- `compute_top_marker_genes`: Rank genes per cluster using t-test and return top hits

### `file_handler.py`

Handles file loading and validation:
- `detect_file_type`: Infers the format based on extension
- `load_file_dynamic`: Dynamically loads `.h5ad`, `.csv`, or `.mtx` with metadata
- `get_csv_download_link`: Provides a sample CSV for users to download and format their data

### `ai_assistant.py`

Contains the `GeminiAssistant` class, which:
- Sets up Gemini 1.5 Flash using API key
- Generates context-aware responses based on uploaded dataset
- Maintains chat history and allows user to ask any RNA-seq or bioinformatics question
- Handles errors, clears chat history, and displays assistant output in `st.chat_message()` UI

---

## 📥 Download Testing Data

- **MTX File Type and Supported Files:** [Download](https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-10371/download/zip?fileType=normalised&accessKey=)
- **CSV File Type:** [Download](https://storage.googleapis.com/kaggle-data-sets/1584326/2606779/compressed/brain_counts.csv.zip?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=gcp-kaggle-com%40kaggle-161607.iam.gserviceaccount.com%2F20250614%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20250614T080124Z&X-Goog-Expires=259200&X-Goog-SignedHeaders=host&X-Goog-Signature=de63917d3fe6236b05a2c7d194fed8be2136daeec7cb2e557bf6804cfb9a2f5388cc13d1f6fc4088c28a1bcff274018380bd23d838c60860494d90d4ffe921d2531d228a983e54d50ba1c675f8186763c9678a8164d6b11504ed089d41d2e4a586744c6fd8edb7af5d0974800d00cbde1909d1a37ab7ebb7cae54a6c64089dc8ff1b07c515a87514e9a7bf7c962c968fb604262d7a4138d88188df8a81463eea78a6eb1fb3d6c77bdf3ea6071e47d2e80e74e4e1fd3269a084d19ad684b608efe6fc7f534bcb30d8d917dde66eeb0424987b81187a9a6273a1a8fab088fb08714a876dfb617db6527da90cdef7210cdb4ef6346a15e0085c4605c1ca620e853b)
- **H5ad File Type:** [Download](https://storage.googleapis.com/kaggle-data-sets/1584326/2606779/compressed/glioblastoma_normalized.h5ad.zip?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=gcp-kaggle-com%40kaggle-161607.iam.gserviceaccount.com%2F20250614%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20250614T062652Z&X-Goog-Expires=259200&X-Goog-SignedHeaders=host&X-Goog-Signature=8ce6b9f4169cd46ab3e9bd0a9bd6bb846cac6a65fbb3f3d27acb1f17c0ffcd2a2e7f6449776ce2b1f73c0a50f63911197bebd291d6d2030a53d876f21aa41949ffef915c97044cc30d6f42cfcbf2d916401a6425c94ee0f42e77f1c60949449198a6c711d14fea6df317c10f283a0d152bcb0072c3a72205ba47eadb6424c6ad4a1cf1ef663827e0030f6c7c3fb0f1185ee1c348d7bb8be297ca756809a0426bbc5571f74e5f2d9047e88bf916f5fb0095102fe47d47c3b953dd347d968d4d3b502a51d5c3494ef32ce5dd43825f72fdfd0d2b896627d5e3d7378a8aefbf6eea07f967de872686da0610d03b399375db157ac3fb295d21139d6147060f6a2e58)

---

## 🔐 API Key Setup

To enable the AI Assistant, you must provide your Gemini API key.

**Option 1: Add to `.streamlit/secrets.toml`**

```toml
GEMINI_API_KEY = "your_api_key_here"
```

**Option 2: Set as environment variable**

```bash
export GEMINI_API_KEY=your_api_key_here
```

---

## ✅ Supported Input Formats

- `.h5ad` (AnnData)
- `.csv` (Genes × Cells or Cells × Genes with transposing)
- `.mtx` (Requires `genes.tsv` and `barcodes.tsv`)

---

## 📦 Requirements

Install all dependencies using:

```bash
pip install -r requirements.txt
```

---

## 💻 Run the App

```bash
streamlit run app.py
```

---

## 🧪 Demo Use Case

Try loading a public dataset like PBMC 3k or Paul15 to:
- Visualize cell clusters
- Identify top marker genes
- Interact with the AI assistant to interpret gene functions

---

## 🙌 Credits

Built for the [Your Hackathon Name]  
Team: GeneGazer  
Powered by Streamlit, Scanpy, and Google Gemini