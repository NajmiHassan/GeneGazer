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
- **CSV File Type:** [Download](https://storage.googleapis.com/kaggle-data-sets/1584326/2606779/compressed/brain_counts.csv.zip)
- **H5ad File Type:** [Download](https://storage.googleapis.com/kaggle-data-sets/1584326/2606779/compressed/glioblastoma_normalized.h5ad.zip)

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