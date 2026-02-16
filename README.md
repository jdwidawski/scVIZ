# scVIZ

**Interactive Single-Cell RNA-Sequencing Data Visualization and Analysis Tool**

[![Tests](https://github.com/jdwidawski/scVIZ/actions/workflows/tests.yml/badge.svg)](https://github.com/jdwidawski/scVIZ/actions/workflows/tests.yml)
[![Documentation](https://github.com/jdwidawski/scVIZ/actions/workflows/docs.yml/badge.svg)](https://jdwidawski.github.io/scVIZ/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

scVIZ is a Streamlit-based interactive application for exploring and analyzing single-cell RNA-sequencing datasets using Scanpy. It provides an intuitive interface for researchers to visualize, quality check, and perform differential expression analysis on their scRNA-seq data.

## Features

### Dataset Loading
Load datasets from local files (h5ad/h5 formats), curated examples, or directly from CELLxGENE Discover.

### Explore Dataset
Inspect dataset structure, metadata, dimensionality reductions, and gene/cell annotations.

### Quality Control
Compute QC metrics, cell cycle scoring, interactive exploration with violin plots and UMAP visualization.

### Visualize Dataset
Interactive UMAP plots colored by metadata or gene expression with customizable settings.

### Differential Expression
Pseudobulk DE analysis with PyDESeq2, interactive group selection, volcano plots, and results export.

## Quick Start

### Prerequisites
- Python 3.9 or higher
- [UV](https://docs.astral.sh/uv/) (recommended) or pip

### Installation

```bash
# Clone the repository
git clone https://github.com/jdwidawski/scVIZ.git
cd scVIZ

# Create virtual environment and install dependencies (using UV)
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e ".[dev]"

# Or using pip
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

### Running the App

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

## Documentation

Full documentation is available at [jdwidawski.github.io/scVIZ/](https://jdwidawski.github.io/scVIZ/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
