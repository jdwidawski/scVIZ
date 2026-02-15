# 🧬 scVIZ

**Interactive Single-Cell RNA-Sequencing Data Visualization and Analysis Tool**

[![Tests](https://github.com/jakub-widawski/scviz/actions/workflows/tests.yml/badge.svg)](https://github.com/jakub-widawski/scviz/actions/workflows/tests.yml)
[![Documentation](https://github.com/jakub-widawski/scviz/actions/workflows/docs.yml/badge.svg)](https://jakub-widawski.github.io/scviz/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

scVIZ is a Streamlit-based interactive application for exploring and analyzing single-cell RNA-sequencing datasets using Scanpy. It provides an intuitive interface for researchers to visualize, quality check, and perform differential expression analysis on their scRNA-seq data.

## ✨ Features

### 📥 Dataset Loading
- **Upload local files**: Load h5ad format AnnData objects
- **Example datasets**: Curated collection of human scRNA-seq datasets with UMAP embeddings
- **CELLxGENE integration**: Browse and load datasets directly from CELLxGENE Discover

### 📊 Explore Dataset Contents
- View dataset structure and metadata
- Explore observation (cell) and variable (gene) annotations
- Inspect dimensionality reductions (PCA, UMAP, t-SNE)
- Browse unstructured data (clustering results, marker genes)

### 🔍 Quality Control
- Compute QC metrics (total counts, genes detected, mitochondrial %)
- Cell cycle phase scoring (S, G2M phases)
- Interactive UMAP visualization with lasso selection
- Violin plots for metric distribution analysis
- Support for both gene symbols and Ensembl IDs

### 🗺️ Visualize Dataset
- Interactive UMAP/t-SNE plots with Plotly
- Color by categorical metadata or gene expression
- Customizable figure settings (colormap, dot size, color bounds)
- Lasso selection for cell subset analysis
- Export-ready visualizations

### 📈 Differential Expression Analysis
- Pseudobulk differential expression with PyDESeq2
- Interactive group selection via lasso on UMAP
- Volcano plots and results tables
- LFC shrinkage options
- Scientific notation for p-values

## 🚀 Quick Start

### Prerequisites
- Python 3.9 or higher
- [UV](https://docs.astral.sh/uv/) (recommended) or pip

### Installation

```bash
# Clone the repository
git clone https://github.com/jakub-widawski/scviz.git
cd scviz

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

## 📖 Documentation

Full documentation is available at [jakub-widawski.github.io/scviz/](https://jakub-widawski.github.io/scviz/)

## 🧪 Running Tests

```bash
pytest
```

For coverage report:
```bash
pytest --cov=utils --cov=modules --cov-report=html
```

## 📁 Project Structure

```
scviz/
├── app.py                 # Main Streamlit application
├── modules/               # Analysis modules
│   ├── explore_dataset.py     # Dataset exploration
│   ├── quality_control.py     # QC metrics and visualization
│   ├── visualize_dataset.py   # UMAP/expression visualization
│   ├── visualize_expression.py # Gene expression plots
│   └── differential_expression.py # DE analysis
├── utils/                 # Utility functions
│   └── data_loader.py         # Dataset I/O functions
├── tests/                 # Test suite
│   └── test_data_loader.py
├── docs/                  # Sphinx documentation
├── pyproject.toml         # Project configuration
└── README.md
```

## 🔧 Development

### Module Development Pattern

New analysis modules should:
1. Be placed in `modules/` directory as separate `.py` files
2. Implement a `render(adata: AnnData) -> None` function
3. Import and register in `app.py`'s module routing logic
4. Include a help modal with usage instructions

```python
def render(adata: AnnData) -> None:
    """Render the module UI."""
    st.title("Module Name")
    # ... module logic
```

### Code Standards
- Follow PEP 8 style guidelines
- Use type hints with `typing` module
- Write docstrings with Parameters/Returns sections
- Test utilities in `tests/` directory

## 📋 Dependencies

- **Streamlit** (>=1.30.0): Web application framework
- **Scanpy** (>=1.9.0): Single-cell analysis toolkit
- **Plotly** (>=5.18.0): Interactive visualizations
- **PyDESeq2** (>=0.4.0): Differential expression analysis
- **Pandas** (>=2.0.0): Data manipulation
- **NumPy** (>=1.24.0): Numerical computing

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👤 Author

**Jakub Widawski**

- GitHub: [@jakub-widawski](https://github.com/jakub-widawski)

## 🙏 Acknowledgments

- [Scanpy](https://scanpy.readthedocs.io/) - Single-cell analysis in Python
- [Streamlit](https://streamlit.io/) - The fastest way to build data apps
- [CELLxGENE](https://cellxgene.cziscience.com/) - Single-cell data portal
- [PyDESeq2](https://pydeseq2.readthedocs.io/) - Python implementation of DESeq2
