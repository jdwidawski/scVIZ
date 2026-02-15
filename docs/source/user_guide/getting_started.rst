Getting Started
===============

This guide will help you get started with scVIZ.

Prerequisites
-------------

- Python 3.9 or higher
- UV package manager (recommended) or pip

Installation
------------

Using UV (Recommended)
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/yourusername/scrnaseq-workflow.git
   cd scrnaseq-workflow

   # Create and activate virtual environment
   uv venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate

   # Install the package in development mode
   uv pip install -e ".[dev]"

Using Pip
~~~~~~~~~

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/yourusername/scrnaseq-workflow.git
   cd scrnaseq-workflow

   # Create and activate virtual environment
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate

   # Install the package in development mode
   pip install -e ".[dev]"

Launching the Application
--------------------------

Once installed, launch scVIZ with:

.. code-block:: bash

   streamlit run app.py

The application will open in your default web browser at http://localhost:8501.

First Steps
-----------

1. **Load a Dataset**
   
   Start by loading a dataset using one of three methods:
   
   - **Upload File**: Click "Upload File" tab and select an h5ad file
   - **Example Datasets**: Choose from pre-configured example datasets
   - **CELLxGENE**: Browse and download datasets from CELLxGENE Discover

2. **Explore Your Data**
   
   After loading, navigate to "Explore Dataset Contents" to:
   
   - View metadata dimensions
   - Examine cell observations (obs)
   - Inspect gene variables (var)
   - Check available dimensionality reductions

3. **Visualize**
   
   Go to "Visualize Dataset" to create UMAP plots with:
   
   - Metadata coloring by categorical variables
   - Gene expression overlays
   - Customizable plot dimensions and colors

4. **Analyze**
   
   Use "Differential Expression" to:
   
   - Aggregate cells into pseudobulk samples
   - Run PyDESeq2 analysis
   - Explore results with volcano plots
   - Download results as CSV

Data Format Requirements
------------------------

scVIZ expects AnnData objects in h5ad format with:

- **adata.obs**: Cell-level metadata (sample IDs, cell types, conditions)
- **adata.var**: Gene-level metadata (gene names, symbols)
- **adata.X**: Expression matrix (cells × genes)
- **adata.obsm** (optional): Dimensionality reductions (X_umap, X_tsne, X_pca)

For differential expression analysis:

- A column identifying biological samples (e.g., "sample_id")
- A column specifying experimental conditions (e.g., "treatment", "disease")

Troubleshooting
---------------

Application Won't Start
~~~~~~~~~~~~~~~~~~~~~~~

**Error**: ``ModuleNotFoundError``

**Solution**: Ensure all dependencies are installed:

.. code-block:: bash

   uv pip install -e ".[dev]"

Dataset Loading Fails
~~~~~~~~~~~~~~~~~~~~~~

**Error**: ``File format not recognized``

**Solution**: Verify your file is in h5ad format. Convert using:

.. code-block:: python

   import scanpy as sc
   adata = sc.read_10x_h5("your_data.h5")  # or other format
   adata.write_h5ad("your_data.h5ad")

Missing UMAP Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Error**: No UMAP available

**Solution**: Compute UMAP before using visualization:

.. code-block:: python

   import scanpy as sc
   sc.pp.neighbors(adata)
   sc.tl.umap(adata)
   adata.write_h5ad("your_data.h5ad")

Next Steps
----------

- Learn about the :doc:`explore_module` to inspect your data
- Create visualizations with the :doc:`visualize_module`
- Perform analysis with the :doc:`de_module`
