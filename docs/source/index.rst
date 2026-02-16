scVIZ Documentation
===================

**scVIZ** is an interactive web application for exploring and analyzing single-cell RNA-sequencing datasets using Streamlit and Scanpy.

.. image:: https://img.shields.io/badge/python-3.9%2B-blue
   :alt: Python Version

.. image:: https://img.shields.io/badge/streamlit-1.30%2B-red
   :alt: Streamlit

Features
--------

� **Dataset Loading**
   Load datasets from local h5ad/h5 files, curated example datasets from publications,
   or browse and download directly from CELLxGENE Discover portal.

📊 **Explore Dataset Contents**
   Examine the internal structure of AnnData objects including cell metadata (obs), 
   gene metadata (var), expression matrices (X), and dimensionality reductions (obsm/varm).

🔍 **Quality Control**
   Compute QC metrics (total counts, genes detected, mitochondrial %), perform cell 
   cycle scoring, and visualize distributions with interactive UMAP and violin plots.

🗺️ **Visualize Dataset**
   Create interactive UMAP/t-SNE visualizations with metadata coloring and gene expression 
   overlays. Support for lasso selection and customizable styling.

📈 **Differential Expression**
   Perform pseudobulk differential expression analysis using PyDESeq2 with 
   interactive volcano plots and downloadable results.

Quick Start
-----------

Installation
~~~~~~~~~~~~

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/jdwidawski/scVIZ.git
   cd scVIZ

   # Create virtual environment
   uv venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate

   # Install dependencies
   uv pip install -e ".[dev]"

Running the Application
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   streamlit run app.py

The application will open in your default web browser at http://localhost:8501.

Loading Data
~~~~~~~~~~~~

scVIZ supports multiple data loading methods:

1. **Upload File**: Upload h5ad or 10x h5 files from your computer
2. **Example Datasets**: Load from curated example datasets from publications
3. **CELLxGENE Datasets**: Browse and download single-cell/nuclei datasets from CELLxGENE Discover

All datasets must be in AnnData h5ad format or 10x Genomics h5 format.

User Guide
----------

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   user_guide/getting_started
   user_guide/explore_module
   user_guide/qc_module
   user_guide/visualize_module
   user_guide/de_module

API Reference
-------------

.. toctree::
   :maxdepth: 2
   :caption: API Documentation:

   api/modules
   api/utils

Configuration
-------------

Theme Customization
~~~~~~~~~~~~~~~~~~~

scVIZ uses a custom theme defined in ``.streamlit/config.toml``:

.. code-block:: toml

   [theme]
   base = "light"
   primaryColor = "#3b82f6"
   backgroundColor = "#ffffff"
   secondaryBackgroundColor = "#f8f9fa"
   textColor = "#1a202c"

You can customize colors and styling by editing this configuration file.

Contributing
------------

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create a feature branch (``git checkout -b feature/amazing-feature``)
3. Commit changes (``git commit -m 'feat: add amazing feature'``)
4. Push to branch (``git push origin feature/amazing-feature``)
5. Open a Pull Request

License
-------

This project is licensed under the MIT License - see the LICENSE file for details.

Copyright (c) 2026 Jakub Widawski

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
