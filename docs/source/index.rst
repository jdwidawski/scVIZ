scVIZ Documentation
===================

**scVIZ** is an interactive web application for exploring and analyzing single-cell RNA-sequencing datasets using Streamlit and Scanpy.

.. image:: https://img.shields.io/badge/python-3.9%2B-blue
   :alt: Python Version

.. image:: https://img.shields.io/badge/streamlit-1.30%2B-red
   :alt: Streamlit

Features
--------

📊 **Explore Dataset Contents**
   Examine the internal structure of AnnData objects including cell metadata (obs), 
   gene metadata (var), expression matrices (X), and dimensionality reductions (obsm/varm).

🗺️ **Visualize Dataset**
   Create interactive UMAP visualizations with metadata coloring and gene expression 
   overlays. Support for multi-gene comparison and customizable styling.

🧬 **Differential Expression**
   Perform pseudobulk differential expression analysis using PyDESeq2 with 
   interactive volcano plots and downloadable results.

Quick Start
-----------

Installation
~~~~~~~~~~~~

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/yourusername/scrnaseq-workflow.git
   cd scrnaseq-workflow

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

1. **Upload File**: Upload an h5ad file from your computer
2. **Example Datasets**: Load from curated example datasets with UMAP embeddings
3. **CELLxGENE Datasets**: Browse and download from CELLxGENE Discover portal

All datasets must be in AnnData h5ad format.

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
