Explore Dataset Contents
========================

The Explore Dataset Contents module allows you to inspect the internal structure of your single-cell RNA-seq dataset.

Overview
--------

This module displays comprehensive information about your AnnData object, including:

- Dataset dimensions (cells × genes)
- Available metadata fields
- Expression matrix format (sparse/dense)
- Dimensionality reductions
- Interactive data tables

Accessing the Module
--------------------

1. Load a dataset from the home screen
2. Navigate to "Explore Dataset Contents" in the sidebar

Dataset Summary
---------------

The module displays key information about your dataset:

**Dimensions**
   Number of cells (observations) and genes (variables) in the dataset

**Expression Matrix**
   Data type and sparsity of the expression matrix (adata.X)

**Available Reductions**
   List of computed dimensionality reductions (UMAP, t-SNE, PCA)

Cell Metadata (obs)
-------------------

The **Cell Metadata** tab shows the ``adata.obs`` DataFrame containing cell-level annotations:

- Sample identifiers
- Cell types
- Experimental conditions
- Quality control metrics
- Cluster assignments

**Features:**

- Interactive table with search and filtering
- Column statistics (unique values, data types)
- Export to CSV

Example metadata columns::

   sample_id        Biological sample identifier
   cell_type        Annotated cell type
   condition        Experimental condition (e.g., control, treatment)
   n_genes          Number of genes detected
   total_counts     Total UMI counts

Gene Metadata (var)
-------------------

The **Gene Metadata** tab displays the ``adata.var`` DataFrame with gene-level information:

- Gene symbols
- Gene IDs (Ensembl, Entrez)
- Feature selection flags
- Gene biotype
- QC metrics

**Features:**

- Search for specific genes
- Sort by columns
- View gene annotations

Example metadata columns::

   gene_ids         Ensembl gene identifier
   feature_types    Gene biotype (protein_coding, lncRNA)
   highly_variable  Boolean flag for variable genes

Expression Matrix (X)
---------------------

The **Expression Matrix** section provides information about ``adata.X``:

- Matrix dimensions (cells × genes)
- Data type (float32, int32, etc.)
- Storage format (sparse CSR, dense array)
- Memory usage

**Sparse vs Dense:**

scVIZ automatically handles both sparse and dense matrices. Sparse matrices are common in scRNA-seq to save memory.

Dimensionality Reductions
--------------------------

The **Available Reductions** section lists computed embeddings in ``adata.obsm``:

- ``X_umap``: UMAP coordinates
- ``X_tsne``: t-SNE coordinates  
- ``X_pca``: PCA coordinates

These reductions are used by the Visualize Dataset module for plotting.

Tips
----

✅ **Search Functionality**
   Use the search box in data tables to quickly find specific cells or genes

✅ **Export Data**
   Download metadata tables as CSV for further analysis

✅ **Check Before Analysis**
   Verify that required metadata columns exist before running differential expression

Common Use Cases
----------------

Verify Sample Identifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~

Before differential expression analysis:

1. Go to Cell Metadata tab
2. Check the sample identifier column (e.g., "sample_id")
3. Verify sufficient replicates per condition

Check Gene Annotations
~~~~~~~~~~~~~~~~~~~~~~

To find genes of interest:

1. Go to Gene Metadata tab
2. Search for gene symbols or IDs
3. Note the index for expression visualization

Inspect Data Quality
~~~~~~~~~~~~~~~~~~~~

Review QC metrics:

1. Check ``n_genes`` and ``total_counts`` in Cell Metadata
2. Look for outliers or low-quality cells
3. Verify preprocessing steps were applied

Next Steps
----------

- Create visualizations with :doc:`visualize_module`
- Run differential expression with :doc:`de_module`
