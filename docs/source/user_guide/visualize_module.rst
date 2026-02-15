Visualize Dataset
=================

The Visualize Dataset module creates interactive UMAP visualizations for exploring cell metadata and gene expression patterns.

Overview
--------

This module provides two main visualization modes:

📊 **Metadata Visualization**
   Color cells by categorical metadata (cell types, conditions, clusters)

🧬 **Gene Expression Visualization**
   Overlay gene expression levels on UMAP plots

Prerequisites
-------------

Your dataset must have a computed UMAP reduction in ``adata.obsm['X_umap']``.

If UMAP is missing, compute it using Scanpy:

.. code-block:: python

   import scanpy as sc
   sc.pp.neighbors(adata)
   sc.tl.umap(adata)

Metadata Visualization
----------------------

Visualize categorical metadata fields on UMAP coordinates.

Selecting Metadata
~~~~~~~~~~~~~~~~~~

1. Navigate to the **Metadata Visualization** tab
2. Choose a categorical column from the dropdown (e.g., "cell_type", "sample_id")
3. The UMAP plot will update automatically

Plot Customization
~~~~~~~~~~~~~~~~~~

**Plot Dimensions**
   Adjust width and height using sliders (300-1200 pixels)

**Color Scheme**
   Select from 8 categorical color palettes:
   
   - Plotly (default)
   - D3
   - G10
   - T10
   - Alphabet
   - Dark24
   - Light24
   - Set3

**Legend**
   Always displayed showing category labels and colors

Category Statistics
~~~~~~~~~~~~~~~~~~~

The module displays summary statistics for each category:

- **Count**: Number of cells in the category
- **Percentage**: Proportion of total cells

These statistics help you understand the composition of your dataset.

Gene Expression Visualization
------------------------------

Overlay gene expression levels on UMAP plots with customizable color scales.

Selecting Genes
~~~~~~~~~~~~~~~

1. Navigate to the **Gene Expression Visualization** tab
2. Enter one or more gene names (comma-separated)
3. Each gene will be displayed in a separate tab

Example input::

   CD3D, CD8A, CD4

Finding Gene Names
~~~~~~~~~~~~~~~~~~

Use the Explore Dataset Contents module to:

1. Go to Gene Metadata tab
2. Search for genes of interest
3. Copy exact gene names from the index

Color Scale Options
~~~~~~~~~~~~~~~~~~~

**Color Mode**

- **Auto**: Automatically scales from 0.0 to 99.5th percentile
- **Percentile**: Custom percentile range (e.g., 25th-95th)
- **Absolute**: Fixed expression value range

**Colormaps**

Choose from 13 continuous colormaps:

- **Reds** (default): Best for sparse gene expression
- Viridis, Plasma, Inferno: Perceptually uniform
- Blues, Greens, Oranges, Purples: Single-hue scales
- Cividis, Turbo: Colorblind-friendly options
- Hot, Jet: Classic scientific colormaps

Plot Customization
~~~~~~~~~~~~~~~~~~

- **Width/Height**: Adjust plot size (300-1200 pixels)
- **Marker Size**: Control point size (1-10)

Expression Statistics
~~~~~~~~~~~~~~~~~~~~~

For each gene, view summary statistics:

- **Mean**: Average expression across all cells
- **Max**: Maximum expression value
- **% Expressing**: Percentage of cells with non-zero expression

Plot Ordering
~~~~~~~~~~~~~

Cells are automatically sorted by expression level, with higher-expressing cells plotted on top for better visibility.

Tips
----

✅ **Use Auto Mode**
   The Auto color scale (0.0 to 99.5th percentile) works well for most genes

✅ **Compare Multiple Genes**
   Enter comma-separated genes to compare expression patterns side-by-side

✅ **Choose Appropriate Colormaps**
   Use Reds for sparse data, Viridis/Plasma for continuous distributions

✅ **Check Expression Statistics**
   Review % Expressing to verify gene is detected in your dataset

Common Use Cases
----------------

Identify Cell Types
~~~~~~~~~~~~~~~~~~~~

1. Select "cell_type" in Metadata Visualization
2. Note the spatial distribution of cell populations
3. Compare with marker gene expression

Validate Marker Genes
~~~~~~~~~~~~~~~~~~~~~~

1. Enter known marker genes (e.g., "CD3D, CD8A" for T cells)
2. Check if expression matches expected cell populations
3. Use statistics to quantify marker specificity

Compare Conditions
~~~~~~~~~~~~~~~~~~

1. Select "condition" or "treatment" in Metadata Visualization
2. Look for spatial separation of groups
3. Overlay genes of interest to find differentially expressed markers

Explore Clusters
~~~~~~~~~~~~~~~~

1. Visualize clustering results (e.g., "leiden", "louvain")
2. Check cluster-specific marker genes
3. Annotate clusters based on expression patterns

Troubleshooting
---------------

Gene Not Found
~~~~~~~~~~~~~~

**Error**: Gene name not in index

**Solution**: Check exact spelling in Gene Metadata tab. Use gene symbols as they appear in ``adata.var_names``.

No Expression Detected
~~~~~~~~~~~~~~~~~~~~~~

**Issue**: All cells show zero expression (% Expressing = 0%)

**Possible causes**:

- Gene not expressed in this dataset
- Gene filtered during preprocessing
- Incorrect gene name

UMAP Not Available
~~~~~~~~~~~~~~~~~~

**Error**: No UMAP reduction found

**Solution**: Compute UMAP using Scanpy:

.. code-block:: python

   sc.pp.neighbors(adata)
   sc.tl.umap(adata)

Next Steps
----------

- Use :doc:`de_module` to find differentially expressed genes
- Return to :doc:`explore_module` to inspect metadata
