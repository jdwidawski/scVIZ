Quality Control Module
======================

The Quality Control (QC) module provides comprehensive tools for assessing
and visualizing single-cell data quality metrics.

Overview
--------

This module enables you to:

- **Compute QC Metrics**: Calculate total counts, genes detected, and mitochondrial percentage
- **Cell Cycle Scoring**: Score cells by their cell cycle phase (S, G2M)
- **Visualize on UMAP**: Color UMAP plots by QC metrics with interactive selection
- **Distribution Analysis**: View violin plots grouped by categorical metadata

.. note::
   The QC module supports both gene symbols (e.g., CDK1, TOP2A) and Ensembl IDs
   (e.g., ENSG00000012048) for gene matching.

Accessing the Module
--------------------

1. Load a dataset from the Home screen
2. Navigate to **Quality Control** in the sidebar

Computing QC Metrics
--------------------

QC metrics are computed using Scanpy's ``sc.pp.calculate_qc_metrics`` function:

- **n_genes_by_counts**: Number of genes with positive counts per cell
- **total_counts**: Total UMI counts per cell  
- **pct_counts_mt**: Percentage of counts mapping to mitochondrial genes

The module automatically detects mitochondrial genes by looking for gene names
starting with "MT-" or "mt-".

Cell Cycle Scoring
------------------

Cell cycle scoring uses canonical marker genes for S-phase and G2M-phase:

**S-phase markers**: MCM5, PCNA, TYMS, FEN1, MCM2, MCM4, RRM1, UNG, GINS2, 
MCM6, CDCA7, DTL, PRIM1, UHRF1, MLF1IP, HELLS, RFC2, RPA2, NASP, RAD51AP1, 
GMNN, WDR76, SLBP, CCNE2, UBR7, POLD3, MSH2, ATAD2, RAD51, RRM2, CDC45, 
CDC6, EXO1, TIPIN, DSCC1, BLM, CASP8AP2, USP1, CLSPN, POLA1, CHAF1B, 
BRIP1, E2F8

**G2M-phase markers**: HMGB2, CDK1, NUSAP1, UBE2C, BIRC5, TPX2, TOP2A, NDC80,
CKS2, NUF2, CKS1B, MKI67, TMPO, CENPF, TACC3, FAM64A, SMC4, CCNB2, CKAP2L, 
CKAP2, AURKB, BUB1, KIF11, ANP32E, TUBB4B, GTSE1, KIF20B, HJURP, CDCA3, 
HN1, CDC20, TTK, CDC25C, KIF2C, RANGAP1, NCAPD2, DLGAP5, CDCA2, CDCA8, 
ECT2, KIF23, HMMR, AURKA, PSRC1, ANLN, LBR, CKAP5, CENPE, CTCF, NEK2, 
G2E3, GAS2L3, CBX5, CENPA

Visualization Options
---------------------

UMAP Visualization
~~~~~~~~~~~~~~~~~~

- **Color metric**: Select which QC metric to display
- **Color scale**: Choose colormap (viridis, plasma, inferno, magma, etc.)
- **Color bounds**: Set as percentiles (p1 to p99.5) or absolute values
- **Z-ordering**: Higher values plotted on top for better visibility
- **Lasso selection**: Select cells interactively to view statistics

Violin Plots
~~~~~~~~~~~~

- **Group by**: Split violins by categorical metadata (e.g., cell type)
- **Box plot**: Optional box plot overlay
- **Category palette**: 12 different color schemes available:
  - Plotly, D3, G10, T10, Alphabet, Dark24, Light24, Set1, Set2, Set3, Pastel, Bold
- **Statistics table**: View summary statistics per group

Figure Settings
~~~~~~~~~~~~~~~

Configure plot appearance:

- Width and height (pixels)
- Dot size for scatter plots
- Line width for violin plots

Example Workflow
----------------

1. **Load data**: Start with a dataset containing raw counts
2. **Compute metrics**: Click "Compute QC Metrics" to add quality metrics
3. **Explore UMAP**: View cells colored by metrics, identify outliers
4. **Group analysis**: Use violin plots to compare metrics across cell types
5. **Filter cells**: Use insights to set thresholds for downstream filtering

API Reference
-------------

The QC module uses these key functions:

.. code-block:: python

   import scanpy as sc
   
   # Compute QC metrics
   sc.pp.calculate_qc_metrics(
       adata, 
       qc_vars=['mt'], 
       percent_top=None, 
       log1p=False, 
       inplace=True
   )
   
   # Score cell cycle
   sc.tl.score_genes_cell_cycle(
       adata, 
       s_genes=S_GENES, 
       g2m_genes=G2M_GENES
   )

See Also
--------

- :doc:`visualize_module` - General dataset visualization
- :doc:`de_module` - Differential expression analysis
- `Scanpy QC tutorial <https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html>`_
