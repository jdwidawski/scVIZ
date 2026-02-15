Differential Expression Analysis
=================================

The Differential Expression module performs pseudobulk differential expression analysis using PyDESeq2 to identify genes with significant expression changes between conditions.

Overview
--------

This module:

1. Aggregates single cells into pseudobulk samples
2. Runs PyDESeq2 analysis (DESeq2 in Python)
3. Generates interactive volcano plots
4. Provides downloadable results tables

The analysis follows established best practices for scRNA-seq differential expression by aggregating cells to the sample level before statistical testing.

Prerequisites
-------------

Your dataset requires:

**Sample Identifiers**
   A column in ``adata.obs`` identifying biological replicates (e.g., "sample_id", "donor")

**Condition Labels**
   A column specifying experimental conditions to compare (e.g., "disease", "treatment")

**Replicates**
   At least 2 samples per condition (3+ recommended)

Configuration
-------------

Sample Identifier Column
~~~~~~~~~~~~~~~~~~~~~~~~

Select the column that defines biological samples:

- Each unique value represents one pseudobulk sample
- Must have multiple cells per sample
- Should reflect true biological replicates

Example::

   sample_id: donor1, donor2, donor3, donor4
   (NOT: cell1, cell2, cell3, cell4)

Condition Column
~~~~~~~~~~~~~~~~

Select the column defining groups to compare:

- Should have exactly 2 unique values (e.g., "control", "treatment")
- Must be associated with sample identifiers
- All cells from same sample must have same condition

Example::

   condition: control, treatment
   disease_status: healthy, diseased

Minimum Cells per Sample
~~~~~~~~~~~~~~~~~~~~~~~~~

Filter samples with too few cells:

- Default: 10 cells minimum
- Higher values (50-100) recommended for robust pseudobulk
- Samples below threshold are excluded

Analysis Parameters
-------------------

Reference Condition
~~~~~~~~~~~~~~~~~~~

Choose the baseline condition for comparison:

- Used as denominator in fold-change calculations
- Typically "control", "healthy", or "untreated"
- Affects sign of log2 fold changes

Log2 Fold Change Shrinkage
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Enable Wald statistic shrinkage:

- **On** (default): Reduces noise in fold-change estimates
- **Off**: Uses raw fold changes
- Recommended for datasets with few replicates

Significance Thresholds
~~~~~~~~~~~~~~~~~~~~~~~~

Set cutoffs for calling significant genes:

**Adjusted p-value**
   Default: 0.05 (Benjamini-Hochberg FDR correction)

**Log2 Fold Change**
   Default: 1.0 (2-fold change)

Running Analysis
----------------

1. Configure sample and condition columns
2. Set analysis parameters
3. Click "Run Differential Expression Analysis"
4. Wait for computation (may take 1-5 minutes)

The module will:

- Create pseudobulk count matrix
- Normalize samples with DESeq2 size factors
- Estimate dispersions
- Test for differential expression
- Apply shrinkage (if enabled)
- Adjust p-values for multiple testing

Results
-------

Volcano Plot
~~~~~~~~~~~~

Interactive Plotly scatter plot showing:

- **X-axis**: Log2 fold change
- **Y-axis**: -log10(adjusted p-value)
- **Colors**: 
  
  - Red: Significantly upregulated
  - Blue: Significantly downregulated  
  - Gray: Not significant

**Features:**

- Hover to see gene names and statistics
- Zoom and pan
- Download as PNG

Results Table
~~~~~~~~~~~~~

Comprehensive table with columns:

- **gene**: Gene name/ID
- **baseMean**: Average normalized expression
- **log2FoldChange**: Effect size (log2 scale)
- **lfcSE**: Standard error of log2 fold change
- **pvalue**: Raw p-value
- **padj**: Adjusted p-value (FDR)

**Filtering:**

- Search for specific genes
- Sort by any column
- Filter by significance

**Export:**

Download results as CSV file for further analysis in R, Python, or Excel.

Interpreting Results
--------------------

Log2 Fold Change
~~~~~~~~~~~~~~~~

Represents the magnitude of expression change:

- **Positive**: Upregulated in test condition vs reference
- **Negative**: Downregulated in test condition vs reference
- **±1**: 2-fold change
- **±2**: 4-fold change

Example: log2FC = 2.5 means gene is ~5.7× higher in test condition.

Adjusted P-value
~~~~~~~~~~~~~~~~

Probability that the result is a false positive:

- **< 0.05**: Statistically significant (5% FDR)
- **< 0.01**: Highly significant (1% FDR)
- Uses Benjamini-Hochberg correction for multiple testing

Base Mean
~~~~~~~~~

Average normalized expression across all samples:

- Higher values indicate more abundant genes
- Low baseMean genes may have noisy fold changes
- Consider filtering lowly expressed genes

Tips
----

✅ **Use Biological Replicates**
   Need ≥3 samples per condition for reliable statistics

✅ **Check Sample Balance**
   Similar number of samples per condition improves power

✅ **Enable Shrinkage**
   Especially important with few replicates (3-5 samples)

✅ **Validate Top Hits**
   Plot expression of significant genes in Visualize module

✅ **Consider Batch Effects**
   PyDESeq2 doesn't handle covariates; pre-correct batches if needed

Common Use Cases
----------------

Disease vs Healthy
~~~~~~~~~~~~~~~~~~

Compare gene expression between disease states:

1. Sample column: "patient_id"
2. Condition column: "disease_status" (healthy, diseased)
3. Reference: "healthy"

Treatment vs Control
~~~~~~~~~~~~~~~~~~~~

Identify treatment effects:

1. Sample column: "sample_id"
2. Condition column: "treatment" (control, drug)
3. Reference: "control"

Time Series
~~~~~~~~~~~

For two time points:

1. Sample column: "animal_id"
2. Condition column: "timepoint" (day0, day7)
3. Reference: "day0"

Troubleshooting
---------------

Too Few Samples
~~~~~~~~~~~~~~~

**Error**: Need at least 2 samples per condition

**Solution**: 

- Lower "Minimum cells per sample" threshold
- Check that sample IDs are correct (not cell IDs)
- Combine technical replicates if appropriate

Singular Matrix Error
~~~~~~~~~~~~~~~~~~~~~

**Error**: Design matrix is singular

**Possible causes**:

- Sample and condition perfectly confounded
- All samples in one condition have same metadata
- Insufficient replicates

**Solution**: Check that multiple samples exist per condition with proper replicates.

All Genes Non-significant
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Issue**: No genes pass significance thresholds

**Possible causes**:

- Weak biological effect
- High variability between replicates
- Insufficient power (too few samples)

**Solutions**:

- Lower significance thresholds (p < 0.1, LFC > 0.5)
- Increase sample size
- Check data quality and normalization

Analysis Takes Too Long
~~~~~~~~~~~~~~~~~~~~~~~

**Issue**: Analysis runs for >10 minutes

**Possible causes**:

- Very large gene count
- Many samples

**Solutions**:

- Filter lowly expressed genes before loading
- Use smaller subset for testing
- Be patient (complex datasets can take time)

Next Steps
----------

- Visualize significant genes in :doc:`visualize_module`
- Export results for pathway analysis
- Return to :doc:`explore_module` to verify sample metadata
