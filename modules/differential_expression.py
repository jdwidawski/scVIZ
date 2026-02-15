"""
Module: Differential Expression Analysis

Provides pseudobulk differential expression analysis using PyDESeq2
for single-cell RNA-seq datasets with two comparison modes:
1. Category-based: Compare metadata categories
2. Selection-based: Lasso select cells on UMAP
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from anndata import AnnData
from typing import List, Optional, Tuple
import io


def get_categorical_columns(adata: AnnData, max_categories: int = 100) -> List[str]:
    """Get list of categorical columns from obs with limited unique values."""
    categorical_cols = []
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category':
            if adata.obs[col].nunique() <= max_categories:
                categorical_cols.append(col)
    return categorical_cols


def create_pseudobulk_from_groups(
    adata: AnnData,
    sample_col: str,
    group_labels: np.ndarray,
    min_cells_per_sample: int = 10
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create pseudobulk counts by aggregating cells per sample within each group.
    
    Parameters:
        adata: AnnData object with raw counts
        sample_col: Column in obs identifying samples
        group_labels: Array of group labels ('Group1', 'Group2') for each cell
        min_cells_per_sample: Minimum cells required per sample
    
    Returns:
        Tuple of (counts_df, metadata_df)
    """
    # Add temporary group column
    adata.obs['_de_group'] = group_labels
    
    # Get unique samples
    samples = adata.obs[sample_col].unique()
    
    # Initialize counts dictionary
    counts_dict = {}
    metadata_list = []
    
    for sample in samples:
        # Get cells belonging to this sample
        sample_mask = adata.obs[sample_col] == sample
        sample_cells = adata[sample_mask]
        
        if sample_cells.n_obs < min_cells_per_sample:
            continue
        
        # Get the majority group for this sample
        group_counts = sample_cells.obs['_de_group'].value_counts()
        
        # Skip samples that don't have a clear group assignment (mixed samples)
        if len(group_counts) == 0:
            continue
        
        # Use the dominant group for this sample
        dominant_group = group_counts.index[0]
        dominant_count = group_counts.iloc[0]
        
        # Only use cells from the dominant group for cleaner pseudobulk
        group_mask = sample_cells.obs['_de_group'] == dominant_group
        group_cells = sample_cells[group_mask]
        
        if group_cells.n_obs < min_cells_per_sample:
            continue
        
        # Sum counts across cells (pseudobulk)
        if hasattr(group_cells.X, 'toarray'):
            sample_counts = np.array(group_cells.X.toarray().sum(axis=0)).flatten()
        else:
            sample_counts = np.array(group_cells.X.sum(axis=0)).flatten()
        
        counts_dict[sample] = sample_counts
        metadata_list.append({
            'sample': sample,
            'condition': dominant_group,
            'n_cells': group_cells.n_obs
        })
    
    # Clean up temporary column
    del adata.obs['_de_group']
    
    if len(counts_dict) == 0:
        raise ValueError("No samples passed the minimum cell count filter")
    
    # Create DataFrames
    counts_df = pd.DataFrame(
        counts_dict,
        index=adata.var_names
    ).T.astype(int)
    
    metadata_df = pd.DataFrame(metadata_list).set_index('sample')
    
    # Ensure condition is categorical
    metadata_df['condition'] = pd.Categorical(
        metadata_df['condition'],
        categories=['Group2', 'Group1']  # Reference first, then test
    )
    
    return counts_df, metadata_df


def run_pydeseq2(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    alpha: float = 0.05,
    lfc_shrink: bool = True,
    min_counts: int = 10
) -> pd.DataFrame:
    """
    Run PyDESeq2 differential expression analysis following the standard workflow.
    
    Parameters:
        counts_df: Pseudobulk counts (samples x genes)
        metadata_df: Sample metadata with 'condition' column
        alpha: Significance threshold
        lfc_shrink: Whether to perform LFC shrinkage
        min_counts: Minimum total counts per gene to include
    
    Returns:
        DataFrame with DE results
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    
    # Filter genes with very low counts
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= min_counts]
    counts_filtered = counts_df[genes_to_keep]
    
    if len(counts_filtered.columns) == 0:
        raise ValueError(f"No genes passed the minimum count filter (>= {min_counts} total counts)")
    
    # Create inference object
    inference = DefaultInference(n_cpus=4)
    
    # Step 1: Create DESeq dataset
    dds = DeseqDataSet(
        counts=counts_filtered,
        metadata=metadata_df,
        design="~condition",
        refit_cooks=True,
        inference=inference
    )
    
    # Step 2: Fit the model (runs all steps: size factors, dispersions, LFCs)
    dds.deseq2()
    
    # Step 3: Statistical analysis with DeseqStats
    # Get the condition levels
    condition_levels = metadata_df['condition'].cat.categories.tolist()
    
    # Contrast: Group1 vs Group2 (Group2 is reference)
    ds = DeseqStats(
        dds,
        contrast=["condition", "Group1", "Group2"],
        alpha=alpha,
        cooks_filter=True,
        independent_filter=True,
        inference=inference
    )
    
    # Run Wald test
    ds.summary()
    
    # Optionally shrink LFCs for better visualization
    if lfc_shrink:
        try:
            ds.lfc_shrink(coeff="condition[T.Group1]")
        except Exception:
            pass  # Continue without shrinkage if it fails
    
    # Get results
    results_df = ds.results_df.copy()
    
    # Add gene names as column
    results_df = results_df.reset_index()
    results_df = results_df.rename(columns={'index': 'gene'})
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('padj', ascending=True)
    
    return results_df


def volcano_plot(
    de_df: pd.DataFrame,
    pvalue_column: str = "padj",
    lfc_column: str = "log2FoldChange",
    pvalue_cutoff: float = 0.05,
    lfc_cutoff: float = 0.58,
    genes_column: str = "gene",
    height: int = 600,
    width: int = 900,
) -> go.Figure:
    """
    Make volcano plot using DE results table.
    """
    if len(de_df) == 0:
        raise ValueError("Provided `de_df` is empty.")
    
    # Ensure gene names column exists
    if genes_column not in de_df.columns:
        de_df = de_df.reset_index()
        if 'index' in de_df.columns:
            de_df['gene'] = de_df['index']
            genes_column = 'gene'
    
    # Remove rows with NaN p-values
    de_df = de_df.dropna(subset=[pvalue_column, lfc_column])
    
    # Split into up, down, and not significant
    upregulated_genes = de_df[
        (de_df[lfc_column] > lfc_cutoff) & (de_df[pvalue_column] < pvalue_cutoff)
    ]
    
    downregulated_genes = de_df[
        (de_df[lfc_column] < -lfc_cutoff) & (de_df[pvalue_column] < pvalue_cutoff)
    ]
    
    rest_genes = de_df[
        ~de_df[genes_column].isin(
            downregulated_genes[genes_column].tolist() + 
            upregulated_genes[genes_column].tolist()
        )
    ]
    
    # Create figure
    fig = go.Figure()
    
    # Upregulated genes (green)
    ys = -np.log10(upregulated_genes[pvalue_column])
    xs = upregulated_genes[lfc_column]
    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        name=f"Upregulated ({len(upregulated_genes)})",
        text=upregulated_genes[genes_column],
        mode="markers",
        marker=dict(color="#4FBA6F"),
        hovertemplate=(
            'gene: <b>%{text}</b><br>' +
            'log2(FC): <b>%{x:.3f}</b><br>' +
            '-log10(padj): <b>%{y:.3f}</b><br>' +
            '<extra></extra>'
        ),
    ))
    
    # Downregulated genes (red)
    ys = -np.log10(downregulated_genes[pvalue_column])
    xs = downregulated_genes[lfc_column]
    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        name=f"Downregulated ({len(downregulated_genes)})",
        text=downregulated_genes[genes_column],
        mode="markers",
        marker=dict(color="#F15C5A"),
        hovertemplate=(
            'gene: <b>%{text}</b><br>' +
            'log2(FC): <b>%{x:.3f}</b><br>' +
            '-log10(padj): <b>%{y:.3f}</b><br>' +
            '<extra></extra>'
        ),
    ))
    
    # Not significant genes (gray)
    ys = -np.log10(rest_genes[pvalue_column])
    xs = rest_genes[lfc_column]
    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        name=f"Not significant ({len(rest_genes)})",
        text=rest_genes[genes_column],
        mode="markers",
        marker=dict(color="Gray"),
        hovertemplate=(
            'gene: <b>%{text}</b><br>' +
            'log2(FC): <b>%{x:.3f}</b><br>' +
            '-log10(padj): <b>%{y:.3f}</b><br>' +
            '<extra></extra>'
        ),
    ))
    
    # Update marker style
    fig.update_traces(marker_line_width=0.5, marker_size=8)
    
    # Define axis limits
    absxlim = de_df[lfc_column].abs().max() * 1.2
    xlim = (-absxlim, absxlim)
    max_y = de_df[pvalue_column].apply(lambda x: -np.log10(x) if x > 0 else 0).max() * 1.1
    ylim = (0, max_y)
    
    # Add cutoff lines
    fig.add_hline(
        y=-np.log10(pvalue_cutoff),
        line_dash="dash",
        line_color="gray",
        line_width=1,
        annotation_text=f"p={pvalue_cutoff}",
        annotation_position="right"
    )
    fig.add_vline(x=lfc_cutoff, line_dash="dash", line_color="gray", line_width=1)
    fig.add_vline(x=-lfc_cutoff, line_dash="dash", line_color="gray", line_width=1)
    
    # Layout styling
    fig.update_layout(
        xaxis_range=xlim,
        yaxis_range=ylim,
        height=height,
        width=width,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(
            font=dict(size=12),
            title="DE Category",
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02
        ),
        xaxis=dict(
            title=dict(text="log2(Fold Change)", font=dict(size=14)),
            tickfont=dict(size=12),
            linecolor="black",
            linewidth=1,
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            title=dict(text="-log10(adjusted p-value)", font=dict(size=14)),
            tickfont=dict(size=12),
            linecolor="black",
            linewidth=1,
            showgrid=False,
            zeroline=False
        )
    )
    
    return fig


@st.dialog("About Differential Expression Analysis", width="large")
def show_de_info_modal():
    """Display information about differential expression analysis."""
    st.markdown("""
    ## What is Differential Expression Analysis?
    
    Differential expression (DE) analysis identifies genes that are significantly 
    up- or down-regulated between two conditions or cell populations.
    
    ### How It Works in scVIZ
    
    1. **Pseudobulk Aggregation**: Single cells are grouped by biological sample 
       (e.g., patient, donor) and their gene counts are summed. This creates 
       "pseudobulk" samples that better represent biological variability.
    
    2. **PyDESeq2 Analysis**: The pseudobulk counts are analyzed using PyDESeq2, 
       the Python implementation of the widely-used DESeq2 method:
       - Normalization with size factors
       - Dispersion estimation
       - Negative binomial model fitting
       - Wald statistical test
       - Multiple testing correction (Benjamini-Hochberg)
    
    3. **LFC Shrinkage** (optional): Log fold change estimates are "shrunk" 
       toward zero for genes with high variance, providing more reliable estimates.
    
    ### Comparison Modes
    
    - **Category-based**: Compare cells based on metadata categories (e.g., 
      disease vs. healthy, treated vs. untreated)
    - **Lasso selection**: Interactively select cell populations on the UMAP 
      to compare
    
    ### Requirements
    
    - **Sample column**: A metadata column identifying biological replicates
    - **At least 2 samples per group**: Statistical testing requires replication
    - **Raw counts**: The expression matrix should contain raw (unnormalized) counts
    
    ### Interpreting Results
    
    - **log2FoldChange**: Positive values = higher in Group 1, negative = lower
    - **padj**: Adjusted p-value (FDR-corrected). Values < 0.05 are typically significant
    - **baseMean**: Average expression level across all samples
    """)
    
    st.info("💡 **Tip**: For best results, ensure your dataset contains raw counts and has proper sample-level metadata.")


def render(adata: AnnData) -> None:
    """
    Render the Differential Expression Analysis module.
    """
    # Title with info button
    col_title, col_info = st.columns([6, 1])
    with col_title:
        st.title("🧬 Differential Expression Analysis")
    with col_info:
        st.markdown("<div style='margin-top: 12px;'></div>", unsafe_allow_html=True)
        if st.button("ℹ️ Help", key="de_help_btn", use_container_width=True):
            show_de_info_modal()
    
    st.markdown("Perform pseudobulk differential expression analysis using **PyDESeq2**")
    st.markdown("---")
    
    # Check if pydeseq2 is installed
    try:
        import pydeseq2
    except ImportError:
        st.error("⚠️ PyDESeq2 is not installed. Please install it using:")
        st.code("pip install pydeseq2")
        return
    
    # Check for UMAP (needed for selection mode)
    has_umap = 'X_umap' in adata.obsm.keys()
    
    # Get categorical columns
    categorical_cols = get_categorical_columns(adata)
    
    if len(categorical_cols) < 1:
        st.warning("⚠️ Dataset needs at least 1 categorical column in `.obs` for sample identification.")
        return
    
    # Step 1: Choose comparison mode
    st.subheader("1️⃣ Select Comparison Mode")
    
    comparison_mode = st.radio(
        "How would you like to define the groups to compare?",
        options=["Category-based comparison", "Lasso selection on UMAP"],
        horizontal=True,
        key="comparison_mode",
        help="Category-based: compare metadata categories. Lasso: select cells on UMAP."
    )
    
    st.markdown("---")
    
    # Initialize group labels
    group_labels = None
    group1_name = "Group1"
    group2_name = "Group2"
    
    if comparison_mode == "Category-based comparison":
        group_labels, group1_name, group2_name = render_category_comparison(adata, categorical_cols)
    else:
        if not has_umap:
            st.warning("⚠️ No UMAP embedding found. Please compute UMAP first in the Visualize Dataset module.")
            return
        group_labels, group1_name, group2_name = render_lasso_selection(adata)
    
    # Check if we have valid groups
    if group_labels is None:
        return
    
    # Count cells in each group
    n_group1 = np.sum(group_labels == 'Group1')
    n_group2 = np.sum(group_labels == 'Group2')
    
    st.markdown("---")
    
    # Step 2: Sample identification
    st.subheader("2️⃣ Select Sample Identifier")
    st.caption("Column that identifies biological replicates for pseudobulk aggregation")
    
    col_sample, col_info = st.columns([2, 1])
    
    with col_sample:
        sample_col = st.selectbox(
            "Sample column:",
            options=categorical_cols,
            key="sample_col",
            help="Each unique value represents a biological replicate"
        )
    
    with col_info:
        if sample_col:
            n_samples = adata.obs[sample_col].nunique()
            st.metric("Unique Samples", n_samples)
    
    # Analysis parameters
    st.markdown("---")
    st.subheader("3️⃣ Analysis Parameters")
    
    col_p1, col_p2, col_p3 = st.columns(3)
    
    with col_p1:
        min_cells = st.number_input(
            "Min cells per sample",
            min_value=1,
            max_value=100,
            value=10,
            key="min_cells",
            help="Samples with fewer cells will be excluded"
        )
    
    with col_p2:
        alpha = st.number_input(
            "Significance threshold (α)",
            min_value=0.001,
            max_value=0.1,
            value=0.05,
            step=0.01,
            key="alpha"
        )
    
    with col_p3:
        lfc_shrink = st.radio(
            "LFC shrinkage",
            options=["Yes (recommended)", "No"],
            index=0,
            key="lfc_shrink_radio",
            help="Shrink log fold change estimates for more reliable results"
        ) == "Yes (recommended)"
    
    st.markdown("---")
    
    # Summary before running
    st.subheader("📊 Comparison Summary")
    
    col_s1, col_s2 = st.columns(2)
    with col_s1:
        st.success(f"**{group1_name}** (Group1): {n_group1:,} cells")
    with col_s2:
        st.info(f"**{group2_name}** (Group2 - Reference): {n_group2:,} cells")
    
    # Run analysis button
    st.markdown("---")
    
    if st.button("🚀 Run Differential Expression Analysis", type="primary", use_container_width=True):
        run_de_analysis(adata, sample_col, group_labels, min_cells, alpha, lfc_shrink, group1_name, group2_name)
    
    # Display results if available
    if 'de_results' in st.session_state and st.session_state.de_results is not None:
        display_de_results(st.session_state.de_results, alpha, group1_name, group2_name)


def render_category_comparison(adata: AnnData, categorical_cols: List[str]) -> Tuple[Optional[np.ndarray], str, str]:
    """Render category-based comparison UI and return group labels."""
    
    st.markdown("**Select category and values to compare:**")
    
    col1, col2 = st.columns(2)
    
    with col1:
        category_col = st.selectbox(
            "Category column:",
            options=categorical_cols,
            key="category_col",
            help="The metadata column containing groups to compare"
        )
    
    if not category_col:
        return None, "", ""
    
    unique_values = adata.obs[category_col].unique().tolist()
    
    with col2:
        comparison_type = st.radio(
            "Comparison type:",
            options=["Two categories", "One category vs rest"],
            key="comparison_type",
            horizontal=True
        )
    
    if comparison_type == "Two categories":
        col_g1, col_g2 = st.columns(2)
        
        with col_g1:
            group1_value = st.selectbox(
                "Group 1 (Test):",
                options=unique_values,
                index=0,
                key="group1_value"
            )
        
        with col_g2:
            remaining_values = [v for v in unique_values if v != group1_value]
            group2_value = st.selectbox(
                "Group 2 (Reference):",
                options=remaining_values,
                index=0 if remaining_values else None,
                key="group2_value"
            )
        
        if not group2_value:
            st.warning("Need at least 2 different category values")
            return None, "", ""
        
        # Create group labels
        group_labels = np.array(['_exclude'] * adata.n_obs)
        group_labels[adata.obs[category_col] == group1_value] = 'Group1'
        group_labels[adata.obs[category_col] == group2_value] = 'Group2'
        
        return group_labels, str(group1_value), str(group2_value)
    
    else:  # One category vs rest
        group1_value = st.selectbox(
            "Select category (vs all other cells):",
            options=unique_values,
            key="group1_vs_rest"
        )
        
        # Create group labels
        group_labels = np.array(['Group2'] * adata.n_obs)  # All cells start as "rest"
        group_labels[adata.obs[category_col] == group1_value] = 'Group1'
        
        return group_labels, str(group1_value), "Rest"


def render_lasso_selection(adata: AnnData) -> Tuple[Optional[np.ndarray], str, str]:
    """Render lasso selection UI on UMAP and return group labels."""
    
    # Selection mode choice
    selection_mode = st.radio(
        "Selection mode:",
        options=["Select Group 1 vs rest", "Select both Group 1 and Group 2"],
        horizontal=True,
        key="lasso_selection_mode"
    )
    
    st.markdown("---")
    
    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']
    
    if selection_mode == "Select Group 1 vs rest":
        st.markdown("**Select cells for Group 1 (will be compared against all other cells):**")
        st.caption("Use lasso or box select to choose cells.")
        
        # Create dataframe for plotting
        df_plot = pd.DataFrame({
            'UMAP1': umap_coords[:, 0],
            'UMAP2': umap_coords[:, 1],
            'cell_idx': range(len(umap_coords))
        })
        
        # Create the plot
        fig = px.scatter(
            df_plot,
            x='UMAP1',
            y='UMAP2',
            title='Select cells for Group 1',
            width=700,
            height=500
        )
        
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(
                showgrid=False,
                zeroline=False,
                showline=True,
                linecolor='black',
                linewidth=1
            ),
            yaxis=dict(
                showgrid=False,
                zeroline=False,
                showline=True,
                linecolor='black',
                linewidth=1
            ),
            dragmode='lasso'
        )
        
        fig.update_traces(marker=dict(size=4, color='#3b82f6'))
        
        # Display chart with selection
        event = st.plotly_chart(fig, use_container_width=False, on_select="rerun", key="de_umap_select_g1")
        
        # Process selection
        selected_indices = []
        if event and event.selection and event.selection.point_indices:
            selected_indices = list(event.selection.point_indices)
        
        if len(selected_indices) == 0:
            st.info("👆 Use lasso or box select on the plot above to select cells for Group 1")
            return None, "", ""
        
        # Create group labels
        group_labels = np.array(['Group2'] * adata.n_obs)  # All cells start as "rest"
        for idx in selected_indices:
            group_labels[idx] = 'Group1'
        
        st.success(f"✅ Selected **{len(selected_indices):,}** cells for Group 1 (vs {adata.n_obs - len(selected_indices):,} remaining cells)")
        
        return group_labels, "Selected", "Rest"
    
    else:  # Select both groups
        col_g1, col_g2 = st.columns(2)
        
        with col_g1:
            st.markdown("**Select cells for Group 1 (Test):**")
            
            df_plot_g1 = pd.DataFrame({
                'UMAP1': umap_coords[:, 0],
                'UMAP2': umap_coords[:, 1],
                'cell_idx': range(len(umap_coords))
            })
            
            fig_g1 = px.scatter(
                df_plot_g1,
                x='UMAP1',
                y='UMAP2',
                title='Group 1 (Test)',
                width=450,
                height=400
            )
            
            fig_g1.update_layout(
                plot_bgcolor='white',
                paper_bgcolor='white',
                xaxis=dict(showgrid=False, zeroline=False, showline=True, linecolor='black', linewidth=1),
                yaxis=dict(showgrid=False, zeroline=False, showline=True, linecolor='black', linewidth=1),
                dragmode='lasso',
                margin=dict(l=20, r=20, t=40, b=20)
            )
            
            fig_g1.update_traces(marker=dict(size=3, color='#22c55e'))
            
            event_g1 = st.plotly_chart(fig_g1, use_container_width=True, on_select="rerun", key="de_umap_select_g1_dual")
            
            selected_g1 = []
            if event_g1 and event_g1.selection and event_g1.selection.point_indices:
                selected_g1 = list(event_g1.selection.point_indices)
            
            if selected_g1:
                st.success(f"✅ Group 1: {len(selected_g1):,} cells")
            else:
                st.info("Select cells for Group 1")
        
        with col_g2:
            st.markdown("**Select cells for Group 2 (Reference):**")
            
            df_plot_g2 = pd.DataFrame({
                'UMAP1': umap_coords[:, 0],
                'UMAP2': umap_coords[:, 1],
                'cell_idx': range(len(umap_coords))
            })
            
            fig_g2 = px.scatter(
                df_plot_g2,
                x='UMAP1',
                y='UMAP2',
                title='Group 2 (Reference)',
                width=450,
                height=400
            )
            
            fig_g2.update_layout(
                plot_bgcolor='white',
                paper_bgcolor='white',
                xaxis=dict(showgrid=False, zeroline=False, showline=True, linecolor='black', linewidth=1),
                yaxis=dict(showgrid=False, zeroline=False, showline=True, linecolor='black', linewidth=1),
                dragmode='lasso',
                margin=dict(l=20, r=20, t=40, b=20)
            )
            
            fig_g2.update_traces(marker=dict(size=3, color='#3b82f6'))
            
            event_g2 = st.plotly_chart(fig_g2, use_container_width=True, on_select="rerun", key="de_umap_select_g2_dual")
            
            selected_g2 = []
            if event_g2 and event_g2.selection and event_g2.selection.point_indices:
                selected_g2 = list(event_g2.selection.point_indices)
            
            if selected_g2:
                st.info(f"✅ Group 2: {len(selected_g2):,} cells")
            else:
                st.info("Select cells for Group 2")
        
        # Check for overlapping selections
        overlap = set(selected_g1) & set(selected_g2)
        if overlap:
            st.warning(f"⚠️ {len(overlap)} cells selected in both groups. These will be excluded.")
        
        if len(selected_g1) == 0 or len(selected_g2) == 0:
            st.info("👆 Select cells in both plots to define the comparison groups")
            return None, "", ""
        
        # Create group labels
        group_labels = np.array(['_exclude'] * adata.n_obs)
        for idx in selected_g1:
            if idx not in overlap:
                group_labels[idx] = 'Group1'
        for idx in selected_g2:
            if idx not in overlap:
                group_labels[idx] = 'Group2'
        
        return group_labels, "Selection 1", "Selection 2"


def run_de_analysis(
    adata: AnnData,
    sample_col: str,
    group_labels: np.ndarray,
    min_cells: int,
    alpha: float,
    lfc_shrink: bool,
    group1_name: str,
    group2_name: str
) -> None:
    """Run the DE analysis pipeline."""
    
    with st.spinner("Creating pseudobulk samples..."):
        try:
            # Filter to only cells in groups (exclude '_exclude' cells)
            valid_mask = group_labels != '_exclude'
            adata_filtered = adata[valid_mask].copy()
            group_labels_filtered = group_labels[valid_mask]
            
            counts_df, metadata_df = create_pseudobulk_from_groups(
                adata_filtered,
                sample_col,
                group_labels_filtered,
                min_cells_per_sample=min_cells
            )
            
            # Check we have enough samples per group
            group_counts = metadata_df['condition'].value_counts()
            
            st.write("**Pseudobulk samples created:**")
            for group, count in group_counts.items():
                group_name = group1_name if group == 'Group1' else group2_name
                st.write(f"  - {group_name}: {count} samples")
            
            if len(group_counts) < 2:
                st.error("⚠️ Need samples in both groups. Check your sample column.")
                return
            
            for group, count in group_counts.items():
                if count < 2:
                    st.error(f"⚠️ Need at least 2 samples per group. {group} has only {count}.")
                    return
            
        except Exception as e:
            st.error(f"❌ Error creating pseudobulk: {str(e)}")
            return
    
    with st.spinner("Running PyDESeq2 analysis..."):
        try:
            results_df = run_pydeseq2(
                counts_df,
                metadata_df,
                alpha=alpha,
                lfc_shrink=lfc_shrink
            )
            
            # Store results in session state
            st.session_state.de_results = results_df
            st.session_state.de_group1_name = group1_name
            st.session_state.de_group2_name = group2_name
            
            st.success(f"✅ Analysis complete! Found {len(results_df)} genes tested.")
            
        except Exception as e:
            st.error(f"❌ Error running PyDESeq2: {str(e)}")
            return


def format_scientific(val):
    """Format a value in scientific notation for display."""
    if pd.isna(val) or val == 0:
        return "0"
    elif val < 0.0001:
        return f"{val:.2e}"
    else:
        return f"{val:.4f}"


def display_de_results(results_df: pd.DataFrame, alpha: float, group1_name: str, group2_name: str) -> None:
    """Display DE results with volcano plot and table."""
    
    st.markdown("---")
    st.subheader("📈 Results")
    
    # Significance cutoffs for volcano
    col_cut1, col_cut2 = st.columns(2)
    with col_cut1:
        lfc_cutoff = st.slider(
            "Log2 Fold Change cutoff",
            min_value=0.0,
            max_value=10.0,
            value=0.58,
            step=0.1,
            key="lfc_cutoff"
        )
    with col_cut2:
        pval_cutoff = st.select_slider(
            "Adjusted p-value cutoff",
            options=[0.001, 0.005, 0.01, 0.025, 0.05, 0.1],
            value=0.05,
            key="pval_cutoff"
        )
    
    # Count significant genes
    sig_up = len(results_df[
        (results_df['log2FoldChange'] > lfc_cutoff) & 
        (results_df['padj'] < pval_cutoff)
    ])
    sig_down = len(results_df[
        (results_df['log2FoldChange'] < -lfc_cutoff) & 
        (results_df['padj'] < pval_cutoff)
    ])
    
    col_m1, col_m2, col_m3 = st.columns(3)
    with col_m1:
        st.metric("Total Genes", len(results_df))
    with col_m2:
        st.metric(f"⬆️ Up in {group1_name}", sig_up)
    with col_m3:
        st.metric(f"⬇️ Down in {group1_name}", sig_down)
    
    # Create tabs for volcano and table
    tab_volcano, tab_table = st.tabs(["🌋 Volcano Plot", "📋 Results Table"])
    
    with tab_volcano:
        try:
            fig = volcano_plot(
                results_df,
                pvalue_cutoff=pval_cutoff,
                lfc_cutoff=lfc_cutoff
            )
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"Positive log2FC = higher expression in **{group1_name}** vs **{group2_name}**")
        except Exception as e:
            st.error(f"Could not create volcano plot: {str(e)}")
    
    with tab_table:
        # Filter options
        st.markdown("**Filter results:**")
        col_f1, col_f2 = st.columns(2)
        
        with col_f1:
            show_significant = st.checkbox("Show only significant genes", value=False, key="show_sig")
        with col_f2:
            gene_search = st.text_input("Search gene:", key="gene_search")
        
        # Apply filters
        display_df = results_df.copy()
        
        if show_significant:
            display_df = display_df[
                (display_df['padj'] < pval_cutoff) & 
                (display_df['log2FoldChange'].abs() > lfc_cutoff)
            ]
        
        if gene_search:
            display_df = display_df[
                display_df['gene'].str.contains(gene_search, case=False, na=False)
            ]
        
        # Prepare display dataframe with scientific notation for p-values
        table_df = display_df[['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj']].copy()
        
        # Keep numeric values for sorting but format for display
        st.dataframe(
            table_df,
            use_container_width=True,
            height=400,
            column_config={
                "gene": st.column_config.TextColumn("Gene"),
                "baseMean": st.column_config.NumberColumn("baseMean", format="%.2f"),
                "log2FoldChange": st.column_config.NumberColumn("log2FC", format="%.3f"),
                "lfcSE": st.column_config.NumberColumn("lfcSE", format="%.3f"),
                "pvalue": st.column_config.NumberColumn("pvalue", format="%.2e"),
                "padj": st.column_config.NumberColumn("padj", format="%.2e")
            }
        )
        
        # Download button
        csv = results_df.to_csv(index=False)
        st.download_button(
            label="📥 Download full results (CSV)",
            data=csv,
            file_name="de_results.csv",
            mime="text/csv",
            use_container_width=True
        )
