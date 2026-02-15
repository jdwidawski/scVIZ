"""
Module: Visualize Dataset

Provides interactive UMAP visualizations for both metadata categories 
and gene expression with customizable styling options.
"""
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from anndata import AnnData
import scanpy as sc
from typing import List, Optional, Tuple


# Available colormaps for visualization
CATEGORICAL_COLORMAPS = [
    'tab10', 'tab20', 'Set1', 'Set2', 'Set3', 'Paired', 
    'Accent', 'Dark2', 'Pastel1', 'Pastel2'
]

CONTINUOUS_COLORMAPS = [
    'Viridis', 'Plasma', 'Inferno', 'Magma', 'Cividis',
    'Blues', 'Reds', 'Greens', 'Oranges', 'Purples',
    'YlOrRd', 'YlGnBu', 'RdBu', 'Spectral'
]


def get_categorical_columns(adata: AnnData, max_categories: int = 50) -> List[str]:
    """Get list of categorical columns from obs with limited unique values."""
    categorical_cols = []
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category':
            if adata.obs[col].nunique() <= max_categories:
                categorical_cols.append(col)
    return categorical_cols


def search_gene_by_symbol(adata: AnnData, search_term: str, symbol_column: Optional[str] = None) -> List[Tuple[str, str]]:
    """
    Search for genes by name or alternative symbol.
    
    Parameters:
        adata: AnnData object
        search_term: Search query
        symbol_column: Alternative gene symbol column in var (optional)
    
    Returns:
        List of tuples (gene_id, display_name)
    """
    results = []
    search_lower = search_term.lower()
    
    # Search in main index (var_names)
    for gene in adata.var_names:
        if search_lower in gene.lower():
            results.append((gene, gene))
    
    # Search in alternative symbol column if provided
    if symbol_column and symbol_column in adata.var.columns:
        for idx, symbol in enumerate(adata.var[symbol_column]):
            if pd.notna(symbol) and search_lower in str(symbol).lower():
                gene_id = adata.var_names[idx]
                if gene_id not in [r[0] for r in results]:
                    results.append((gene_id, f"{gene_id} ({symbol})"))
    
    return results[:100]  # Limit results


def get_gene_expression(adata: AnnData, gene_name: str) -> np.ndarray:
    """Get expression values for a gene."""
    if gene_name not in adata.var_names:
        return None
    
    gene_idx = adata.var_names.get_loc(gene_name)
    
    if hasattr(adata.X, 'toarray'):
        expression = adata.X[:, gene_idx].toarray().flatten()
    else:
        expression = np.array(adata.X[:, gene_idx]).flatten()
    
    return expression


@st.dialog("About Visualize Dataset", width="large")
def show_visualize_dataset_info_modal():
    """Display information about the Visualize Dataset module."""
    st.markdown("""
    ## What is Visualize Dataset?
    
    This module provides interactive UMAP visualizations for exploring your 
    single-cell dataset by metadata categories or gene expression.
    
    ### UMAP Visualization
    
    **UMAP** (Uniform Manifold Approximation and Projection) is a dimensionality 
    reduction technique that projects high-dimensional gene expression data into 
    2D space, placing similar cells close together.
    
    ### Visualization Tabs
    
    #### 📊 Metadata Tab
    - Color cells by categorical metadata (cell type, batch, condition)
    - Explore cluster distributions and sample composition
    - Customize colormap and point styling
    
    #### 🧬 Gene Expression Tab
    - Search for genes by name or symbol
    - Visualize expression levels on UMAP
    - Identify expression patterns across cell populations
    
    ### Interactive Selection (Lasso Tool)
    
    Both tabs support **lasso selection**:
    1. Use the lasso or box select tool on the plot
    2. Selected cells are highlighted and summarized
    3. View statistics for the selected population
    
    ### Customization Options
    
    - **Point size**: Adjust visibility for large/small datasets
    - **Colormap**: Choose from categorical or continuous palettes
    - **Opacity**: Balance visibility with density perception
    """)
    
    st.info("💡 **Tip**: Use lasso selection to interactively explore cell populations and identify interesting subgroups.")


def render(adata: AnnData) -> None:
    """
    Render the Visualize Dataset module.
    
    Parameters:
        adata (AnnData): The loaded AnnData object to visualize.
    """
    # Title with info button
    col_title, col_info = st.columns([6, 1])
    with col_title:
        st.title("🗺️ Visualize Dataset")
    with col_info:
        st.markdown("<div style='margin-top: 12px;'></div>", unsafe_allow_html=True)
        if st.button("ℹ️ Help", key="viz_dataset_help_btn", use_container_width=True):
            show_visualize_dataset_info_modal()
    
    st.markdown("---")
    
    # Check for UMAP
    has_umap = 'X_umap' in adata.obsm.keys()
    
    if not has_umap:
        st.warning("⚠️ No UMAP embedding found in the dataset.")
        
        if st.button("🔄 Compute UMAP", type="primary"):
            with st.spinner("Computing PCA and UMAP..."):
                if 'X_pca' not in adata.obsm.keys():
                    sc.pp.pca(adata)
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)
            st.success("✅ UMAP computed successfully!")
            st.rerun()
        return
    
    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']
    
    # Information about dataset
    st.info(f"📊 Dataset: **{adata.n_obs:,}** cells × **{adata.n_vars:,}** genes")
    
    st.markdown("---")
    
    # Create tabs for metadata and gene expression visualization
    tab1, tab2 = st.tabs(["📋 Metadata Visualization", "🧬 Gene Expression"])
    
    with tab1:
        render_metadata_visualization(adata, umap_coords)
    
    with tab2:
        render_gene_expression_visualization(adata, umap_coords)


def render_metadata_visualization(adata: AnnData, umap_coords: np.ndarray) -> None:
    """Render metadata visualization on UMAP."""
    st.subheader("Metadata Visualization")
    
    # Get categorical columns
    categorical_cols = get_categorical_columns(adata)
    
    if not categorical_cols:
        st.warning("No categorical metadata columns available for visualization.")
        return
    
    # Settings sidebar in expander
    with st.expander("⚙️ Figure Settings", expanded=False):
        col_s1, col_s2, col_s3 = st.columns(3)
        
        with col_s1:
            width_meta = st.slider("Width (px)", 300, 1200, 700, 50, key="width_meta")
            height_meta = st.slider("Height (px)", 300, 1200, 600, 50, key="height_meta")
        
        with col_s2:
            dot_size_meta = st.slider("Dot size", 1, 20, 4, 1, key="dotsize_meta")
            colormap_meta = st.selectbox("Colormap", CATEGORICAL_COLORMAPS, index=1, key="cmap_meta")
    
    # Category selection
    selected_category = st.selectbox(
        "Select metadata category to visualize:",
        categorical_cols,
        key="meta_category"
    )
    
    if selected_category:
        # Create the plot
        df_plot = pd.DataFrame({
            'UMAP1': umap_coords[:, 0],
            'UMAP2': umap_coords[:, 1],
            'category': adata.obs[selected_category].astype(str).values
        })
        # Add original index for selection tracking
        df_plot['cell_idx'] = range(len(df_plot))
        
        # Get unique categories and assign colors
        unique_cats = df_plot['category'].unique()
        n_colors = len(unique_cats)
        
        # Create color sequence from colormap
        try:
            cmap = plt.cm.get_cmap(colormap_meta)
            colors = [f'rgb({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)})' 
                     for c in [cmap(i/max(n_colors-1, 1)) for i in range(n_colors)]]
        except:
            colors = None
        
        fig = px.scatter(
            df_plot,
            x='UMAP1',
            y='UMAP2',
            color='category',
            color_discrete_sequence=colors,
            title=f'UMAP - {selected_category}',
            width=width_meta,
            height=height_meta,
            custom_data=['cell_idx']
        )
        
        # Update layout for white background, no grid
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
            showlegend=True,
            legend=dict(
                title=selected_category,
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02
            ),
            dragmode='lasso'  # Enable lasso selection by default
        )
        
        fig.update_traces(marker=dict(size=dot_size_meta))
        
        # Display chart with selection enabled
        st.caption("💡 Use lasso or box select to analyze a subset of cells")
        event = st.plotly_chart(fig, use_container_width=False, on_select="rerun", key="meta_umap_select")
        
        # Process selection
        selected_indices = []
        if event and event.selection and event.selection.point_indices:
            # point_indices is a list of indices
            selected_indices = list(event.selection.point_indices)
        
        # Show category statistics
        st.markdown("### 📊 Category Statistics")
        
        if selected_indices:
            # Filter to selected cells
            selected_categories = df_plot.iloc[selected_indices]['category']
            cat_counts = selected_categories.value_counts()
            st.info(f"📌 **Showing statistics for {len(selected_indices):,} selected cells** (use lasso/box select on the plot)")
        else:
            cat_counts = adata.obs[selected_category].value_counts()
        
        st.dataframe(
            pd.DataFrame({
                'Category': cat_counts.index,
                'Count': cat_counts.values,
                'Percentage': (cat_counts.values / cat_counts.sum() * 100).round(2)
            }),
            use_container_width=True
        )


def render_gene_expression_visualization(adata: AnnData, umap_coords: np.ndarray) -> None:
    """Render gene expression visualization on UMAP."""
    st.subheader("Gene Expression Visualization")
    
    # Settings in expander
    with st.expander("⚙️ Figure Settings", expanded=False):
        col_s1, col_s2, col_s3 = st.columns(3)
        
        with col_s1:
            width_gene = st.slider("Width (px)", 300, 1200, 500, 50, key="width_gene")
            height_gene = st.slider("Height (px)", 300, 1200, 450, 50, key="height_gene")
        
        with col_s2:
            dot_size_gene = st.slider("Dot size", 1, 20, 4, 1, key="dotsize_gene")
            colormap_gene = st.selectbox("Colormap", CONTINUOUS_COLORMAPS, index=6, key="cmap_gene")
        
        with col_s3:
            # Color range mode
            color_mode = st.radio(
                "Color range mode:",
                ["Auto", "Percentile", "Absolute"],
                key="color_mode",
                help="Auto: use data min/max. Percentile: use expression percentiles. Absolute: specify exact values."
            )
            
            if color_mode == "Percentile":
                vmin_pct = st.number_input("Min percentile", 0.0, 100.0, 1.0, 0.5, key="vmin_pct",
                                          help="e.g., 1 for 1st percentile")
                vmax_pct = st.number_input("Max percentile", 0.0, 100.0, 99.0, 0.5, key="vmax_pct",
                                          help="e.g., 99 for 99th percentile")
                vmin_abs, vmax_abs = None, None
            elif color_mode == "Absolute":
                vmin_abs = st.number_input("Color min", value=0.0, key="vmin_gene")
                vmax_abs = st.number_input("Color max", value=0.0, key="vmax_gene")
                vmin_pct, vmax_pct = None, None
            else:  # Auto
                vmin_abs, vmax_abs, vmin_pct, vmax_pct = None, None, None, None
    
    # Alternative gene symbol column selection
    col_opt1, col_opt2 = st.columns([1, 1])
    
    with col_opt1:
        var_columns = ['(None - use gene index)'] + list(adata.var.columns)
        alt_symbol_col = st.selectbox(
            "Gene symbol column for display:",
            var_columns,
            index=0,
            key="alt_symbol",
            help="Select a column from var to display alternative gene symbols"
        )
        
        if alt_symbol_col == '(None - use gene index)':
            alt_symbol_col = None
    
    # Button to view var table
    with col_opt2:
        st.markdown("<div style='margin-top: 26px;'></div>", unsafe_allow_html=True)
        if st.button("📋 View & Search Gene Table (.var)", use_container_width=True):
            show_var_search_modal(adata)
    
    st.markdown("---")
    
    # Gene selection - multiselect from gene list
    st.markdown("**Select genes to visualize:**")
    
    # Create gene options with optional symbols
    gene_options = []
    gene_display_map = {}  # Maps display name to gene_id
    
    for gene_id in adata.var_names:
        if alt_symbol_col and alt_symbol_col in adata.var.columns:
            symbol = adata.var.loc[gene_id, alt_symbol_col]
            if pd.notna(symbol) and str(symbol).strip():
                display_name = f"{gene_id} ({symbol})"
            else:
                display_name = gene_id
        else:
            display_name = gene_id
        gene_options.append(display_name)
        gene_display_map[display_name] = gene_id
    
    selected_genes_display = st.multiselect(
        "Choose one or more genes:",
        options=gene_options,
        default=None,
        key="gene_multiselect",
        help="Type to search, select multiple genes to compare"
    )
    
    # Convert display names back to gene IDs
    selected_genes = [gene_display_map[d] for d in selected_genes_display]
    
    if selected_genes:
        st.markdown("---")
        
        # Create tabs for each gene
        gene_tabs = st.tabs([d for d in selected_genes_display])
        
        for idx, (gene_id, tab) in enumerate(zip(selected_genes, gene_tabs)):
            with tab:
                render_single_gene_plot(
                    adata, umap_coords, gene_id, alt_symbol_col,
                    width_gene, height_gene, dot_size_gene, colormap_gene,
                    color_mode,
                    vmin_pct if color_mode == "Percentile" else None,
                    vmax_pct if color_mode == "Percentile" else None,
                    vmin_abs if color_mode == "Absolute" else None,
                    vmax_abs if color_mode == "Absolute" else None,
                    show_stats=True,
                    plot_key=f"gene_plot_{gene_id}_{idx}"
                )
    else:
        st.info("👆 Select one or more genes above to visualize expression on UMAP")


@st.dialog("Gene Table (.var)", width="large")
def show_var_search_modal(adata: AnnData):
    """Display searchable var table in a modal."""
    st.markdown(f"**{adata.n_vars:,} genes** in dataset")
    st.caption("Use the search box below to filter genes. Copy gene names to use in the selection above.")
    
    # Search filter
    search_filter = st.text_input("🔍 Search genes:", placeholder="Type to filter...", key="var_search_filter")
    
    var_df = adata.var.reset_index()
    if 'index' in var_df.columns:
        var_df = var_df.rename(columns={'index': 'Gene ID'})
    
    # Apply filter
    if search_filter:
        mask = var_df.apply(lambda row: row.astype(str).str.contains(search_filter, case=False).any(), axis=1)
        var_df = var_df[mask]
        st.markdown(f"**{len(var_df)}** genes matching '{search_filter}'")
    
    st.dataframe(
        var_df,
        use_container_width=True,
        height=500
    )


def render_single_gene_plot(
    adata: AnnData,
    umap_coords: np.ndarray,
    gene_id: str,
    alt_symbol_col: Optional[str],
    width: int,
    height: int,
    dot_size: int,
    colormap: str,
    color_mode: str,
    vmin_pct: Optional[float] = None,
    vmax_pct: Optional[float] = None,
    vmin_abs: Optional[float] = None,
    vmax_abs: Optional[float] = None,
    show_stats: bool = False,
    plot_key: str = "gene_plot"
) -> None:
    """Render a single gene expression UMAP plot with lasso selection support."""
    expression = get_gene_expression(adata, gene_id)
    
    if expression is None:
        st.error(f"Could not retrieve expression for gene '{gene_id}'")
        return
    
    # Create dataframe for plotting - keep original indices for selection
    df_plot = pd.DataFrame({
        'UMAP1': umap_coords[:, 0],
        'UMAP2': umap_coords[:, 1],
        'expression': expression,
        'original_idx': range(len(expression))
    })
    
    # Sort by expression so that cells with higher expression are plotted on top
    df_plot = df_plot.sort_values('expression', ascending=True).reset_index(drop=True)
    
    # Determine color range based on mode
    range_color = None
    if color_mode == "Percentile" and vmin_pct is not None and vmax_pct is not None:
        range_color = [
            np.percentile(expression, vmin_pct),
            np.percentile(expression, vmax_pct)
        ]
    elif color_mode == "Absolute" and (vmin_abs is not None or vmax_abs is not None):
        range_color = [
            vmin_abs if vmin_abs and vmin_abs != 0 else expression.min(),
            vmax_abs if vmax_abs and vmax_abs != 0 else expression.max()
        ]
    else:  # Auto mode: 0.0 to 99.5 percentile
        range_color = [
            0.0,
            np.percentile(expression, 99.5)
        ]
    
    # Get gene display name with symbol if available
    display_name = gene_id
    if alt_symbol_col and alt_symbol_col in adata.var.columns:
        symbol = adata.var.loc[gene_id, alt_symbol_col]
        if pd.notna(symbol):
            display_name = f"{gene_id} ({symbol})"
    
    fig = px.scatter(
        df_plot,
        x='UMAP1',
        y='UMAP2',
        color='expression',
        color_continuous_scale=colormap,
        range_color=range_color,
        title=f'{display_name}',
        width=width,
        height=height,
        custom_data=['original_idx']
    )
    
    # Update layout for white background, no grid
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
        coloraxis_colorbar=dict(
            title="Expr"
        ),
        margin=dict(l=20, r=20, t=40, b=20),
        dragmode='lasso'  # Enable lasso selection by default
    )
    
    fig.update_traces(marker=dict(size=dot_size))
    
    # Display chart with selection enabled
    st.caption("💡 Use lasso or box select to analyze a subset of cells")
    event = st.plotly_chart(fig, use_container_width=False, on_select="rerun", key=plot_key)
    
    # Process selection
    selected_indices = []
    if event and event.selection and event.selection.point_indices:
        # point_indices is a list of indices
        selected_indices = list(event.selection.point_indices)
    
    # Expression statistics (optional)
    if show_stats:
        if selected_indices:
            # Get expression values for selected cells (map back through sorted dataframe)
            selected_expression = df_plot.iloc[selected_indices]['expression'].values
            st.info(f"📌 **Showing statistics for {len(selected_indices):,} selected cells** (use lasso/box select on the plot)")
        else:
            selected_expression = expression
        
        col_stat1, col_stat2, col_stat3 = st.columns(3)
        
        with col_stat1:
            st.metric("Mean", f"{np.mean(selected_expression):.3f}", help="Average expression across selected cells")
        with col_stat2:
            st.metric("Max", f"{np.max(selected_expression):.3f}", help="Maximum expression value in selection")
        with col_stat3:
            pct_expressing = (np.sum(selected_expression > 0) / len(selected_expression) * 100)
            st.metric("% Expr", f"{pct_expressing:.1f}%", help="Percentage of selected cells with non-zero expression")
