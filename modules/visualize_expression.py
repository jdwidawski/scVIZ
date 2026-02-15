"""
Module 2: Visualize Gene Expression

Provides interactive visualizations for gene expression analysis including
UMAP plots, violin plots, and heatmaps.
"""
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from anndata import AnnData
import scanpy as sc
from typing import List, Optional


@st.dialog("About Visualize Gene Expression", width="large")
def show_expression_info_modal():
    """Display information about the Visualize Gene Expression module."""
    st.markdown("""
    ## What is Visualize Gene Expression?
    
    This module provides multiple visualization types to explore gene expression 
    patterns across your single-cell dataset.
    
    ### Visualization Types
    
    #### 🗺️ Dimensionality Reduction (UMAP/t-SNE)
    - Project cells in 2D space colored by gene expression
    - Compare expression levels across cell populations
    - Identify marker gene patterns
    
    #### 🎻 Violin Plots
    - Show expression distribution across groups
    - Compare expression between cell types or conditions
    - Visualize statistical differences
    
    #### 🔥 Heatmap
    - Display multiple genes simultaneously
    - Group cells by metadata categories
    - Identify co-expression patterns
    
    #### 🔴 Dot Plot
    - Show both expression level and percent expressing
    - Compact view for many genes across groups
    - Commonly used for marker gene visualization
    
    ### Gene Search
    
    - Search by gene name, symbol, or ID
    - Supports partial matching for discovery
    - Shows alternative symbols when available
    
    ### Tips for Analysis
    
    - **Marker genes**: Use dot plots to visualize known markers
    - **Expression gradients**: UMAP shows spatial expression patterns
    - **Group comparisons**: Violin plots reveal distribution differences
    - **Co-expression**: Heatmaps identify genes with similar patterns
    """)
    
    st.info("💡 **Tip**: Start with known marker genes to validate cell type annotations, then explore novel genes.")


def render(adata: AnnData) -> None:
    """
    Render the Visualize Gene Expression module.
    
    Parameters:
        adata (AnnData): The loaded AnnData object to visualize.
    
    This module provides:
    - UMAP/t-SNE scatter plots colored by metadata or gene expression
    - Violin plots for gene expression across groups
    - Heatmaps for multiple genes across cell types
    """
    # Title with info button
    col_title, col_info = st.columns([6, 1])
    with col_title:
        st.title("🔬 Visualize Gene Expression")
    with col_info:
        st.markdown("<div style='margin-top: 12px;'></div>", unsafe_allow_html=True)
        if st.button("ℹ️ Help", key="expression_help_btn", use_container_width=True):
            show_expression_info_modal()
    
    st.markdown("---")
    
    # Check for dimensional reductions
    has_umap = 'X_umap' in adata.obsm.keys()
    has_tsne = 'X_tsne' in adata.obsm.keys()
    has_pca = 'X_pca' in adata.obsm.keys()
    
    if not (has_umap or has_tsne or has_pca):
        st.warning("No dimensionality reduction found. Computing UMAP...")
        with st.spinner("Computing PCA and UMAP..."):
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
        st.success("UMAP computed!")
        has_umap = True
    
    # Visualization type selector
    viz_type = st.selectbox(
        "Select Visualization Type:",
        ["Dimensionality Reduction (UMAP/t-SNE)", "Violin Plots", "Heatmap", "Dot Plot"]
    )
    
    st.markdown("---")
    
    if viz_type == "Dimensionality Reduction (UMAP/t-SNE)":
        render_scatter_plot(adata, has_umap, has_tsne, has_pca)
    
    elif viz_type == "Violin Plots":
        render_violin_plot(adata)
    
    elif viz_type == "Heatmap":
        render_heatmap(adata)
    
    elif viz_type == "Dot Plot":
        render_dot_plot(adata)


def render_scatter_plot(adata: AnnData, has_umap: bool, has_tsne: bool, has_pca: bool) -> None:
    """
    Render dimensionality reduction scatter plots.
    
    Parameters:
        adata (AnnData): The AnnData object.
        has_umap (bool): Whether UMAP coordinates are available.
        has_tsne (bool): Whether t-SNE coordinates are available.
        has_pca (bool): Whether PCA coordinates are available.
    """
    st.header("Dimensionality Reduction")
    
    # Select reduction type
    reduction_options = []
    if has_umap:
        reduction_options.append("UMAP")
    if has_tsne:
        reduction_options.append("t-SNE")
    if has_pca:
        reduction_options.append("PCA")
    
    reduction_type = st.selectbox("Select reduction:", reduction_options)
    
    # Map to obsm key
    reduction_map = {"UMAP": "X_umap", "t-SNE": "X_tsne", "PCA": "X_pca"}
    coords_key = reduction_map[reduction_type]
    
    # Color by options
    color_by = st.radio("Color by:", ["Metadata", "Gene Expression"])
    
    col1, col2 = st.columns([2, 1])
    
    if color_by == "Metadata":
        categorical_cols = [col for col in adata.obs.columns 
                           if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category']
        numerical_cols = [col for col in adata.obs.columns 
                         if adata.obs[col].dtype in ['int64', 'float64']]
        
        all_cols = categorical_cols + numerical_cols
        
        if all_cols:
            selected_col = col2.selectbox("Select metadata column:", all_cols)
            
            # Create scatter plot
            coords = adata.obsm[coords_key]
            df_plot = pd.DataFrame({
                f'{reduction_type}1': coords[:, 0],
                f'{reduction_type}2': coords[:, 1],
                'color': adata.obs[selected_col].values
            })
            
            with col1:
                fig = px.scatter(
                    df_plot,
                    x=f'{reduction_type}1',
                    y=f'{reduction_type}2',
                    color='color',
                    title=f'{reduction_type} colored by {selected_col}',
                    width=700,
                    height=600
                )
                fig.update_traces(marker=dict(size=3))
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No metadata available for coloring")
    
    else:  # Gene Expression
        # Gene selection
        gene_name = col2.text_input("Enter gene name:", value=adata.var_names[0] if len(adata.var_names) > 0 else "")
        
        # Optional: Add gene search
        with col2.expander("Search for genes"):
            search_term = st.text_input("Search:", "")
            if search_term:
                matching_genes = [g for g in adata.var_names if search_term.lower() in g.lower()]
                if matching_genes:
                    st.write(f"Found {len(matching_genes)} genes:")
                    selected_from_search = st.selectbox("Select gene:", matching_genes[:50])
                    if st.button("Use this gene"):
                        gene_name = selected_from_search
                        st.rerun()
        
        if gene_name in adata.var_names:
            gene_idx = adata.var_names.get_loc(gene_name)
            
            if hasattr(adata.X, 'toarray'):
                expression = adata.X[:, gene_idx].toarray().flatten()
            else:
                expression = adata.X[:, gene_idx].flatten()
            
            coords = adata.obsm[coords_key]
            df_plot = pd.DataFrame({
                f'{reduction_type}1': coords[:, 0],
                f'{reduction_type}2': coords[:, 1],
                'expression': expression
            })
            
            with col1:
                fig = px.scatter(
                    df_plot,
                    x=f'{reduction_type}1',
                    y=f'{reduction_type}2',
                    color='expression',
                    title=f'{reduction_type} - {gene_name} expression',
                    color_continuous_scale='Viridis',
                    width=700,
                    height=600
                )
                fig.update_traces(marker=dict(size=3))
                st.plotly_chart(fig, use_container_width=True)
            
            # Expression statistics
            with col2:
                st.markdown("**Expression Statistics:**")
                st.metric("Mean", f"{np.mean(expression):.3f}")
                st.metric("Max", f"{np.max(expression):.3f}")
                st.metric("% Expressing", f"{(np.sum(expression > 0) / len(expression) * 100):.1f}%")
        
        elif gene_name:
            col1.error(f"Gene '{gene_name}' not found in dataset")


def render_violin_plot(adata: AnnData) -> None:
    """
    Render violin plots for gene expression across groups.
    
    Parameters:
        adata (AnnData): The AnnData object.
    """
    st.header("Violin Plots")
    
    col1, col2 = st.columns([3, 1])
    
    with col2:
        # Select genes
        genes_input = st.text_area(
            "Enter gene names (one per line):",
            value="\n".join(adata.var_names[:3].tolist()) if len(adata.var_names) >= 3 else ""
        )
        genes = [g.strip() for g in genes_input.split('\n') if g.strip()]
        
        # Select grouping variable
        categorical_cols = [col for col in adata.obs.columns 
                           if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category']
        
        if categorical_cols:
            groupby = st.selectbox("Group by:", categorical_cols)
        else:
            st.warning("No categorical metadata available for grouping")
            groupby = None
    
    if genes and groupby:
        valid_genes = [g for g in genes if g in adata.var_names]
        
        if valid_genes:
            with col1:
                # Create violin plot using scanpy
                fig, axes = plt.subplots(len(valid_genes), 1, figsize=(12, 4 * len(valid_genes)))
                
                if len(valid_genes) == 1:
                    axes = [axes]
                
                for idx, gene in enumerate(valid_genes):
                    sc.pl.violin(adata, gene, groupby=groupby, ax=axes[idx], show=False)
                
                plt.tight_layout()
                st.pyplot(fig)
                plt.close()
        else:
            col1.error("None of the specified genes found in dataset")


def render_heatmap(adata: AnnData) -> None:
    """
    Render heatmap for multiple genes across cell types.
    
    Parameters:
        adata (AnnData): The AnnData object.
    """
    st.header("Expression Heatmap")
    
    col1, col2 = st.columns([3, 1])
    
    with col2:
        # Select genes
        genes_input = st.text_area(
            "Enter gene names (one per line):",
            value="\n".join(adata.var_names[:10].tolist()) if len(adata.var_names) >= 10 else ""
        )
        genes = [g.strip() for g in genes_input.split('\n') if g.strip()]
        
        # Select grouping variable
        categorical_cols = [col for col in adata.obs.columns 
                           if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category']
        
        if categorical_cols:
            groupby = st.selectbox("Group by:", categorical_cols)
        else:
            groupby = None
    
    if genes and groupby:
        valid_genes = [g for g in genes if g in adata.var_names]
        
        if valid_genes:
            with col1:
                # Create heatmap using scanpy
                fig = sc.pl.heatmap(
                    adata,
                    valid_genes,
                    groupby=groupby,
                    show=False,
                    return_fig=True,
                    figsize=(10, len(valid_genes) * 0.4)
                )
                st.pyplot(fig)
                plt.close()
        else:
            col1.error("None of the specified genes found in dataset")


def render_dot_plot(adata: AnnData) -> None:
    """
    Render dot plot for gene expression across groups.
    
    Parameters:
        adata (AnnData): The AnnData object.
    """
    st.header("Dot Plot")
    
    col1, col2 = st.columns([3, 1])
    
    with col2:
        # Select genes
        genes_input = st.text_area(
            "Enter gene names (one per line):",
            value="\n".join(adata.var_names[:10].tolist()) if len(adata.var_names) >= 10 else ""
        )
        genes = [g.strip() for g in genes_input.split('\n') if g.strip()]
        
        # Select grouping variable
        categorical_cols = [col for col in adata.obs.columns 
                           if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category']
        
        if categorical_cols:
            groupby = st.selectbox("Group by:", categorical_cols)
        else:
            groupby = None
    
    if genes and groupby:
        valid_genes = [g for g in genes if g in adata.var_names]
        
        if valid_genes:
            with col1:
                # Create dot plot using scanpy
                fig = sc.pl.dotplot(
                    adata,
                    valid_genes,
                    groupby=groupby,
                    show=False,
                    return_fig=True
                )
                st.pyplot(fig)
                plt.close()
        else:
            col1.error("None of the specified genes found in dataset")
