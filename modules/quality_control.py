"""
Module: Quality Control

Provides interactive quality control visualizations for single-cell RNA-seq data.
Includes metrics like total counts, total genes, mitochondrial/ribosomal percentages,
and cell cycle scoring.
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from anndata import AnnData
import scanpy as sc
from typing import List, Optional, Tuple


# Cell cycle genes for humans (from Tirosh et al. 2016, used by Seurat/Scanpy)
# Gene symbols
S_GENES = [
    'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2',
    'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2',
    'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7',
    'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1',
    'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B',
    'BRIP1', 'E2F8'
]

G2M_GENES = [
    'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80',
    'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A',
    'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E',
    'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK',
    'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8',
    'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5',
    'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA'
]

# Ensembl IDs for S phase genes (human)
S_GENES_ENSEMBL = [
    'ENSG00000100297', 'ENSG00000132646', 'ENSG00000176890', 'ENSG00000168496',
    'ENSG00000073111', 'ENSG00000104738', 'ENSG00000167325', 'ENSG00000076248',
    'ENSG00000131153', 'ENSG00000076003', 'ENSG00000136982', 'ENSG00000143476',
    'ENSG00000198056', 'ENSG00000276043', 'ENSG00000101868', 'ENSG00000119969',
    'ENSG00000049541', 'ENSG00000117748', 'ENSG00000132780', 'ENSG00000111247',
    'ENSG00000112312', 'ENSG00000171848', 'ENSG00000129173', 'ENSG00000175305',
    'ENSG00000094804', 'ENSG00000077514', 'ENSG00000095002', 'ENSG00000092470',
    'ENSG00000051180', 'ENSG00000136492', 'ENSG00000093009', 'ENSG00000162607',
    'ENSG00000163950', 'ENSG00000118412', 'ENSG00000075131', 'ENSG00000197299',
    'ENSG00000174371', 'ENSG00000144354', 'ENSG00000156802', 'ENSG00000151725',
    'ENSG00000012963', 'ENSG00000092853', 'ENSG00000159259'
]

# Ensembl IDs for G2M phase genes (human)
G2M_GENES_ENSEMBL = [
    'ENSG00000164104', 'ENSG00000170312', 'ENSG00000137804', 'ENSG00000175063',
    'ENSG00000089685', 'ENSG00000088325', 'ENSG00000131747', 'ENSG00000080986',
    'ENSG00000123485', 'ENSG00000143228', 'ENSG00000117399', 'ENSG00000148773',
    'ENSG00000120802', 'ENSG00000117650', 'ENSG00000013810', 'ENSG00000178999',
    'ENSG00000113810', 'ENSG00000157456', 'ENSG00000169679', 'ENSG00000138182',
    'ENSG00000175216', 'ENSG00000169607', 'ENSG00000138778', 'ENSG00000100401',
    'ENSG00000188229', 'ENSG00000075218', 'ENSG00000087586', 'ENSG00000102974',
    'ENSG00000111665', 'ENSG00000184661', 'ENSG00000117724', 'ENSG00000112742',
    'ENSG00000158402', 'ENSG00000126787', 'ENSG00000011426', 'ENSG00000092140',
    'ENSG00000137807', 'ENSG00000134222', 'ENSG00000134690', 'ENSG00000114346',
    'ENSG00000072571', 'ENSG00000173207', 'ENSG00000087586', 'ENSG00000115163',
    'ENSG00000189159', 'ENSG00000143815', 'ENSG00000094916', 'ENSG00000138160',
    'ENSG00000010292', 'ENSG00000129195', 'ENSG00000139354', 'ENSG00000142945',
    'ENSG00000136108', 'ENSG00000143401'
]

# Mitochondrial gene Ensembl IDs (human MT genes)
MT_GENES_ENSEMBL_PREFIXES = ['ENSG00000198888', 'ENSG00000198763', 'ENSG00000198804',
                              'ENSG00000198712', 'ENSG00000198899', 'ENSG00000198938',
                              'ENSG00000198840', 'ENSG00000212907', 'ENSG00000198886',
                              'ENSG00000198786', 'ENSG00000198695', 'ENSG00000198727',
                              'ENSG00000198747']

# Available colormaps for continuous data
CONTINUOUS_COLORMAPS = [
    'Reds', 'Viridis', 'Plasma', 'Inferno', 'Magma', 'Cividis',
    'Blues', 'Greens', 'Oranges', 'Purples', 'YlOrRd', 'YlGnBu', 'RdBu'
]

# Available color palettes for categorical data (violin plots)
CATEGORICAL_PALETTES = {
    'Plotly': px.colors.qualitative.Plotly,
    'D3': px.colors.qualitative.D3,
    'Set1': px.colors.qualitative.Set1,
    'Set2': px.colors.qualitative.Set2,
    'Set3': px.colors.qualitative.Set3,
    'Pastel1': px.colors.qualitative.Pastel1,
    'Pastel2': px.colors.qualitative.Pastel2,
    'Dark2': px.colors.qualitative.Dark2,
    'Bold': px.colors.qualitative.Bold,
    'Vivid': px.colors.qualitative.Vivid,
    'Safe': px.colors.qualitative.Safe,
    'Alphabet': px.colors.qualitative.Alphabet,
}


def get_categorical_columns(adata: AnnData, max_categories: int = 50) -> List[str]:
    """Get list of categorical columns from obs with limited unique values."""
    categorical_cols = []
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category':
            if adata.obs[col].nunique() <= max_categories:
                categorical_cols.append(col)
    return categorical_cols


def find_genes_in_dataset(adata: AnnData, gene_symbols: List[str], 
                          ensembl_ids: List[str]) -> List[str]:
    """
    Find genes in dataset by matching either gene symbols or Ensembl IDs.
    
    Checks both var_names and common var columns for gene symbols.
    
    Parameters:
        adata: AnnData object
        gene_symbols: List of gene symbols to search for
        ensembl_ids: List of Ensembl IDs to search for
        
    Returns:
        List of gene names as they appear in adata.var_names
    """
    found_genes = []
    var_names_upper = adata.var_names.str.upper()
    
    # Check for gene symbol column in var
    symbol_col = None
    for col in ['gene_symbol', 'gene_name', 'symbol', 'gene', 'feature_name']:
        if col in adata.var.columns:
            symbol_col = col
            break
    
    # Try to match gene symbols directly in var_names
    for symbol in gene_symbols:
        # Direct match
        if symbol in adata.var_names:
            found_genes.append(symbol)
            continue
        
        # Case-insensitive match
        matches = adata.var_names[var_names_upper == symbol.upper()]
        if len(matches) > 0:
            found_genes.append(matches[0])
            continue
        
        # Check symbol column
        if symbol_col is not None:
            col_upper = adata.var[symbol_col].astype(str).str.upper()
            matches = adata.var_names[col_upper == symbol.upper()]
            if len(matches) > 0:
                found_genes.append(matches[0])
                continue
    
    # Try to match Ensembl IDs
    for ensembl in ensembl_ids:
        if ensembl in found_genes:
            continue
            
        # Direct match
        if ensembl in adata.var_names:
            found_genes.append(ensembl)
            continue
        
        # Partial match (e.g., ENSG00000123456.1 -> ENSG00000123456)
        base_id = ensembl.split('.')[0]
        matches = adata.var_names[adata.var_names.str.startswith(base_id)]
        if len(matches) > 0:
            if matches[0] not in found_genes:
                found_genes.append(matches[0])
    
    return list(set(found_genes))


def identify_mt_genes(adata: AnnData) -> pd.Series:
    """
    Identify mitochondrial genes by symbol or Ensembl ID.
    
    Returns:
        Boolean Series indicating mitochondrial genes
    """
    # Check by gene symbol prefix (MT-)
    mt_mask = adata.var_names.str.upper().str.startswith('MT-')
    
    # Also check for Ensembl IDs
    for ensembl in MT_GENES_ENSEMBL_PREFIXES:
        mt_mask = mt_mask | adata.var_names.str.startswith(ensembl)
    
    # Check gene symbol column if available
    for col in ['gene_symbol', 'gene_name', 'symbol', 'gene', 'feature_name']:
        if col in adata.var.columns:
            col_mask = adata.var[col].astype(str).str.upper().str.startswith('MT-')
            mt_mask = mt_mask | col_mask
            break
    
    return mt_mask


def identify_ribo_genes(adata: AnnData) -> pd.Series:
    """
    Identify ribosomal genes by symbol or Ensembl ID.
    
    Returns:
        Boolean Series indicating ribosomal genes
    """
    # Check by gene symbol prefix (RPS, RPL)
    ribo_mask = (
        adata.var_names.str.upper().str.startswith('RPS') | 
        adata.var_names.str.upper().str.startswith('RPL')
    )
    
    # Check gene symbol column if available
    for col in ['gene_symbol', 'gene_name', 'symbol', 'gene', 'feature_name']:
        if col in adata.var.columns:
            col_upper = adata.var[col].astype(str).str.upper()
            col_mask = col_upper.str.startswith('RPS') | col_upper.str.startswith('RPL')
            ribo_mask = ribo_mask | col_mask
            break
    
    return ribo_mask


def compute_qc_metrics(adata: AnnData) -> AnnData:
    """
    Compute quality control metrics for the dataset.
    
    Computes:
    - total_counts: Total UMI counts per cell
    - n_genes_by_counts: Number of genes with positive counts
    - pct_counts_mt: Percentage of mitochondrial gene counts
    - pct_counts_ribo: Percentage of ribosomal gene counts
    - S_score, G2M_score, phase: Cell cycle scores
    
    Parameters:
        adata: AnnData object
        
    Returns:
        AnnData with QC metrics added to obs
    """
    # Identify mitochondrial genes
    if 'mt' not in adata.var.columns:
        adata.var['mt'] = identify_mt_genes(adata)
    
    # Identify ribosomal genes
    if 'ribo' not in adata.var.columns:
        adata.var['ribo'] = identify_ribo_genes(adata)
    
    # Calculate QC metrics if not already present
    if 'total_counts' not in adata.obs.columns or 'pct_counts_mt' not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(
            adata, 
            qc_vars=['mt', 'ribo'],
            percent_top=None,
            log1p=False,
            inplace=True
        )
    
    return adata


def compute_cell_cycle_scores(adata: AnnData) -> AnnData:
    """
    Compute cell cycle phase scores using S and G2M marker genes.
    
    Supports both gene symbols and Ensembl IDs.
    
    Parameters:
        adata: AnnData object
        
    Returns:
        AnnData with cell cycle scores added to obs
    """
    # Check if already computed
    if 'S_score' in adata.obs.columns and 'G2M_score' in adata.obs.columns:
        return adata
    
    # Find S phase genes (try both symbols and Ensembl IDs)
    s_genes_found = find_genes_in_dataset(adata, S_GENES, S_GENES_ENSEMBL)
    g2m_genes_found = find_genes_in_dataset(adata, G2M_GENES, G2M_GENES_ENSEMBL)
    
    if len(s_genes_found) < 5 or len(g2m_genes_found) < 5:
        # Set default values if not enough genes found
        adata.obs['S_score'] = 0.0
        adata.obs['G2M_score'] = 0.0
        adata.obs['phase'] = 'Unknown'
        return adata
    
    # Score cell cycle using scanpy
    try:
        sc.tl.score_genes_cell_cycle(
            adata,
            s_genes=s_genes_found,
            g2m_genes=g2m_genes_found
        )
    except Exception as e:
        adata.obs['S_score'] = 0.0
        adata.obs['G2M_score'] = 0.0
        adata.obs['phase'] = 'Unknown'
    
    return adata


@st.dialog("About Quality Control", width="large")
def show_qc_info_modal():
    """Display information about the Quality Control module."""
    st.markdown("""
    ## What is Quality Control?
    
    Quality control (QC) is an essential step in single-cell RNA-seq analysis 
    to identify and filter out low-quality cells that may introduce noise.
    
    ### QC Metrics
    
    | Metric | Description | Typical Threshold |
    |--------|-------------|-------------------|
    | **Total Counts** | Total UMI counts per cell | > 500-1000 |
    | **Genes Detected** | Number of genes with counts | > 200-500 |
    | **% Mitochondrial** | Percentage of MT gene counts | < 10-20% |
    | **% Ribosomal** | Percentage of ribosomal counts | Dataset-dependent |
    
    ### Interpreting Metrics
    
    #### Low Counts / Few Genes
    - May indicate dying cells or poor capture
    - Filter cells with very low counts/genes
    
    #### High Mitochondrial %
    - Often indicates dying or stressed cells
    - Mitochondrial genes are retained when cytoplasmic mRNA leaks out
    - Typical cutoff: 10-20% depending on tissue
    
    #### Ribosomal Content
    - High in actively translating cells
    - Dataset/cell-type dependent
    - Consider in context of biology
    
    ### Cell Cycle Scoring
    
    Cells are assigned to cell cycle phases based on expression of marker genes:
    - **G1**: Cells in gap/growth phase 1
    - **S**: Cells actively replicating DNA
    - **G2M**: Cells in gap 2 or mitosis
    
    Uses marker genes from [Tirosh et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4944528/).
    
    ### Gene Matching
    
    This module supports both **gene symbols** (e.g., MT-CO1, RPS6) and 
    **Ensembl IDs** (e.g., ENSG00000198888) for identifying mitochondrial, 
    ribosomal, and cell cycle genes.
    """)
    
    st.info("💡 **Tip**: Visualize metrics on UMAP to identify spatial patterns of quality issues.")


def render(adata: AnnData) -> None:
    """
    Render the Quality Control module.
    
    Parameters:
        adata (AnnData): The loaded AnnData object to analyze.
    """
    # Title with info button
    col_title, col_info = st.columns([6, 1])
    with col_title:
        st.title("🔍 Quality Control")
    with col_info:
        st.markdown("<div style='margin-top: 12px;'></div>", unsafe_allow_html=True)
        if st.button("ℹ️ Help", key="qc_help_btn", use_container_width=True):
            show_qc_info_modal()
    
    st.markdown("Explore quality control metrics to assess dataset quality.")
    st.markdown("---")
    
    # Check for UMAP
    has_umap = 'X_umap' in adata.obsm.keys()
    
    if not has_umap:
        st.warning("⚠️ No UMAP embedding found. UMAP visualization requires computing UMAP first.")
        if st.button("🔄 Compute UMAP", type="primary"):
            with st.spinner("Computing PCA and UMAP..."):
                if 'X_pca' not in adata.obsm.keys():
                    sc.pp.pca(adata)
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)
            st.success("✅ UMAP computed successfully!")
            st.rerun()
        return
    
    # Compute QC metrics and cell cycle scores automatically
    with st.spinner("Computing QC metrics..."):
        adata = compute_qc_metrics(adata)
        adata = compute_cell_cycle_scores(adata)
    
    # Available QC metrics
    qc_metrics = []
    qc_metric_names = {
        'total_counts': 'Total Counts',
        'n_genes_by_counts': 'Genes Detected',
        'pct_counts_mt': '% Mitochondrial',
        'pct_counts_ribo': '% Ribosomal',
        'S_score': 'S Phase Score',
        'G2M_score': 'G2/M Phase Score'
    }
    
    for metric, name in qc_metric_names.items():
        if metric in adata.obs.columns:
            qc_metrics.append((metric, name))
    
    if len(qc_metrics) == 0:
        st.error("❌ No QC metrics available. Please ensure the dataset has count data.")
        return
    
    # Metric selection
    st.subheader("🎯 Select Metric")
    
    metric_options = [name for _, name in qc_metrics]
    metric_lookup = {name: key for key, name in qc_metrics}
    
    selected_metric_name = st.selectbox(
        "Choose a QC metric to visualize:",
        options=metric_options,
        key="qc_metric_select"
    )
    
    selected_metric = metric_lookup[selected_metric_name]
    
    # Show summary statistics
    metric_values = adata.obs[selected_metric].values
    col_s1, col_s2, col_s3, col_s4 = st.columns(4)
    
    with col_s1:
        st.metric("Mean", f"{np.mean(metric_values):.2f}")
    with col_s2:
        st.metric("Median", f"{np.median(metric_values):.2f}")
    with col_s3:
        st.metric("Min", f"{np.min(metric_values):.2f}")
    with col_s4:
        st.metric("Max", f"{np.max(metric_values):.2f}")
    
    st.markdown("---")
    
    # UMAP visualization section
    st.subheader("🗺️ UMAP Visualization")
    st.markdown(f"**UMAP colored by {selected_metric_name}**")
    st.caption("Use lasso or box select to choose cells for violin plot comparison below.")
    
    # UMAP Figure settings expander
    with st.expander("⚙️ UMAP Figure Settings", expanded=False):
        col_s1, col_s2, col_s3 = st.columns(3)
        
        with col_s1:
            width = st.slider("Width (px)", 300, 1200, 700, 50, key="qc_umap_width")
            height = st.slider("Height (px)", 300, 1200, 600, 50, key="qc_umap_height")
        
        with col_s2:
            dot_size = st.slider("Dot size", 1, 20, 4, 1, key="qc_umap_dotsize")
            colormap = st.selectbox("Colormap", CONTINUOUS_COLORMAPS, index=0, key="qc_umap_cmap")
        
        with col_s3:
            # Color range mode - default to Percentile
            color_mode = st.radio(
                "Color range mode:",
                ["Percentile", "Auto", "Absolute"],
                index=0,
                key="qc_color_mode",
                help="Percentile: use percentiles. Auto: use data min/max. Absolute: specify exact values."
            )
            
            if color_mode == "Percentile":
                vmin_pct = st.number_input("Min percentile", 0.0, 100.0, 1.0, 0.5, key="qc_vmin_pct",
                                          help="e.g., 1 for 1st percentile")
                vmax_pct = st.number_input("Max percentile", 0.0, 100.0, 99.5, 0.5, key="qc_vmax_pct",
                                          help="e.g., 99.5 for 99.5th percentile")
            elif color_mode == "Absolute":
                vmin_abs = st.number_input("Color min", value=float(np.min(metric_values)), key="qc_vmin_abs")
                vmax_abs = st.number_input("Color max", value=float(np.max(metric_values)), key="qc_vmax_abs")
    
    # Compute color range
    if color_mode == "Percentile":
        color_min = np.percentile(metric_values, vmin_pct)
        color_max = np.percentile(metric_values, vmax_pct)
    elif color_mode == "Absolute":
        color_min = vmin_abs
        color_max = vmax_abs
    else:  # Auto
        color_min = np.min(metric_values)
        color_max = np.max(metric_values)
    
    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']
    
    # Create dataframe with original indices for selection tracking
    df_plot = pd.DataFrame({
        'UMAP1': umap_coords[:, 0],
        'UMAP2': umap_coords[:, 1],
        selected_metric_name: metric_values,
        'original_idx': range(len(umap_coords))
    })
    
    # Sort by metric value so cells with higher values are plotted on top
    df_plot = df_plot.sort_values(selected_metric_name, ascending=True).reset_index(drop=True)
    
    # Create figure
    fig = px.scatter(
        df_plot,
        x='UMAP1',
        y='UMAP2',
        color=selected_metric_name,
        color_continuous_scale=colormap.lower(),
        range_color=[color_min, color_max],
        hover_data={selected_metric_name: ':.2f', 'original_idx': True},
        title=f'{selected_metric_name} on UMAP',
        width=width,
        height=height,
        custom_data=['original_idx']
    )
    
    fig.update_traces(marker=dict(size=dot_size))
    
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        dragmode='lasso',
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
            title=selected_metric_name,
            thickness=15
        )
    )
    
    # Display chart with selection
    event = st.plotly_chart(fig, use_container_width=False, on_select="rerun", key="qc_umap_select")
    
    # Process selection - extract original indices from custom_data
    selected_indices = []
    if event and event.selection and event.selection.point_indices:
        # Get the original indices from the sorted dataframe
        for pt_idx in event.selection.point_indices:
            original_idx = int(df_plot.iloc[pt_idx]['original_idx'])
            selected_indices.append(original_idx)
    
    st.markdown("---")
    
    # Violin plot section
    st.subheader("🎻 Violin Plot")
    
    # Grouping options
    grouping_mode = st.radio(
        "Group by:",
        options=["Full dataset", "Metadata category", "UMAP selection (above)"],
        key="qc_violin_grouping",
        horizontal=True
    )
    
    # Category selection (shown inline only when needed)
    selected_category = None
    if grouping_mode == "Metadata category":
        categorical_cols = get_categorical_columns(adata)
        if len(categorical_cols) == 0:
            st.warning("⚠️ No categorical columns found in metadata.")
        else:
            selected_category = st.selectbox(
                "Select category:",
                options=categorical_cols,
                key="qc_violin_category"
            )
    
    # Violin figure settings
    with st.expander("⚙️ Violin Figure Settings", expanded=False):
        col_v1, col_v2, col_v3, col_v4 = st.columns(4)
        
        with col_v1:
            violin_width = st.slider("Width (px)", 300, 1200, 700, 50, key="qc_violin_width")
        with col_v2:
            violin_height = st.slider("Height (px)", 200, 800, 400, 50, key="qc_violin_height")
        with col_v3:
            violin_palette = st.selectbox("Color palette", list(CATEGORICAL_PALETTES.keys()), 
                                         index=0, key="qc_violin_palette")
        with col_v4:
            use_boxplot = st.checkbox("Show box plot", value=False, key="qc_violin_box",
                                     help="Overlay a box plot on the violin")
    
    # Handle different grouping modes - violin plot always below
    if grouping_mode == "Full dataset":
        render_violin_full_dataset(metric_values, selected_metric_name, violin_width, violin_height, 
                                   violin_palette, use_boxplot)
    
    elif grouping_mode == "Metadata category":
        if selected_category is not None:
            render_violin_by_category(adata, metric_values, selected_metric_name, selected_category,
                                     violin_width, violin_height, violin_palette, use_boxplot)
    
    else:  # UMAP selection
        if len(selected_indices) == 0:
            st.info("👆 Use lasso or box select on the UMAP above to select cells for comparison")
        else:
            st.success(f"✅ Selected **{len(selected_indices):,}** cells from UMAP")
            render_violin_by_selection(metric_values, selected_metric_name, selected_indices,
                                       violin_width, violin_height, violin_palette, use_boxplot)


def render_violin_full_dataset(metric_values: np.ndarray, metric_name: str,
                               width: int = 700, height: int = 400,
                               palette: str = 'Plotly', use_boxplot: bool = True) -> None:
    """Render violin plot for full dataset."""
    df_plot = pd.DataFrame({
        metric_name: metric_values,
        'Group': 'All Cells'
    })
    
    # Get colors from categorical palette
    colors = CATEGORICAL_PALETTES.get(palette, px.colors.qualitative.Plotly)
    
    fig = px.violin(
        df_plot,
        y=metric_name,
        x='Group',
        box=use_boxplot,
        points='outliers',
        title=f'{metric_name} Distribution',
        color_discrete_sequence=[colors[0]]
    )
    
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=width,
        height=height,
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=False)
    
    # Show distribution stats as table
    st.markdown("**Distribution Statistics:**")
    stats_data = {
        'Statistic': ['Count', 'Mean', 'Std', 'Min', 'Q1 (25%)', 'Median (50%)', 'Q3 (75%)', 'Max'],
        'Value': [
            f"{len(metric_values):,}",
            f"{np.mean(metric_values):.2f}",
            f"{np.std(metric_values):.2f}",
            f"{np.min(metric_values):.2f}",
            f"{np.percentile(metric_values, 25):.2f}",
            f"{np.percentile(metric_values, 50):.2f}",
            f"{np.percentile(metric_values, 75):.2f}",
            f"{np.max(metric_values):.2f}"
        ]
    }
    st.dataframe(pd.DataFrame(stats_data), use_container_width=True, hide_index=True)


def render_violin_by_category(adata: AnnData, metric_values: np.ndarray, 
                              metric_name: str, category: str,
                              width: int = 700, height: int = 400,
                              palette: str = 'Plotly', use_boxplot: bool = True) -> None:
    """Render violin plot grouped by metadata category."""
    df_plot = pd.DataFrame({
        metric_name: metric_values,
        category: adata.obs[category].values
    })
    
    # Order by median value
    order = df_plot.groupby(category)[metric_name].median().sort_values(ascending=False).index.tolist()
    
    # Get colors from categorical palette
    colors = CATEGORICAL_PALETTES.get(palette, px.colors.qualitative.Plotly)
    
    fig = px.violin(
        df_plot,
        y=metric_name,
        x=category,
        box=use_boxplot,
        points='outliers',
        title=f'{metric_name} by {category}',
        category_orders={category: order},
        color=category,
        color_discrete_sequence=colors
    )
    
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=width,
        height=height,
        xaxis_tickangle=-45,
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=False)
    
    # Show category statistics
    st.markdown("**Category Statistics:**")
    stats_df = df_plot.groupby(category)[metric_name].agg(['mean', 'median', 'std', 'count'])
    stats_df = stats_df.round(2).sort_values('median', ascending=False)
    st.dataframe(stats_df, use_container_width=True, height=200)


def render_violin_by_selection(metric_values: np.ndarray, metric_name: str, 
                               selected_indices: List[int],
                               width: int = 700, height: int = 400,
                               palette: str = 'Plotly', use_boxplot: bool = True) -> None:
    """Render violin plot comparing selected vs rest cells."""
    # Create comparison groups
    group_labels = np.array(['Rest'] * len(metric_values))
    for idx in selected_indices:
        group_labels[idx] = 'Selected'
    
    df_plot = pd.DataFrame({
        metric_name: metric_values,
        'Group': group_labels
    })
    
    # Get colors from categorical palette
    colors = CATEGORICAL_PALETTES.get(palette, px.colors.qualitative.Plotly)
    
    fig = px.violin(
        df_plot,
        y=metric_name,
        x='Group',
        color='Group',
        box=use_boxplot,
        points='outliers',
        title=f'{metric_name}: Selected vs Rest',
        color_discrete_map={'Rest': colors[0], 'Selected': colors[1]}
    )
    
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=width,
        height=height,
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=False)
    
    # Show comparison statistics as table
    st.markdown("**Comparison Statistics:**")
    
    selected_values = metric_values[np.array(selected_indices)]
    rest_mask = np.ones(len(metric_values), dtype=bool)
    rest_mask[selected_indices] = False
    rest_values = metric_values[rest_mask]
    
    stats_data = {
        'Statistic': ['Count', 'Mean', 'Median', 'Std', 'Min', 'Max'],
        'Selected': [
            f"{len(selected_values):,}",
            f"{selected_values.mean():.2f}",
            f"{np.median(selected_values):.2f}",
            f"{selected_values.std():.2f}",
            f"{np.min(selected_values):.2f}",
            f"{np.max(selected_values):.2f}"
        ],
        'Rest': [
            f"{len(rest_values):,}",
            f"{rest_values.mean():.2f}",
            f"{np.median(rest_values):.2f}",
            f"{rest_values.std():.2f}",
            f"{np.min(rest_values):.2f}",
            f"{np.max(rest_values):.2f}"
        ]
    }
    st.dataframe(pd.DataFrame(stats_data), use_container_width=True, hide_index=True)
