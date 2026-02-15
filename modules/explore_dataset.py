"""
Module 1: Explore Dataset Contents

Provides interactive visualizations and summaries of dataset metadata,
exploring the AnnData structure visually.
"""
import streamlit as st
import pandas as pd
import numpy as np
import scipy
import random
import json
from anndata import AnnData
from typing import List


MAX_N_CATEGORIES = 50  # Maximum number of unique categories to display

# Tooltips for each AnnData component (from official documentation)
TOOLTIPS = {
    'var': 'One-dimensional annotation of variables/features (DataFrame). Contains gene metadata.',
    'obs': 'One-dimensional annotation of observations (DataFrame). Contains cell/sample metadata.',
    'X': 'Data matrix of shape n_obs × n_vars. The main data matrix (e.g., gene expression counts).',
    'obsm': 'Multi-dimensional annotation of observations. Stores embeddings like UMAP, t-SNE, PCA coordinates.',
    'obsp': 'Pairwise annotation of observations. Stores square matrices representing graphs (e.g., cell-cell distances).',
    'varm': 'Multi-dimensional annotation of variables/features. Stores multi-dimensional gene annotations.',
    'varp': 'Pairwise annotation of variables/features. Stores square matrices for gene-gene relationships.',
    'uns': 'Unstructured annotation (dictionary). Stores additional metadata that doesn\'t fit into other slots.',
    'layers': 'Dictionary-like object with values of the same dimensions as X. Stores alternative data matrices (e.g., raw counts, normalized).'
}


@st.dialog("About Explore Dataset", width="large")
def show_explore_info_modal():
    """Display information about the Explore Dataset module."""
    st.markdown("""
    ## What is Explore Dataset?
    
    This module provides an interactive overview of your single-cell dataset's 
    structure, helping you understand what data is available for analysis.
    
    ### AnnData Structure
    
    Single-cell datasets are stored in **AnnData** (Annotated Data) format, 
    which organizes data into several components:
    
    | Component | Description |
    |-----------|-------------|
    | **X** | The main expression matrix (cells × genes) |
    | **obs** | Cell/observation metadata (e.g., cell type, batch) |
    | **var** | Gene/variable metadata (e.g., gene symbols, IDs) |
    | **obsm** | Embeddings like UMAP, t-SNE, PCA coordinates |
    | **obsp** | Cell-cell relationships (e.g., neighbor graphs) |
    | **uns** | Unstructured data (e.g., color palettes, analysis params) |
    | **layers** | Alternative matrices (e.g., raw counts, normalized) |
    
    ### How to Use This Module
    
    1. **Click on panels** to expand and explore each data component
    2. **Hover over buttons** to see detailed descriptions
    3. **Browse tables** to inspect metadata columns and values
    4. **Check available embeddings** in obsm before visualization
    
    ### Tips
    
    - The **obs** panel shows cell-level annotations useful for grouping
    - The **var** panel shows gene metadata for filtering or mapping
    - Check **layers** for raw vs. normalized expression matrices
    """)
    
    st.info("💡 **Tip**: Understanding your data structure helps choose the right analysis parameters in other modules.")


def render(adata: AnnData) -> None:
    """
    Render the Explore Dataset module.
    
    Parameters:
        adata (AnnData): The loaded AnnData object to explore.
    
    This module displays an interactive visual representation of the AnnData
    structure with clickable panels for each data layer (var, obs, X, obsm, etc.)
    """
    # Title with info button
    col_title, col_info = st.columns([6, 1])
    with col_title:
        st.title("📊 Explore Dataset Contents")
    with col_info:
        st.markdown("<div style='margin-top: 12px;'></div>", unsafe_allow_html=True)
        if st.button("ℹ️ Help", key="explore_help_btn", use_container_width=True):
            show_explore_info_modal()
    
    st.markdown("---")
    
    # Show info notification
    st.info("💡 Hover over the buttons to see descriptions. Click on the panels below to explore different layers of the AnnData object.")
    
    # Dataset Overview Cards
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div style='background-color: #f0f2f6; padding: 20px; border-radius: 10px; margin-bottom: 10px;'>
            <h4>📦 Loaded AnnData Object</h4>
            <p><b>Number of observations:</b> {}<br>
            <b>Number of variables:</b> {}</p>
        </div>
        """.format(adata.n_obs, adata.n_vars), unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div style='background-color: #f0f2f6; padding: 20px; border-radius: 10px; margin-bottom: 10px;'>
            <h4>📚 AnnData Object Structure</h4>
            <p>For more information, visit the 
            <a href='https://anndata.readthedocs.io/' target='_blank'>official documentation</a></p>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("---")
    
    # Check which layers are available
    has_obsm = len(adata.obsm) > 0
    has_obsp = len(adata.obsp) > 0
    has_varm = len(adata.varm) > 0
    has_varp = len(adata.varp) > 0
    has_uns = len(adata.uns) > 0
    has_layers = len(adata.layers) > 0
    
    # Custom CSS for the AnnData structure layout and button heights
    # Colors match PyShiny reference app styling
    st.markdown("""
    <style>
    /* Individual button colors matching PyShiny reference */
    .st-key-btn_var button {
        min-height: 6rem !important;
        height: 6rem !important;
        font-size: 1.3rem !important;
        font-weight: bold !important;
        background-color: #2C96C0 !important;
        border-color: #2C96C0 !important;
        color: white !important;
    }
    
    .st-key-btn_obsp button {
        min-height: 10rem !important;
        height: 10rem !important;
        font-size: 1.3rem !important;
        font-weight: bold !important;
        background-color: #F15C5A !important;
        border-color: #F15C5A !important;
        color: white !important;
    }
    
    .st-key-btn_obsm button {
        min-height: 10rem !important;
        height: 10rem !important;
        font-size: 1.3rem !important;
        font-weight: bold !important;
        background-color: #EF9021 !important;
        border-color: #EF9021 !important;
        color: white !important;
    }
    
    .st-key-btn_X button {
        min-height: 10rem !important;
        height: 10rem !important;
        font-size: 1.3rem !important;
        font-weight: bold !important;
        background-color: #4FBA6F !important;
        border-color: #4FBA6F !important;
        color: white !important;
    }
    
    .st-key-btn_obs button {
        min-height: 10rem !important;
        height: 10rem !important;
        font-size: 1.3rem !important;
        font-weight: bold !important;
        background-color: #EFC41C !important;
        border-color: #EFC41C !important;
        color: white !important;
    }
    
    .st-key-btn_varm button {
        min-height: 6rem !important;
        height: 6rem !important;
        font-size: 1.3rem !important;
        font-weight: bold !important;
        background-color: #194C61 !important;
        border-color: #194C61 !important;
        color: white !important;
    }
    
    .st-key-btn_varp button {
        min-height: 10rem !important;
        height: 10rem !important;
        font-size: 1.3rem !important;
        font-weight: bold !important;
        background-color: #965BA5 !important;
        border-color: #965BA5 !important;
        color: white !important;
    }
    
    /* Hover states - slightly darker */
    .st-key-btn_var button:hover { background-color: #247fa3 !important; border-color: #247fa3 !important; }
    .st-key-btn_obsp button:hover { background-color: #d94a48 !important; border-color: #d94a48 !important; }
    .st-key-btn_obsm button:hover { background-color: #d6801d !important; border-color: #d6801d !important; }
    .st-key-btn_X button:hover { background-color: #45a561 !important; border-color: #45a561 !important; }
    .st-key-btn_obs button:hover { background-color: #d6b019 !important; border-color: #d6b019 !important; }
    .st-key-btn_varm button:hover { background-color: #133c4d !important; border-color: #133c4d !important; }
    .st-key-btn_varp button:hover { background-color: #7d4c8a !important; border-color: #7d4c8a !important; }
    
    /* Disabled button states */
    .st-key-btn_obsp button:disabled,
    .st-key-btn_obsm button:disabled,
    .st-key-btn_varm button:disabled,
    .st-key-btn_varp button:disabled {
        background-color: #adb5bd !important;
        border-color: #adb5bd !important;
        opacity: 0.65 !important;
    }
    
    /* Uns button - special styling with lightgray background */
    .st-key-btn_uns button {
        min-height: 6rem !important;
        height: 6rem !important;
        font-size: 1.3rem !important;
        font-weight: bold !important;
        background-color: lightgray !important;
        border-color: lightgray !important;
        color: #333 !important;
    }
    
    .st-key-btn_uns button:hover {
        background-color: #b0b0b0 !important;
        border-color: #b0b0b0 !important;
    }
    
    .st-key-btn_uns button:disabled {
        background-color: #ddd !important;
        border-color: #ddd !important;
        opacity: 0.65 !important;
    }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown("### AnnData Structure")
    st.caption("Click on each component to explore its contents")
    
    # Create the visual AnnData structure layout
    # Following PyShiny reference: var/varm/varp vertically aligned, obs elements horizontally aligned
    # Column widths from PyShiny: (6, 4, 2) for rows 1,3,4 and (4, 2, 4, 2) for row 2
    
    # Row 1: Empty(6) | var(4) | Empty(2) - var aligns with X column
    col1, col2, col3 = st.columns([6, 4, 2])
    with col2:
        if st.button("var", key="btn_var", use_container_width=True, 
                    help=TOOLTIPS['var']):
            show_var_modal(adata)
    
    # Row 2: obsp(4) | obsm(2) | X(4) | obs(2) - horizontally aligned
    col1, col2, col3, col4 = st.columns([4, 2, 4, 2])
    with col1:
        if st.button("obsp", key="btn_obsp", use_container_width=True, 
                    disabled=not has_obsp,
                    help=TOOLTIPS['obsp'] if has_obsp else "No obsp data available"):
            if has_obsp:
                show_obsp_modal(adata)
    with col2:
        if st.button("obsm", key="btn_obsm", use_container_width=True,
                    disabled=not has_obsm,
                    help=TOOLTIPS['obsm'] if has_obsm else "No obsm data available"):
            if has_obsm:
                show_obsm_modal(adata)
    with col3:
        if st.button("X", key="btn_X", use_container_width=True, 
                    help=TOOLTIPS['X']):
            show_X_modal(adata, has_layers)
    with col4:
        if st.button("obs", key="btn_obs", use_container_width=True, 
                    help=TOOLTIPS['obs']):
            show_obs_modal(adata)
    
    # Row 3: uns(6) | varm(4) | Empty(2) - uns separate, varm aligns with X
    col1, col2, col3 = st.columns([6, 4, 2])
    with col1:
        if st.button("{uns}", key="btn_uns", use_container_width=True,
                    disabled=not has_uns,
                    help=TOOLTIPS['uns'] if has_uns else "No uns data available"):
            if has_uns:
                show_uns_modal(adata)
    with col2:
        if st.button("varm", key="btn_varm", use_container_width=True,
                    disabled=not has_varm,
                    help=TOOLTIPS['varm'] if has_varm else "No varm data available"):
            if has_varm:
                show_varm_modal(adata)
    
    # Row 4: Empty(6) | varp(4) | Empty(2) - varp aligns with X/var/varm
    col1, col2, col3 = st.columns([6, 4, 2])
    with col2:
        if st.button("varp", key="btn_varp", use_container_width=True,
                    disabled=not has_varp,
                    help=TOOLTIPS['varp'] if has_varp else "No varp data available"):
            if has_varp:
                show_varp_modal(adata)


###########################################################
##########                  VAR                  ##########
###########################################################

@st.dialog("Table of variables (.var)", width="large")
def show_var_modal(adata: AnnData):
    """Display the var table in a modal dialog."""
    st.markdown(f"**Shape:** {adata.n_vars} variables (genes/features)")
    st.caption(TOOLTIPS['var'])
    st.markdown("---")
    
    var_df = adata.var.reset_index()
    if 'index' in var_df.columns:
        var_df = var_df.rename(columns={'index': 'Gene identifier used in app'})
    
    st.dataframe(
        var_df,
        use_container_width=True,
        height=500
    )


###########################################################
##########                  OBS                  ##########
###########################################################

def get_categorical_columns(adata: AnnData) -> List[str]:
    """Get list of categorical columns from obs."""
    categorical_cols = []
    for col in adata.obs.columns:
        if adata.obs[col].nunique() <= MAX_N_CATEGORIES:
            categorical_cols.append(col)
    return categorical_cols


@st.dialog("Table of observations (.obs) - summarized", width="large")
def show_obs_modal(adata: AnnData):
    """Display the obs table summary in a modal dialog."""
    st.markdown(f"**Shape:** {adata.n_obs} observations (cells/samples)")
    st.caption(TOOLTIPS['obs'])
    st.markdown("---")
    
    data_list = []
    categorical_cols = get_categorical_columns(adata)
    
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        
        if col in categorical_cols:
            unique = ", ".join([str(x) for x in adata.obs[col].unique().tolist()])
        else:
            unique = f"Too many unique values to display (limit = {MAX_N_CATEGORIES})"
        
        data_list.append([col, n_unique, unique])
    
    obs_df = pd.DataFrame(data_list, columns=["Column", "No. unique", "Unique categories"])
    
    st.dataframe(
        obs_df,
        use_container_width=True,
        height=500
    )


###########################################################
##########                   X                   ##########
###########################################################

def get_count_sample(adata: AnnData, layer_name: str = None):
    """Get a sample of the count matrix."""
    if layer_name is None or layer_name == "AnnData.X":
        matrix = adata.X
    else:
        matrix = adata.layers[layer_name]
    
    # Get sample of up to 1000x1000
    max_size = min(1000, matrix.shape[0], matrix.shape[1])
    
    if scipy.sparse.issparse(matrix):
        X_sample = np.array(matrix[:max_size, :max_size].todense())
    else:
        X_sample = np.array(matrix[:max_size, :max_size])
    
    # Sort and select random subset
    X_sample.sort()
    X_sample = X_sample[:, ::-1]
    
    # Select up to 5 random rows and 5 columns
    n_rows = min(5, max_size)
    n_cols = min(5, max_size)
    
    if max_size > 5:
        row_indices = random.sample(range(max_size), n_rows)
        X_sample = X_sample[row_indices, :n_cols]
    else:
        X_sample = X_sample[:n_rows, :n_cols]
    
    return pd.DataFrame(X_sample)


@st.dialog("Count matrices available", width="large")
def show_X_modal(adata: AnnData, has_layers: bool):
    """Display count matrix samples in a modal dialog."""
    st.markdown(f"**Data matrix shape:** {adata.n_obs} observations × {adata.n_vars} variables")
    st.caption(TOOLTIPS['X'])
    
    if has_layers:
        st.info(f"💡 This dataset has {len(adata.layers)} additional layer(s): {', '.join(adata.layers.keys())}")
    
    st.markdown("---")
    
    # Get available layers
    layer_options = ["AnnData.X"] + list(adata.layers.keys())
    
    selected_layer = st.radio(
        "Choose from available count layers",
        options=layer_options,
        horizontal=True,
        key="layer_selector",
        help=TOOLTIPS['layers'] if has_layers else "Main data matrix"
    )
    
    st.markdown("---")
    st.subheader("Sample of counts in chosen layer")
    st.caption("(Maximum of 5 columns and 5 rows shown)")
    
    count_sample = get_count_sample(adata, selected_layer if selected_layer != "AnnData.X" else None)
    st.code(count_sample.to_string(), language=None)


###########################################################
##########                 OBSP                  ##########
###########################################################

def get_matrix_sample(matrix):
    """Get a sample of a matrix (obsp/obsm/varm/varp)."""
    max_size = min(1000, matrix.shape[0], matrix.shape[1] if matrix.ndim > 1 else 1000)
    
    if scipy.sparse.issparse(matrix):
        if matrix.ndim > 1:
            element_sample = np.array(matrix[:max_size, :max_size].todense())
        else:
            element_sample = np.array(matrix[:max_size].todense())
    else:
        if matrix.ndim > 1:
            element_sample = np.array(matrix[:max_size, :max_size])
        else:
            element_sample = np.array(matrix[:max_size])
    
    # Handle 1D arrays
    if element_sample.ndim == 1:
        element_sample = element_sample.reshape(-1, 1)
    
    # Sort and select random subset
    element_sample.sort()
    element_sample = element_sample[:, ::-1]
    
    # Select up to 5 random rows and 5 columns
    n_rows = min(5, element_sample.shape[0])
    n_cols = min(5, element_sample.shape[1])
    
    if element_sample.shape[0] > 5:
        row_indices = random.sample(range(element_sample.shape[0]), n_rows)
        element_sample = element_sample[row_indices, :n_cols]
    else:
        element_sample = element_sample[:n_rows, :n_cols]
    
    return pd.DataFrame(element_sample)


@st.dialog("Obsp elements available", width="large")
def show_obsp_modal(adata: AnnData):
    """Display obsp elements in a modal dialog."""
    st.markdown(f"**Number of elements:** {len(adata.obsp)}")
    st.caption(TOOLTIPS['obsp'])
    st.markdown("---")
    
    obsp_options = list(adata.obsp.keys())
    
    selected_obsp = st.radio(
        "Choose from available obsp elements",
        options=obsp_options,
        horizontal=True,
        key="obsp_selector"
    )
    
    if selected_obsp:
        element_shape = adata.obsp[selected_obsp].shape
        st.info(f"**Shape:** {element_shape[0]} × {element_shape[1]}")
    
    st.markdown("---")
    st.subheader("Sample of values in chosen obsp element")
    st.caption("(Maximum of 5 columns and 5 rows shown)")
    
    element_sample = get_matrix_sample(adata.obsp[selected_obsp])
    st.code(element_sample.to_string(), language=None)


###########################################################
##########                 OBSM                  ##########
###########################################################

@st.dialog("Obsm elements available", width="large")
def show_obsm_modal(adata: AnnData):
    """Display obsm elements in a modal dialog."""
    st.markdown(f"**Number of elements:** {len(adata.obsm)}")
    st.caption(TOOLTIPS['obsm'])
    st.markdown("---")
    
    obsm_options = list(adata.obsm.keys())
    
    selected_obsm = st.radio(
        "Choose from available obsm elements",
        options=obsm_options,
        horizontal=True,
        key="obsm_selector"
    )
    
    if selected_obsm:
        element_shape = adata.obsm[selected_obsm].shape
        st.info(f"**Shape:** {element_shape[0]} observations × {element_shape[1]} dimensions")
    
    st.markdown("---")
    st.subheader("Sample of values in chosen obsm element")
    st.caption("(Maximum of 5 columns and 5 rows shown)")
    
    element_sample = get_matrix_sample(adata.obsm[selected_obsm])
    st.code(element_sample.to_string(), language=None)


###########################################################
##########                 VARM                  ##########
###########################################################

@st.dialog("Varm elements available", width="large")
def show_varm_modal(adata: AnnData):
    """Display varm elements in a modal dialog."""
    st.markdown(f"**Number of elements:** {len(adata.varm)}")
    st.caption(TOOLTIPS['varm'])
    st.markdown("---")
    
    varm_options = list(adata.varm.keys())
    
    selected_varm = st.radio(
        "Choose from available varm elements",
        options=varm_options,
        horizontal=True,
        key="varm_selector"
    )
    
    if selected_varm:
        element_shape = adata.varm[selected_varm].shape
        st.info(f"**Shape:** {element_shape[0]} variables × {element_shape[1]} dimensions")
    
    st.markdown("---")
    st.subheader("Sample of values in chosen varm element")
    st.caption("(Maximum of 5 columns and 5 rows shown)")
    
    element_sample = get_matrix_sample(adata.varm[selected_varm])
    st.code(element_sample.to_string(), language=None)


###########################################################
##########                 VARP                  ##########
###########################################################

@st.dialog("Varp elements available", width="large")
def show_varp_modal(adata: AnnData):
    """Display varp elements in a modal dialog."""
    st.markdown(f"**Number of elements:** {len(adata.varp)}")
    st.caption(TOOLTIPS['varp'])
    st.markdown("---")
    
    varp_options = list(adata.varp.keys())
    
    selected_varp = st.radio(
        "Choose from available varp elements",
        options=varp_options,
        horizontal=True,
        key="varp_selector"
    )
    
    if selected_varp:
        element_shape = adata.varp[selected_varp].shape
        st.info(f"**Shape:** {element_shape[0]} × {element_shape[1]}")
    
    st.markdown("---")
    st.subheader("Sample of values in chosen varp element")
    st.caption("(Maximum of 5 columns and 5 rows shown)")
    
    element_sample = get_matrix_sample(adata.varp[selected_varp])
    st.code(element_sample.to_string(), language=None)


###########################################################
##########                 UNS                   ##########
###########################################################

@st.dialog("Uns elements available", width="large")
def show_uns_modal(adata: AnnData):
    """Display uns elements in a modal dialog."""
    st.markdown(f"**Number of elements:** {len(adata.uns)}")
    st.caption(TOOLTIPS['uns'])
    st.markdown("---")
    
    uns_options = list(adata.uns.keys())
    
    selected_uns = st.radio(
        "Choose from available uns elements",
        options=uns_options,
        horizontal=True,
        key="uns_selector"
    )
    
    st.markdown("---")
    st.subheader("Values in chosen uns element")
    
    element = adata.uns[selected_uns]
    
    # Display type information
    st.caption(f"**Type:** {type(element).__name__}")
    
    # Try to format as JSON if it's a dict
    if isinstance(element, dict):
        try:
            element_str = json.dumps(element, indent=4)
        except TypeError:
            element_str = str(element)
    else:
        element_str = str(element)
    
    st.code(element_str, language=None)
