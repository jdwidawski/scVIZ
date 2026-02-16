"""
scVIZ - Single-cell RNA-Sequencing Data Visualization Tool

Main Streamlit application entry point.
"""
import streamlit as st
import pandas as pd
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, DataReturnMode
from utils.data_loader import (
    load_dataset_from_file, 
    get_available_datasets,
    fetch_cellxgene_datasets,
    load_cellxgene_dataset
)
from modules import explore_dataset, visualize_dataset, differential_expression, quality_control

# Load example datasets from CSV file
# Paths in the CSV are relative to the root directory
def load_example_datasets() -> pd.DataFrame:
    """Load example datasets from CSV file."""
    csv_path = './data/available_datasets.csv'
    try:
        df = pd.read_csv(csv_path)
        
        # Map CSV columns to app columns
        column_mapping = {
            'dataset_title': 'Name',
            'dataset_description': 'Description',
            'tissue': 'Tissue',
            'n_cells': 'Cells',
            'n_genes': 'Genes',
            'data_path': 'URL'
        }
        
        # Rename columns if they exist
        df = df.rename(columns={k: v for k, v in column_mapping.items() if k in df.columns})
        
        # Add Source column if not present (extract from description or use default)
        if 'Source' not in df.columns:
            # Try to extract source from description (look for DOI or author names)
            def extract_source(desc):
                if 'DOI:' in str(desc):
                    # Extract paper reference before DOI
                    return str(desc).split('.')[0] + '.'
                return 'Local Dataset'
            df['Source'] = df['Description'].apply(extract_source) if 'Description' in df.columns else 'Local Dataset'
        
        # Format cell counts with commas
        if 'Cells' in df.columns:
            df['Cells'] = df['Cells'].apply(
                lambda x: f"~{int(x):,}" if pd.notna(x) and str(x).replace('.','').isdigit() else str(x)
            )
        
        # Format gene counts with commas
        if 'Genes' in df.columns:
            df['Genes'] = df['Genes'].apply(
                lambda x: f"{int(x):,}" if pd.notna(x) and str(x).replace('.','').isdigit() else str(x)
            )
        
        # Ensure all required columns exist
        required_cols = ['Name', 'Cells', 'Tissue', 'Description', 'Source', 'URL']
        for col in required_cols:
            if col not in df.columns:
                df[col] = 'N/A'
        
        return df
    except FileNotFoundError:
        st.warning(f"Could not find {csv_path}. Using empty dataset list.")
        return pd.DataFrame(columns=['Name', 'Cells', 'Tissue', 'Genes', 'Description', 'Source', 'URL'])
    except Exception as e:
        st.error(f"Error loading example datasets: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame(columns=['Name', 'Cells', 'Tissue', 'Genes', 'Description', 'Source', 'URL'])

EXAMPLE_DATASETS = load_example_datasets()


def apply_custom_css():
    """Apply custom CSS styling inspired by modern web design."""
    st.markdown("""
    <style>
    /* Import Inter font */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    
    /* Root variables */
    :root {
        --accent-primary: #3b82f6;
        --accent-secondary: #8b5cf6;
        --accent-hover: #2563eb;
        --text-primary: #1a202c;
        --text-secondary: #4a5568;
        --text-muted: #718096;
        --border-color: #e2e8f0;
        --bg-secondary: #f8f9fa;
    }
    
    /* Global font styling */
    html, body, [class*="css"] {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    }
    
    /* Main content area */
    .main .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
    }
    
    /* Headers styling */
    h1 {
        font-weight: 700 !important;
        color: var(--text-primary) !important;
    }
    
    h2, h3 {
        font-weight: 600 !important;
        color: var(--text-primary) !important;
    }
    
    /* Sidebar styling */
    [data-testid="stSidebar"] {
        background-color: var(--bg-secondary);
        border-right: 1px solid var(--border-color);
    }
    
    [data-testid="stSidebar"] .stMarkdown h1 {
        background: linear-gradient(135deg, var(--accent-primary), var(--accent-secondary));
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
    }
    
    /* Primary buttons */
    .stButton > button[kind="primary"],
    .stButton > button[data-testid="baseButton-primary"] {
        background-color: var(--accent-primary) !important;
        border: none !important;
        border-radius: 8px !important;
        font-weight: 600 !important;
        transition: all 0.3s ease !important;
    }
    
    .stButton > button[kind="primary"]:hover,
    .stButton > button[data-testid="baseButton-primary"]:hover {
        background-color: var(--accent-hover) !important;
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(59, 130, 246, 0.3) !important;
    }
    
    /* Secondary buttons */
    .stButton > button[kind="secondary"],
    .stButton > button[data-testid="baseButton-secondary"] {
        background-color: transparent !important;
        color: var(--accent-primary) !important;
        border: 2px solid var(--accent-primary) !important;
        border-radius: 8px !important;
        font-weight: 600 !important;
        transition: all 0.3s ease !important;
    }
    
    .stButton > button[kind="secondary"]:hover,
    .stButton > button[data-testid="baseButton-secondary"]:hover {
        background-color: var(--accent-primary) !important;
        color: white !important;
    }
    
    /* Metrics styling */
    [data-testid="stMetric"] {
        background-color: var(--bg-secondary);
        padding: 1rem;
        border-radius: 12px;
        border: 1px solid var(--border-color);
    }
    
    [data-testid="stMetricLabel"] {
        font-weight: 500 !important;
        color: var(--text-muted) !important;
    }
    
    [data-testid="stMetricValue"] {
        font-weight: 700 !important;
        color: var(--accent-primary) !important;
    }
    
    /* Expander styling */
    .streamlit-expanderHeader {
        font-weight: 600 !important;
        color: var(--text-primary) !important;
        background-color: var(--bg-secondary) !important;
        border-radius: 8px !important;
    }
    
    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
        background-color: transparent;
    }
    
    .stTabs [data-baseweb="tab"] {
        height: auto;
        padding: 0.75rem 1.5rem;
        background-color: var(--bg-secondary);
        border-radius: 8px;
        font-weight: 500;
        color: var(--text-secondary);
        border: 1px solid var(--border-color);
    }
    
    .stTabs [aria-selected="true"] {
        background-color: var(--accent-primary) !important;
        color: white !important;
        border-color: var(--accent-primary) !important;
    }
    
    /* Input fields */
    .stTextInput > div > div > input,
    .stSelectbox > div > div > div,
    .stMultiSelect > div > div > div {
        border-radius: 8px !important;
        border-color: var(--border-color) !important;
    }
    
    .stTextInput > div > div > input:focus,
    .stSelectbox > div > div > div:focus {
        border-color: var(--accent-primary) !important;
        box-shadow: 0 0 0 2px rgba(59, 130, 246, 0.2) !important;
    }
    
    /* Dataframe styling */
    .stDataFrame {
        border-radius: 12px !important;
        overflow: hidden;
    }
    
    /* Info/Success/Warning/Error boxes */
    .stAlert {
        border-radius: 8px !important;
    }
    
    /* Download button */
    .stDownloadButton > button {
        border-radius: 8px !important;
        font-weight: 600 !important;
    }
    
    /* Divider */
    hr {
        border-color: var(--border-color) !important;
        margin: 1.5rem 0 !important;
    }
    
    /* Radio buttons in sidebar */
    [data-testid="stSidebar"] .stRadio > div {
        gap: 0.5rem;
    }
    
    [data-testid="stSidebar"] .stRadio label {
        padding: 0.75rem 1rem !important;
        border-radius: 8px !important;
        transition: all 0.2s ease !important;
    }
    
    [data-testid="stSidebar"] .stRadio label:hover {
        background-color: rgba(59, 130, 246, 0.1) !important;
    }
    
    /* Success message */
    .element-container .stSuccess {
        background-color: rgba(79, 186, 111, 0.1) !important;
        border-left: 4px solid #4FBA6F !important;
    }
    
    /* Info message */
    .element-container .stInfo {
        background-color: rgba(59, 130, 246, 0.1) !important;
        border-left: 4px solid var(--accent-primary) !important;
    }
    
    /* Dialog/Modal styling */
    [data-testid="stModal"] {
        border-radius: 12px !important;
    }
    
    /* Number input */
    .stNumberInput > div > div > input {
        border-radius: 8px !important;
    }
    
    /* Slider */
    .stSlider > div > div > div > div {
        background-color: var(--accent-primary) !important;
    }
    </style>
    """, unsafe_allow_html=True)


def main() -> None:
    """
    Main application function.
    
    Displays the home screen and handles navigation between modules.
    """
    st.set_page_config(
        page_title="scVIZ",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Apply custom CSS styling
    apply_custom_css()
    
    # Initialize session state
    if 'dataset' not in st.session_state:
        st.session_state.dataset = None
    if 'dataset_name' not in st.session_state:
        st.session_state.dataset_name = None
    
    # Sidebar navigation
    with st.sidebar:
        st.title("🧬 scVIZ")
        st.markdown("---")
        
        # Show dataset info if loaded
        if st.session_state.dataset is not None:
            #st.success(f"**Loaded Dataset**")
            #st.info(f"**{st.session_state.dataset_name}**")
            st.metric("Cells", f"{st.session_state.dataset.n_obs:,}")
            st.metric("Genes", f"{st.session_state.dataset.n_vars:,}")
            
            if st.button("🔄 Load Different Dataset", use_container_width=True):
                st.session_state.dataset = None
                st.session_state.dataset_name = None
                st.rerun()
            
            st.markdown("---")
        
        # Module navigation (only show if dataset is loaded)
        if st.session_state.dataset is not None:
            st.header("Analysis Modules")
            module = st.radio(
                "Select Module:",
                options=["Home", "Explore Dataset Contents", "Quality Control", "Visualize Dataset", "Differential Expression"],
                label_visibility="collapsed"
            )
        else:
            module = "Home"
            st.markdown("""
            **scVIZ** is an interactive tool for exploring and analyzing 
            single-cell RNA-sequencing datasets.
            
            **What you can do:**
            - 📊 Explore dataset structure and metadata
            - 🔍 Quality control and cell filtering
            - 🗺️ Visualize cells on UMAP plots
            - 🧬 Analyze gene expression patterns
            - 📈 Perform differential expression analysis
            
            **Get started** by loading a dataset using one of the options on the right.
            """)
            st.markdown("---")
            st.info("👉 Load a dataset to access analysis modules")
        
        # Copyright footer
        st.markdown("---")
        st.markdown(
            """
            <div style='text-align: center; color: #666; font-size: 0.75rem;'>
            <p>© 2026 <b>Jakub Widawski</b></p>
            <p>Built with Streamlit & Scanpy</p>
            </div>
            """,
            unsafe_allow_html=True
        )
    
    # Main content area
    if st.session_state.dataset is None:
        # Show dataset loading interface
        display_dataset_loading()
    else:
        # Route to appropriate module
        if module == "Home":
            display_home_screen()
        elif module == "Explore Dataset Contents":
            explore_dataset.render(st.session_state.dataset)
        elif module == "Quality Control":
            quality_control.render(st.session_state.dataset)
        elif module == "Visualize Dataset":
            visualize_dataset.render(st.session_state.dataset)
        elif module == "Differential Expression":
            differential_expression.render(st.session_state.dataset)


def display_dataset_loading() -> None:
    """
    Display the dataset loading interface in the main area.
    """
    st.title("🧬 scVIZ - Load Dataset")
    st.markdown("### Choose a dataset to begin your analysis")
    st.markdown("---")
    
    # Create tabs for different loading methods
    tab1, tab2, tab3 = st.tabs([
        "📁 Upload File",
        "📋 Example Datasets",
        "🌐 CELLxGENE Datasets"
    ])
    
    # Tab 1: Upload File
    with tab1:
        st.subheader("Upload Local File")
        st.markdown("Upload an AnnData object in h5ad format from your computer.")
        
        uploaded_file = st.file_uploader(
            "Choose an h5ad or h5 file",
            type=['h5ad', 'h5'],
            help="Upload an AnnData object in h5ad format or 10x h5 format"
        )
        
        if uploaded_file is not None:
            st.success(f"File selected: **{uploaded_file.name}**")
            col1, col2, col3 = st.columns([1, 1, 2])
            with col1:
                if st.button("📥 Load Uploaded File", use_container_width=True, type="primary"):
                    load_dataset(uploaded_file, uploaded_file.name)
    
    # Tab 2: Example Datasets
    with tab2:
        st.subheader("Example Datasets")
        st.markdown("""
        Select from curated example datasets from tutorials and publications.
        All datasets include pre-computed UMAP embeddings and are ready to explore.
        Click a row to select, then click the Load button.
        """)
        
        # Display columns for the table (exclude URL from display)
        display_cols = ['Name', 'Cells', 'Genes', 'Tissue', 'Description', 'Source']
        
        # Configure AgGrid for example datasets
        gb = GridOptionsBuilder.from_dataframe(EXAMPLE_DATASETS[display_cols])
        gb.configure_selection(selection_mode='single', use_checkbox=False)
        gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
        gb.configure_default_column(resizable=True, filterable=True, sortable=True)
        
        # Configure column widths
        gb.configure_column('Name', minWidth=200)
        gb.configure_column('Cells', maxWidth=110)
        gb.configure_column('Genes', maxWidth=110)
        gb.configure_column('Tissue', maxWidth=150)
        gb.configure_column('Description', minWidth=250)
        gb.configure_column('Source', minWidth=180)
        
        grid_options = gb.build()
        
        # Display grid
        grid_response = AgGrid(
            EXAMPLE_DATASETS[display_cols],
            gridOptions=grid_options,
            update_mode=GridUpdateMode.SELECTION_CHANGED,
            data_return_mode=DataReturnMode.FILTERED_AND_SORTED,
            fit_columns_on_grid_load=True,
            height=350,
            theme='streamlit'
        )
        
        # Handle selection
        selected_rows = grid_response['selected_rows']
        
        if selected_rows is not None and len(selected_rows) > 0:
            selected_name = selected_rows.iloc[0]['Name']
            # Get full row from EXAMPLE_DATASETS to retrieve URL
            selected_row = EXAMPLE_DATASETS[EXAMPLE_DATASETS['Name'] == selected_name].iloc[0]
            
            st.markdown("---")
            st.subheader("Selected Dataset")
            
            # Display dataset info
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Cells", selected_row['Cells'])
            with col2:
                st.metric("Genes", selected_row.get('Genes', 'N/A'))
            with col3:
                st.metric("Tissue", selected_row['Tissue'])
            
            st.caption(f"**Description:** {selected_row['Description']}")
            
            col1, col2, col3 = st.columns([1, 1, 2])
            with col1:
                if st.button("📥 Load Selected Dataset", use_container_width=True, type="primary"):
                    load_dataset(selected_row['URL'], selected_row['Name'])
            with col2:
                st.info(f"Source: `{selected_row['Source']}`")
    
    # Tab 3: CELLxGENE Datasets
    with tab3:
        st.subheader("Browse CELLxGENE Discover Datasets")
        st.markdown("""
        Explore and load datasets from the [CELLxGENE Discover](https://cellxgene.cziscience.com/datasets) portal.
        Click a row to select, then click the Load button.
        """)
        
        # Fetch datasets button
        col1, col2 = st.columns([1, 3])
        with col1:
            if st.button("🔍 Fetch Datasets", use_container_width=True, type="secondary"):
                with st.spinner("Fetching datasets from CELLxGENE..."):
                    cellxgene_df = fetch_cellxgene_datasets(limit=100)
                    if cellxgene_df is not None:
                        st.session_state.cellxgene_datasets = cellxgene_df
                        st.success(f"✅ Found {len(cellxgene_df)} datasets!")
                    else:
                        st.error("❌ Failed to fetch datasets. Please try again.")
        
        # Display datasets if available
        if 'cellxgene_datasets' in st.session_state and st.session_state.cellxgene_datasets is not None:
            df = st.session_state.cellxgene_datasets
            
            # Determine which columns to display
            display_cols = []
            # Always include ID columns first (we'll hide them)
            if 'Dataset ID' in df.columns:
                display_cols.append('Dataset ID')
            if 'Collection ID' in df.columns:
                display_cols.append('Collection ID')
            
            # Add other columns
            for col in ['Collection', 'Title', 'Assay', 'Tissue', 'Disease', 'Organism', 'Cell Count', 'Mean Genes/Cell', 'Suspension Type']:
                if col in df.columns:
                    display_cols.append(col)
            
            # Configure AgGrid
            gb = GridOptionsBuilder.from_dataframe(df[display_cols])
            gb.configure_selection(selection_mode='single', use_checkbox=False)
            gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)
            gb.configure_default_column(resizable=True, filterable=True, sortable=True)
            
            # Hide ID columns
            if 'Dataset ID' in display_cols:
                gb.configure_column('Dataset ID', hide=True)
            if 'Collection ID' in display_cols:
                gb.configure_column('Collection ID', hide=True)
            
            # Format number columns
            if 'Cell Count' in display_cols:
                gb.configure_column('Cell Count', type=['numericColumn', 'numberColumnFilter'], valueFormatter="value.toLocaleString()")
            if 'Mean Genes/Cell' in display_cols:
                gb.configure_column('Mean Genes/Cell', type=['numericColumn', 'numberColumnFilter'])
            
            grid_options = gb.build()
            
            # Display grid
            grid_response = AgGrid(
                df[display_cols],
                gridOptions=grid_options,
                update_mode=GridUpdateMode.SELECTION_CHANGED,
                data_return_mode=DataReturnMode.FILTERED_AND_SORTED,
                fit_columns_on_grid_load=True,
                height=400,
                theme='streamlit'
            )
            
            # Handle selection
            selected_rows = grid_response['selected_rows']
            
            if selected_rows is not None and len(selected_rows) > 0:
                selected_row = selected_rows.iloc[0]
                
                st.markdown("---")
                st.subheader("Selected Dataset")
                
                # Display dataset info
                col1, col2, col3 = st.columns(3)
                with col1:
                    cell_count = selected_row.get('Cell Count', 0)
                    st.metric("Cells", f"{cell_count:,}" if isinstance(cell_count, (int, float)) else str(cell_count))
                with col2:
                    gene_count = selected_row.get('Mean Genes/Cell', 'N/A')
                    st.metric("Mean Genes/Cell", f"{gene_count:,}" if isinstance(gene_count, (int, float)) else str(gene_count))
                with col3:
                    organism = str(selected_row.get('Organism', 'Homo sapiens'))
                    st.metric("Organism", organism)
                
                # Get IDs from selected row
                dataset_id = selected_row.get('Dataset ID')
                collection_id = selected_row.get('Collection ID')
                
                # Build display name
                display_name = selected_row.get('Collection', 'Dataset')
                if 'Title' in selected_row and selected_row.get('Title') not in ['N/A', None, '']:
                    display_name = f"{display_name} - {selected_row.get('Title')}"
                elif 'Tissue' in selected_row and selected_row.get('Tissue') not in ['N/A', None, '']:
                    display_name = f"{display_name} - {selected_row.get('Tissue')}"
                
                col1, col2, col3 = st.columns([1, 1, 2])
                with col1:
                    if st.button("📥 Load Selected Dataset", use_container_width=True, type="primary"):
                        load_cellxgene_dataset_ui(
                            collection_id,
                            dataset_id,
                            display_name
                        )
                with col2:
                    st.info(f"ID: `{dataset_id}`")
        else:
            st.info("👆 Click 'Fetch Datasets' to browse available datasets from CELLxGENE")


def load_dataset(path, name: str) -> None:
    """
    Load a dataset and update session state.
    
    Parameters:
        path: Path, URL, or file-like object to load from
        name (str): Display name for the dataset
    """
    with st.spinner(f"Loading {name}..."):
        adata = load_dataset_from_file(path)
        if adata is not None:
            st.session_state.dataset = adata
            st.session_state.dataset_name = name
            st.success(f"✅ Successfully loaded **{name}**!")
            st.rerun()
        else:
            st.error(f"❌ Failed to load dataset from: {path}")


def load_cellxgene_dataset_ui(collection_id: str, dataset_id: str, name: str) -> None:
    """
    Load a CELLxGENE dataset and update session state.
    
    Parameters:
        collection_id (str): CELLxGENE collection ID
        dataset_id (str): CELLxGENE dataset ID
        name (str): Display name for the dataset
    """
    with st.spinner(f"Loading {name} from CELLxGENE..."):
        st.info("⏳ Downloading dataset from CELLxGENE...")
        adata = load_cellxgene_dataset(collection_id, dataset_id)
        if adata is not None:
            st.session_state.dataset = adata
            st.session_state.dataset_name = name
            st.success(f"✅ Successfully loaded **{name}**!")
            st.rerun()
        else:
            st.error(f"❌ Failed to load dataset from CELLxGENE")


def display_home_screen() -> None:
    """
    Display the home screen with app instructions and dataset info.
    """
    st.title("🧬 Welcome to scVIZ")
    st.markdown("### Interactive Single-Cell RNA-Sequencing Data Visualization")
    
    # Show current dataset info
    st.success(f"✅ **Current Dataset:** {st.session_state.dataset_name}")
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Total Cells", f"{st.session_state.dataset.n_obs:,}")
    with col2:
        st.metric("Total Genes", f"{st.session_state.dataset.n_vars:,}")
    
    st.markdown("---")
    
    st.markdown("## Available Analysis Modules")
    st.markdown("Use the **sidebar** to navigate between different analysis workflows:")
    
    # Module cards
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        ### 📊 Explore Dataset Contents
        
        Examine the internal structure of your AnnData object:
        
        - **obs (Cell Metadata)**: View annotations for each cell (e.g., cell type, sample ID, quality metrics)
        - **var (Gene Metadata)**: Browse gene-level information (e.g., gene symbols, highly variable genes)
        - **X (Expression Matrix)**: Preview the raw or normalized expression values
        - **obsm/varm**: Inspect dimensionality reductions (PCA, UMAP) and other embeddings
        - **uns**: Access unstructured metadata (e.g., analysis parameters, color palettes)
        """)
        
        st.markdown("""
        ### 🔍 Quality Control
        
        Assess dataset quality with key QC metrics:
        
        - **Total counts/genes**: Identify low-quality cells with sparse data
        - **Mitochondrial %**: Detect dying or stressed cells
        - **Ribosomal %**: Assess ribosomal gene content
        - **Cell cycle scoring**: Phase assignment (G1, S, G2/M)
        - **Interactive filtering**: Use lasso selection on UMAP
        """)
        
        st.markdown("""
        ### 🧬 Differential Expression
        
        Perform pseudobulk differential expression analysis using PyDESeq2:
        
        - **Select sample groups**: Choose which cell populations to compare
        - **Create pseudobulk**: Aggregate single-cell data by sample
        - **Run statistical analysis**: Identify differentially expressed genes
        - **Volcano plot**: Visualize results with interactive plots
        - **Export results**: Download DE tables for downstream analysis
        """)
    
    with col2:
        st.markdown("""
        ### 🗺️ Visualize Dataset
        
        Create interactive UMAP visualizations to explore your data:
        
        - **Metadata Visualization**: Color cells by categorical annotations (cell type, condition, batch)
        - **Gene Expression**: Overlay expression values for any gene on the UMAP
        - **Multi-gene comparison**: Select multiple genes and view them side-by-side in tabs
        - **Customization**: Adjust figure size, dot size, colormaps, and color ranges
        - **Gene search**: Find genes by ID or alternative symbols from the var table
        """)
        
        st.markdown("""
        ### 💡 Tips
        
        - All visualizations are **interactive** — hover for details, zoom, and pan
        - Use **percentile-based** color scaling to handle outliers in gene expression
        - Downloaded datasets are **cached** for faster subsequent loads
        - Your loaded dataset **persists** across all modules — no need to reload
        """)


if __name__ == "__main__":
    main()
