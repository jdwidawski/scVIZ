"""
Dataset loading utilities for scVIZ.

Handles loading single-cell RNA-seq datasets from various sources.
"""
import os
import requests
from typing import Dict, Optional, Union, List
from pathlib import Path
from urllib.parse import urlparse
import pandas as pd
import scanpy as sc
from anndata import AnnData
import tempfile


def get_available_datasets() -> Dict[str, Dict[str, str]]:
    """
    Get a dictionary of available pre-configured datasets.
    
    Returns:
        Dict[str, Dict[str, str]]: Dictionary mapping dataset names to their
            metadata (path, description).
    
    Example:
        >>> datasets = get_available_datasets()
        >>> print(datasets['PBMC3k']['path'])
        'data/pbmc3k.h5ad'
    """
    datasets = {}
    
    # Check for datasets in the data directory
    data_dir = Path("data")
    if data_dir.exists():
        for h5ad_file in data_dir.glob("*.h5ad"):
            dataset_name = h5ad_file.stem
            datasets[dataset_name] = {
                'path': str(h5ad_file),
                'description': f'Dataset: {dataset_name}'
            }
    
    # Add example datasets from scanpy if available
    # Note: These are downloaded on-demand by scanpy
    example_datasets = {
        'PBMC3k (Example)': {
            'path': 'scanpy:pbmc3k',
            'description': '3k PBMCs from 10x Genomics'
        },
        'PBMC68k (Example)': {
            'path': 'scanpy:pbmc68k_reduced',
            'description': '68k PBMCs (reduced) from 10x Genomics'
        }
    }
    
    datasets.update(example_datasets)
    
    return datasets


def load_dataset_from_file(file_path: Union[str, Path, object]) -> Optional[AnnData]:
    """
    Load an AnnData object from a file, URL, or file-like object.
    
    Parameters:
        file_path (Union[str, Path, object]): Path to the h5ad file, URL, or a 
            file-like object (e.g., from Streamlit file uploader). Can also be a 
            special string like 'scanpy:pbmc3k' to load example datasets.
    
    Returns:
        Optional[AnnData]: Loaded AnnData object, or None if loading failed.
    
    Example:
        >>> adata = load_dataset_from_file('data/my_dataset.h5ad')
        >>> if adata is not None:
        ...     print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")
        
        >>> adata = load_dataset_from_file('https://example.com/data.h5ad')
    """
    try:
        # Handle special scanpy dataset loading
        if isinstance(file_path, str) and file_path.startswith('scanpy:'):
            dataset_name = file_path.split(':')[1]
            if dataset_name == 'pbmc3k':
                adata = sc.datasets.pbmc3k()
            elif dataset_name == 'pbmc68k_reduced':
                adata = sc.datasets.pbmc68k_reduced()
            else:
                raise ValueError(f"Unknown scanpy dataset: {dataset_name}")
        
        # Handle URLs
        elif isinstance(file_path, str) and _is_url(file_path):
            adata = sc.read_h5ad(file_path, backup_url=file_path)
        
        # Handle file-like objects (from Streamlit uploader)
        elif hasattr(file_path, 'read'):
            # Check file name extension if available
            file_name = getattr(file_path, 'name', '')
            if file_name.endswith('.h5'):
                adata = sc.read_10x_h5(file_path)
            else:
                adata = sc.read_h5ad(file_path)
        
        # Handle file paths
        else:
            file_path = Path(file_path)
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")
            # Determine file type from extension
            if str(file_path).endswith('.h5'):
                adata = sc.read_10x_h5(file_path)
            else:
                adata = sc.read_h5ad(file_path)
        
        return adata
    
    except Exception as e:
        print(f"Error loading dataset: {e}")
        return None


def _is_url(path: str) -> bool:
    """
    Check if a string is a valid URL.
    
    Parameters:
        path (str): The string to check.
    
    Returns:
        bool: True if the string is a URL, False otherwise.
    """
    try:
        result = urlparse(path)
        return all([result.scheme, result.netloc])
    except:
        return False


def save_dataset(adata: AnnData, file_path: Union[str, Path]) -> bool:
    """
    Save an AnnData object to a file.
    
    Parameters:
        adata (AnnData): The AnnData object to save.
        file_path (Union[str, Path]): Path where the file should be saved.
    
    Returns:
        bool: True if save was successful, False otherwise.
    
    Example:
        >>> success = save_dataset(adata, 'data/processed_dataset.h5ad')
        >>> if success:
        ...     print("Dataset saved successfully")
    """
    try:
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(file_path)
        return True
    
    except Exception as e:
        print(f"Error saving dataset: {e}")
        return False


def get_dataset_summary(adata: AnnData) -> Dict[str, any]:
    """
    Generate a summary of key dataset statistics.
    
    Parameters:
        adata (AnnData): The AnnData object to summarize.
    
    Returns:
        Dict[str, any]: Dictionary containing summary statistics.
    
    Example:
        >>> summary = get_dataset_summary(adata)
        >>> print(f"Total cells: {summary['n_cells']}")
    """
    summary = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'obs_keys': list(adata.obs.columns),
        'var_keys': list(adata.var.columns),
        'uns_keys': list(adata.uns.keys()) if adata.uns else [],
        'obsm_keys': list(adata.obsm.keys()) if adata.obsm else [],
    }
    
    return summary


def load_dataset_table(file_path: Union[str, Path] = "datasets.csv") -> pd.DataFrame:
    """
    Load the dataset table from a CSV file.
    
    Parameters:
        file_path (Union[str, Path]): Path to the CSV file containing dataset 
            information. Defaults to "datasets.csv" in the project root.
    
    Returns:
        pd.DataFrame: DataFrame with columns 'Dataset Name' and 'Path/URL'.
            Returns empty DataFrame with correct schema if file doesn't exist.
    
    Example:
        >>> df = load_dataset_table()
        >>> print(df.columns)
        Index(['Dataset Name', 'Path/URL'], dtype='object')
    """
    file_path = Path(file_path)
    
    if file_path.exists():
        try:
            df = pd.read_csv(file_path, dtype=str)
            # Ensure correct columns exist
            if 'Dataset Name' not in df.columns or 'Path/URL' not in df.columns:
                print(f"Warning: {file_path} missing required columns. Creating empty table.")
                return pd.DataFrame({
                    'Dataset Name': pd.Series([], dtype='str'),
                    'Path/URL': pd.Series([], dtype='str')
                })
            return df
        except Exception as e:
            print(f"Error loading dataset table from {file_path}: {e}")
            return pd.DataFrame({
                'Dataset Name': pd.Series([], dtype='str'),
                'Path/URL': pd.Series([], dtype='str')
            })
    else:
        # Return empty DataFrame with correct schema
        return pd.DataFrame({
            'Dataset Name': pd.Series([], dtype='str'),
            'Path/URL': pd.Series([], dtype='str')
        })


def save_dataset_table(df: pd.DataFrame, file_path: Union[str, Path] = "datasets.csv") -> bool:
    """
    Save the dataset table to a CSV file.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing dataset information with columns
            'Dataset Name' and 'Path/URL'.
        file_path (Union[str, Path]): Path where the CSV file should be saved.
            Defaults to "datasets.csv" in the project root.
    
    Returns:
        bool: True if save was successful, False otherwise.
    
    Example:
        >>> df = pd.DataFrame({
        ...     'Dataset Name': ['My Dataset'],
        ...     'Path/URL': ['path/to/data.h5ad']
        ... })
        >>> success = save_dataset_table(df)
        >>> if success:
        ...     print("Table saved successfully")
    """
    try:
        file_path = Path(file_path)
        # Ensure parent directory exists
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Remove empty rows before saving
        df_clean = df.dropna(how='all')
        df_clean = df_clean[df_clean['Dataset Name'].str.strip() != '']
        
        df_clean.to_csv(file_path, index=False)
        return True
    
    except Exception as e:
        print(f"Error saving dataset table to {file_path}: {e}")
        return False


def fetch_cellxgene_datasets(limit: int = 100) -> Optional[pd.DataFrame]:
    """
    Fetch available datasets from the CELLxGENE Curation API.
    
    Uses the CELLxGENE Curation REST API to retrieve dataset metadata including
    dataset IDs, collection names, assays, tissue types, and cell counts.
    
    Parameters:
        limit (int): Maximum number of datasets to retrieve. Defaults to 100.
    
    Returns:
        Optional[pd.DataFrame]: DataFrame with columns: 'Collection ID', 'Dataset ID', 
            'Collection', 'Title', 'Assay', 'Tissue', 'Disease', 'Cell Count'. 
            Returns None if the fetch fails.
    
    Example:
        >>> df = fetch_cellxgene_datasets(limit=10)
        >>> if df is not None:
        ...     print(f"Found {len(df)} datasets")
    """
    try:
        # CELLxGENE Curation API base URL
        api_url_base = "https://api.cellxgene.cziscience.com"
        datasets_path = "/curation/v1/datasets"
        datasets_url = f"{api_url_base}{datasets_path}"
        
        # Fetch all public datasets
        headers = {"Content-Type": "application/json"}
        response = requests.get(url=datasets_url, headers=headers)
        response.raise_for_status()
        datasets = response.json()
        
        # Limit the number of results
        if len(datasets) > limit:
            datasets = datasets[:limit]
        
        # Helper function to safely format list columns
        def format_list_column(value):
            if isinstance(value, (list, tuple)):
                # Extract labels from ontology objects if present
                if len(value) > 0 and isinstance(value[0], dict):
                    return ', '.join(item.get('label', str(item)) for item in value)
                return ', '.join(str(v) for v in value)
            return str(value) if value is not None else 'N/A'
        
        # Extract dataset information
        datasets_list = []
        
        for dataset in datasets:
            collection_id = dataset.get('collection_id', 'N/A')
            dataset_id = dataset.get('dataset_id', 'N/A')
            collection_name = dataset.get('collection_name', 'N/A')
            dataset_title = dataset.get('title', 'N/A')
            
            # Extract metadata fields
            assay = format_list_column(dataset.get('assay', []))
            tissue = format_list_column(dataset.get('tissue', []))
            disease = format_list_column(dataset.get('disease', []))
            organism = format_list_column(dataset.get('organism', []))
            cell_count = dataset.get('cell_count', 0)
            mean_genes = dataset.get('mean_genes_per_cell', 0)
            suspension_type = format_list_column(dataset.get('suspension_type', []))
            
            datasets_list.append({
                'Collection ID': collection_id,
                'Dataset ID': dataset_id,
                'Collection': collection_name,
                'Title': dataset_title,
                'Assay': assay,
                'Tissue': tissue,
                'Disease': disease,
                'Cell Count': cell_count,
                'Mean Genes/Cell': int(mean_genes) if mean_genes else 'N/A',
                'Organism': organism,
                'Suspension Type': suspension_type
            })
        
        # Create DataFrame
        display_df = pd.DataFrame(datasets_list)
        
        # Filter to only include Homo sapiens datasets
        if 'Organism' in display_df.columns and len(display_df) > 0:
            display_df = display_df[
                display_df['Organism'].str.contains('Homo sapiens', case=False, na=False)
            ]
        
        # Sort by cell count descending
        if 'Cell Count' in display_df.columns and len(display_df) > 0:
            display_df = display_df.sort_values('Cell Count', ascending=False)
        
        return display_df
    
    except Exception as e:
        print(f"Error fetching CELLxGENE datasets: {e}")
        import traceback
        traceback.print_exc()
        return None


def load_cellxgene_dataset(
    collection_id: str,
    dataset_id: str
) -> Optional[AnnData]:
    """
    Load a dataset from CELLxGENE by downloading the h5ad file.
    
    Uses the CELLxGENE Curation API to get dataset metadata and download
    the h5ad asset file.
    
    Parameters:
        collection_id (str): The CELLxGENE collection ID.
        dataset_id (str): The CELLxGENE dataset ID to load.
    
    Returns:
        Optional[AnnData]: Loaded AnnData object, or None if loading failed.
    
    Example:
        >>> adata = load_cellxgene_dataset('col-123', 'dataset-456')
        >>> if adata is not None:
        ...     print(f"Loaded {adata.n_obs} cells")
    """
    try:
        # CELLxGENE Curation API base URL
        api_url_base = "https://api.cellxgene.cziscience.com"
        dataset_path = f"/curation/v1/collections/{collection_id}/datasets/{dataset_id}"
        dataset_url = f"{api_url_base}{dataset_path}"
        
        # Fetch dataset metadata
        response = requests.get(url=dataset_url)
        response.raise_for_status()
        dataset_metadata = response.json()
        
        # Get the h5ad asset URL from the assets list
        assets = dataset_metadata.get('assets', [])
        h5ad_url = None
        
        for asset in assets:
            if asset.get('filetype') == 'H5AD':
                h5ad_url = asset.get('url')
                break
        
        if not h5ad_url:
            print(f"No H5AD file found for dataset {dataset_id}")
            return None
        
        # Download the h5ad file to a temporary location
        print(f"Downloading dataset from {h5ad_url}...")
        with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as tmp_file:
            tmp_path = tmp_file.name
            
            # Stream download to handle large files
            with requests.get(h5ad_url, stream=True) as r:
                r.raise_for_status()
                for chunk in r.iter_content(chunk_size=8192):
                    tmp_file.write(chunk)
        
        # Load the AnnData object from the temporary file
        print(f"Loading dataset from temporary file...")
        adata = sc.read_h5ad(tmp_path)
        
        # Clean up temporary file
        os.unlink(tmp_path)
        
        return adata
    
    except Exception as e:
        print(f"Error loading CELLxGENE dataset: {e}")
        import traceback
        traceback.print_exc()
        return None

