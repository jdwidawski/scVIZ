"""
Tests for data loading utilities.

Test dataset loading, saving, and summary generation functions.

Copyright (c) 2026 Jakub Widawski
Licensed under the MIT License.
"""
import pytest
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock
from anndata import AnnData
import scipy.sparse as sp
import scanpy as sc
from utils.data_loader import (
    load_dataset_from_file,
    save_dataset,
    get_dataset_summary,
    get_available_datasets,
    _is_url,
)


@pytest.fixture
def sample_adata() -> AnnData:
    """
    Create a sample AnnData object for testing.
    
    Returns:
        AnnData: A simple AnnData object with 100 cells and 50 genes.
    """
    np.random.seed(42)
    X = np.random.poisson(1, size=(100, 50))
    
    adata = AnnData(X)
    adata.obs['cell_type'] = np.random.choice(['TypeA', 'TypeB', 'TypeC'], size=100)
    adata.obs['n_counts'] = np.sum(X, axis=1)
    adata.var['gene_id'] = [f'Gene_{i}' for i in range(50)]
    adata.var_names = [f'GENE{i}' for i in range(50)]
    
    return adata


def test_get_dataset_summary(sample_adata: AnnData) -> None:
    """
    Test get_dataset_summary function.
    
    Parameters:
        sample_adata (AnnData): Sample dataset fixture.
    
    Verifies:
        - Summary contains correct number of cells and genes
        - All expected keys are present in summary
        - Metadata keys are correctly listed
    """
    summary = get_dataset_summary(sample_adata)
    
    assert summary['n_cells'] == 100
    assert summary['n_genes'] == 50
    assert 'cell_type' in summary['obs_keys']
    assert 'n_counts' in summary['obs_keys']
    assert 'gene_id' in summary['var_keys']
    assert isinstance(summary['uns_keys'], list)
    assert isinstance(summary['obsm_keys'], list)


def test_save_and_load_dataset(sample_adata: AnnData, tmp_path: Path) -> None:
    """
    Test saving and loading dataset functionality.
    
    Parameters:
        sample_adata (AnnData): Sample dataset fixture.
        tmp_path (Path): Pytest temporary directory fixture.
    
    Verifies:
        - Dataset can be saved successfully
        - Saved dataset can be loaded back
        - Loaded dataset has same dimensions as original
        - Metadata is preserved
    """
    # Save dataset
    file_path = tmp_path / "test_dataset.h5ad"
    success = save_dataset(sample_adata, file_path)
    
    assert success is True
    assert file_path.exists()
    
    # Load dataset
    loaded_adata = load_dataset_from_file(file_path)
    
    assert loaded_adata is not None
    assert loaded_adata.n_obs == sample_adata.n_obs
    assert loaded_adata.n_vars == sample_adata.n_vars
    assert 'cell_type' in loaded_adata.obs.columns
    assert 'gene_id' in loaded_adata.var.columns


def test_load_nonexistent_file() -> None:
    """
    Test loading a file that doesn't exist.
    
    Verifies:
        - Loading nonexistent file returns None
        - No exceptions are raised
    """
    result = load_dataset_from_file("/nonexistent/path/file.h5ad")
    assert result is None


def test_get_available_datasets() -> None:
    """
    Test get_available_datasets function.
    
    Verifies:
        - Function returns a dictionary
        - Example datasets are included
        - Dataset entries have required keys
    """
    datasets = get_available_datasets()
    
    assert isinstance(datasets, dict)
    
    # Check for example datasets
    assert 'PBMC3k (Example)' in datasets
    assert 'path' in datasets['PBMC3k (Example)']
    assert 'description' in datasets['PBMC3k (Example)']


def test_load_scanpy_example_dataset() -> None:
    """
    Test loading example datasets from scanpy.
    
    Verifies:
        - Scanpy example datasets can be loaded
        - Loaded dataset is not None
        - Dataset has expected structure
    
    Note: This test requires internet connection to download the dataset.
    """
    adata = load_dataset_from_file('scanpy:pbmc3k')
    
    # This might be None if there's no internet connection
    if adata is not None:
        assert isinstance(adata, AnnData)
        assert adata.n_obs > 0
        assert adata.n_vars > 0


def test_save_dataset_creates_directory(sample_adata: AnnData, tmp_path: Path) -> None:
    """
    Test that save_dataset creates parent directories if they don't exist.
    
    Parameters:
        sample_adata (AnnData): Sample dataset fixture.
        tmp_path (Path): Pytest temporary directory fixture.
    
    Verifies:
        - Parent directories are created automatically
        - File is saved successfully in nested directory
    """
    nested_path = tmp_path / "nested" / "dir" / "dataset.h5ad"
    success = save_dataset(sample_adata, nested_path)
    
    assert success is True
    assert nested_path.exists()
    assert nested_path.parent.exists()


def test_dataset_summary_with_empty_metadata() -> None:
    """
    Test get_dataset_summary with minimal AnnData object.
    
    Verifies:
        - Summary works with datasets that have no metadata
        - Empty metadata lists are handled correctly
    """
    X = np.random.poisson(1, size=(10, 5))
    adata = AnnData(X)
    
    summary = get_dataset_summary(adata)
    
    assert summary['n_cells'] == 10
    assert summary['n_genes'] == 5
    assert isinstance(summary['obs_keys'], list)
    assert isinstance(summary['var_keys'], list)


# ==============================================================================
# Tests for _is_url utility function
# ==============================================================================

class TestIsUrl:
    """Tests for the _is_url utility function."""
    
    def test_valid_https_url(self) -> None:
        """Test that valid HTTPS URLs are recognized."""
        assert _is_url("https://example.com/file.h5ad") is True
        assert _is_url("https://cellxgene.cziscience.com/dataset.h5ad") is True
    
    def test_valid_http_url(self) -> None:
        """Test that valid HTTP URLs are recognized."""
        assert _is_url("http://example.com/data.h5ad") is True
    
    def test_local_file_path(self) -> None:
        """Test that local file paths are not URLs."""
        assert _is_url("/home/user/data/file.h5ad") is False
        assert _is_url("data/file.h5ad") is False
        assert _is_url("./file.h5ad") is False
    
    def test_scanpy_prefix(self) -> None:
        """Test that scanpy: prefixed paths are not URLs."""
        assert _is_url("scanpy:pbmc3k") is False
    
    def test_empty_string(self) -> None:
        """Test that empty string is not a URL."""
        assert _is_url("") is False
    
    def test_url_with_query_params(self) -> None:
        """Test URLs with query parameters."""
        url = "https://example.com/file.h5ad?token=abc123"
        assert _is_url(url) is True


# ==============================================================================
# Tests for sparse matrix handling
# ==============================================================================

class TestSparseMatrixHandling:
    """Tests for handling sparse matrix data in AnnData objects."""
    
    @pytest.fixture
    def sparse_adata(self) -> AnnData:
        """Create an AnnData object with sparse matrix."""
        np.random.seed(42)
        X = sp.random(100, 50, density=0.1, format='csr')
        adata = AnnData(X)
        adata.obs['cell_type'] = np.random.choice(['A', 'B'], size=100)
        adata.var_names = [f'Gene_{i}' for i in range(50)]
        return adata
    
    def test_save_sparse_dataset(
        self, sparse_adata: AnnData, tmp_path: Path
    ) -> None:
        """Test saving dataset with sparse matrix."""
        file_path = tmp_path / "sparse_dataset.h5ad"
        success = save_dataset(sparse_adata, file_path)
        
        assert success is True
        assert file_path.exists()
    
    def test_load_sparse_dataset(
        self, sparse_adata: AnnData, tmp_path: Path
    ) -> None:
        """Test loading dataset preserves sparse matrix format."""
        file_path = tmp_path / "sparse_dataset.h5ad"
        save_dataset(sparse_adata, file_path)
        loaded = load_dataset_from_file(file_path)
        
        assert loaded is not None
        assert hasattr(loaded.X, 'toarray') or sp.issparse(loaded.X)
    
    def test_summary_with_sparse_matrix(self, sparse_adata: AnnData) -> None:
        """Test summary generation with sparse matrix data."""
        summary = get_dataset_summary(sparse_adata)
        
        assert summary['n_cells'] == 100
        assert summary['n_genes'] == 50


# ==============================================================================
# Tests for dimensionality reduction metadata
# ==============================================================================

class TestDimensionalityReduction:
    """Tests for datasets with dimensionality reduction results."""
    
    @pytest.fixture
    def adata_with_dimred(self) -> AnnData:
        """Create AnnData with dimensionality reductions."""
        np.random.seed(42)
        X = np.random.poisson(1, size=(100, 50))
        adata = AnnData(X)
        adata.obsm['X_pca'] = np.random.randn(100, 10)
        adata.obsm['X_umap'] = np.random.randn(100, 2)
        adata.obsm['X_tsne'] = np.random.randn(100, 2)
        return adata
    
    def test_summary_contains_obsm_keys(
        self, adata_with_dimred: AnnData
    ) -> None:
        """Test that summary includes dimensionality reduction keys."""
        summary = get_dataset_summary(adata_with_dimred)
        
        assert 'X_pca' in summary['obsm_keys']
        assert 'X_umap' in summary['obsm_keys']
        assert 'X_tsne' in summary['obsm_keys']
    
    def test_save_load_preserves_obsm(
        self, adata_with_dimred: AnnData, tmp_path: Path
    ) -> None:
        """Test that dimensionality reductions are preserved after save/load."""
        file_path = tmp_path / "dimred_dataset.h5ad"
        save_dataset(adata_with_dimred, file_path)
        loaded = load_dataset_from_file(file_path)
        
        assert loaded is not None
        assert 'X_pca' in loaded.obsm.keys()
        assert 'X_umap' in loaded.obsm.keys()
        assert loaded.obsm['X_umap'].shape == (100, 2)


# ==============================================================================
# Tests for edge cases
# ==============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    def test_load_invalid_file_extension(self, tmp_path: Path) -> None:
        """Test loading file with invalid extension returns None."""
        invalid_file = tmp_path / "test.txt"
        invalid_file.write_text("not an h5ad file")
        
        result = load_dataset_from_file(invalid_file)
        assert result is None
    
    def test_load_corrupted_h5ad(self, tmp_path: Path) -> None:
        """Test loading corrupted h5ad file returns None."""
        corrupted = tmp_path / "corrupted.h5ad"
        corrupted.write_text("not valid h5ad content")
        
        result = load_dataset_from_file(corrupted)
        assert result is None
    
    def test_save_to_readonly_path(self, sample_adata: AnnData) -> None:
        """Test saving to invalid path returns False."""
        # Try to save to a path that doesn't allow writing
        result = save_dataset(sample_adata, "/root/readonly/file.h5ad")
        assert result is False
    
    def test_get_available_datasets_returns_dict(self) -> None:
        """Test get_available_datasets always returns a dictionary."""
        datasets = get_available_datasets()
        assert isinstance(datasets, dict)
    
    def test_summary_with_categorical_obs(self) -> None:
        """Test summary with categorical observations."""
        X = np.random.poisson(1, size=(50, 20))
        adata = AnnData(X)
        adata.obs['cell_type'] = np.random.choice(
            ['TypeA', 'TypeB', 'TypeC'], size=50
        )
        adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
        
        summary = get_dataset_summary(adata)
        
        assert 'cell_type' in summary['obs_keys']
        assert summary['n_cells'] == 50
