# scVIZ - Single-Cell RNA-Seq Visualization Tool

## Project Overview

scVIZ is a Streamlit-based interactive application for exploring and analyzing single-cell RNA-sequencing datasets using Scanpy. The app follows a modular architecture where analysis modules are independent Python scripts in the `modules/` directory.

## Architecture

- **app.py**: Main Streamlit application with home screen, dataset loading, and module routing
- **utils/data_loader.py**: Dataset I/O (load from files/Scanpy examples, save, summarize)
- **modules/**: Analysis modules (explore_dataset.py, visualize_expression.py)
- **tests/**: Pytest tests for core utilities (data loading functions)

Session state (`st.session_state`) stores the loaded AnnData object and dataset name across module navigation.

## Key Technologies

- **Streamlit**: UI framework (session state for data persistence, sidebar for navigation)
- **Scanpy**: scRNA-seq analysis (AnnData objects, dimensionality reduction, plotting)
- **UV**: Dependency management (use `uv pip install -e ".[dev]"` for setup)
- **Pytest**: Testing framework for utilities (not UI modules)

## Development Workflow

```bash
# Setup
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"

# Run app
streamlit run app.py

# Run tests
pytest
```

## Module Development Pattern

New analysis modules should:
1. Be placed in `modules/` directory as separate `.py` files
2. Implement a `render(adata: AnnData) -> None` function
3. Import and call from app.py's module routing logic
4. Use Streamlit components for user interaction
5. Handle edge cases (missing metadata, no dim reduction, etc.)

Example module structure:
```python
def render(adata: AnnData) -> None:
    """Render the module UI."""
    st.title("Module Name")
    # ... module logic
```

## Data Loading Conventions

- Supports h5ad files (AnnData format)
- Scanpy example datasets via `scanpy:dataset_name` syntax (e.g., `scanpy:pbmc3k`)
- Local file uploads through Streamlit file_uploader
- Pre-configured datasets from `data/` directory (gitignored)

## Testing Practices

- Test **utilities** in `utils/` (data loading, transformations)
- Use pytest fixtures for sample AnnData objects
- Mock file I/O when appropriate
- Don't test Streamlit UI components or visualization modules directly

## Common Patterns

- **Sparse matrix handling**: Check `hasattr(adata.X, 'toarray')` before array operations
- **Dimensionality reduction**: Check `adata.obsm.keys()` for `X_umap`, `X_tsne`, `X_pca`
- **Categorical metadata**: Filter by `dtype == 'object' or dtype.name == 'category'`
- **Error handling**: Return `None` from loaders on failure, display `st.error()` in UI

## Python Code Standards

- Follow PEP 8 (4 spaces, 79 char lines)
- Use type hints with `typing` module (`List[str]`, `Dict[str, int]`, `Optional[T]`)
- Write docstrings (PEP 257) with Parameters/Returns sections
- Handle edge cases explicitly (empty inputs, missing metadata, sparse vs dense matrices)