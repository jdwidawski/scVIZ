# scVIZ Documentation

This directory contains the Sphinx documentation for scVIZ.

## Building Locally

### Prerequisites

Install documentation dependencies:

```bash
pip install -r requirements.txt
```

### Build HTML Documentation

```bash
cd docs
sphinx-build -b html source build/html
```

The generated documentation will be in `build/html/`. Open `build/html/index.html` in your browser.

### Quick Build Command

From the project root:

```bash
sphinx-build -b html docs/source docs/build/html
```

## Cleaning Build Files

Remove generated files:

```bash
rm -rf build/
```

## Documentation Structure

```
docs/
├── source/
│   ├── index.rst              # Main documentation page
│   ├── conf.py                # Sphinx configuration
│   ├── user_guide/            # User guides
│   │   ├── getting_started.rst
│   │   ├── explore_module.rst
│   │   ├── visualize_module.rst
│   │   └── de_module.rst
│   └── api/                   # API reference
│       ├── modules.rst
│       └── utils.rst
├── build/                     # Generated HTML (gitignored)
└── requirements.txt           # Sphinx dependencies
```

## GitHub Pages Deployment

Documentation is automatically built and deployed to GitHub Pages on every push to `main` branch via GitHub Actions.

### Manual Deployment Setup

If setting up for the first time:

1. Go to your repository Settings → Pages
2. Set Source to "GitHub Actions"
3. Push to `main` branch to trigger deployment
4. Documentation will be available at: `https://jdwidawski.github.io/scVIZ/`

### Workflow

The `.github/workflows/docs.yml` workflow:

1. Checks out the repository
2. Installs Python and Sphinx
3. Builds HTML documentation
4. Deploys to GitHub Pages

## Adding New Pages

1. Create a new `.rst` file in the appropriate directory
2. Add the file to a `toctree` directive in the parent page
3. Build locally to verify
4. Commit and push to deploy

Example:

```rst
.. toctree::
   :maxdepth: 2
   
   new_page
```

## Writing Documentation

### reStructuredText Syntax

Headers:
```rst
Title
=====

Section
-------

Subsection
~~~~~~~~~~
```

Code blocks:
```rst
.. code-block:: python

   import scanpy as sc
   adata = sc.read_h5ad("data.h5ad")
```

Links:
```rst
:doc:`other_page`
:ref:`section-label`
```

### Autodoc

API documentation uses Sphinx autodoc to extract docstrings:

```rst
.. automodule:: modules.explore_dataset
   :members:
```

Ensure Python modules have proper docstrings following PEP 257.

## Theme Customization

The documentation uses the Read the Docs theme. Customize in `source/conf.py`:

```python
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth': 4,
}
```

## Troubleshooting

### ModuleNotFoundError

If autodoc fails to import modules, ensure the project root is in the Python path:

```python
# conf.py
import sys
import os
sys.path.insert(0, os.path.abspath('../..'))
```

### Build Warnings

Fix warnings by:
- Checking all file paths in `toctree` directives
- Ensuring proper reStructuredText syntax
- Verifying all cross-references exist

### GitHub Pages Not Updating

1. Check Actions tab for build errors
2. Verify Pages settings (Source: GitHub Actions)
3. Check that `.nojekyll` file is created in build output
4. Wait a few minutes for cache to clear
