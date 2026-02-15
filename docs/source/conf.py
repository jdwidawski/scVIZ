# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
# Copyright (c) 2026 Jakub Widawski
# Licensed under the MIT License.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'scVIZ'
copyright = '2026, Jakub Widawski'
author = 'Jakub Widawski'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# Support for Markdown files
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# GitHub Pages configuration
html_baseurl = 'https://jdwidawski.github.io/scVIZ/'

# Theme options for ReadTheDocs theme
html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

# Custom CSS
html_css_files = [
    'custom.css',
]

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}
autodoc_mock_imports = ['streamlit', 'st_aggrid', 'plotly', 'pydeseq2']

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'scanpy': ('https://scanpy.readthedocs.io/en/stable/', None),
    'anndata': ('https://anndata.readthedocs.io/en/stable/', None),
}
