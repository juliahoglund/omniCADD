# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
project = 'omniCADD'
copyright = '2026, omniCADD developers'
author = 'omniCADD developers'
release = '1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'myst_parser',  # For markdown support
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

# Support both .rst and .md files
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'  # ReadTheDocs theme
html_static_path = ['_static']

# Custom CSS to override theme defaults
html_css_files = [
    'custom.css',
]

# Copy additional files to the build directory
html_extra_path = []

# Allow raw HTML in RST files
html_context = {
    'display_github': True,
}

# -- MyST markdown options ---------------------------------------------------
myst_enable_extensions = [
    "deflist",
    "colon_fence",
]
