# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Hybrid Quantum Solver'
copyright = '2026, Benjamin Hess'
author = 'Benjamin Hess'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
import os
import sys
# Expose the root repository layout context to Sphinx's parsing engine
sys.path.insert(0, os.path.abspath('../../'))

# Core Sphinx extension requirements
extensions = [
    'sphinx.ext.autodoc',       # Automatically extract module docstrings
    'sphinx.ext.napoleon',      # Parse Google/NumPy style docstring architectures
    'sphinx.ext.viewcode',      # Add direct hyperlinked source code highlights
    'myst_parser'               # Native Markdown integration engine support
]

# Theme configuration options
html_theme = 'sphinx_rtd_theme'

# General settings
templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
