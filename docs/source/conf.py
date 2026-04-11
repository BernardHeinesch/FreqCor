# Configuration file for the Sphinx documentation builder.

import os
import sys

_HERE = os.path.abspath(os.path.dirname(__file__))
_SRC = os.path.abspath(os.path.join(_HERE, "..", "..", "src"))
sys.path.insert(0, _SRC)

# -- Project information -----------------------------------------------------

project = 'FREQCOR'
author = 'FREQCOR contributors'
release = '1.4'

# -- General configuration ---------------------------------------------------

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
]

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_attr_annotations = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
}

templates_path = ['_templates']
exclude_patterns = []

suppress_warnings = [
    'ref.python',
]

autosummary_generate = True

autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': False,
    'no-index': True,
    'exclude-members': '__weakref__',
}

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = [
    'custom.css',
]
_EXAMPLES_OUTPUT = os.path.abspath(os.path.join(_HERE, "..", "..", "examples", "output"))
html_extra_path = [_EXAMPLES_OUTPUT] if os.path.exists(_EXAMPLES_OUTPUT) else []

myst_enable_extensions = [
    'colon_fence',
    'deflist',
    'dollarmath',
    'fieldlist',
    'html_admonition',
    'html_image',
    'replacements',
    'smartquotes',
    'substitution',
    'tasklist',
]
