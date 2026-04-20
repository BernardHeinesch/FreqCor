# Configuration file for the Sphinx documentation builder.

import os
import sys
import glob
import re

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


def _latest_example_output(pattern):
    """Return the newest matching example output filename (basename).

    Selection strategy:
    - Prefer parsing the run timestamp suffix '__YYMMDDTHHMM' if present.
    - Fallback to file modification time.
    """
    if not os.path.isdir(_EXAMPLES_OUTPUT):
        return ''
    matches = glob.glob(os.path.join(_EXAMPLES_OUTPUT, pattern))
    if not matches:
        return ''

    ts_re = re.compile(r'__(\d{6}T\d{4})\.')

    def _sort_key(p):
        base = os.path.basename(p)
        m = ts_re.search(base)
        ts = m.group(1) if m else ''
        mtime = 0.0
        try:
            mtime = os.path.getmtime(p)
        except Exception:
            pass
        return (ts, mtime)

    return os.path.basename(sorted(matches, key=_sort_key)[-1])


myst_substitutions = {
    'ex_2_all_individual_co2': _latest_example_output('2_all_individual_co2__BE-Lon__co2__cosp__all__*.png'),
    'ex_3_filtering_cof_gas': _latest_example_output('3_filtering_cof_gas__BE-Lon__co2__cosp__all__*.png'),
    'ex_3_filtering_cf_h_unst': _latest_example_output('3_filtering_CF_H_unst__BE-Lon__co2__cosp__all__*.png'),
    'ex_4_mean_tf_all_classes': _latest_example_output('4_mean_TF_all_classes__wd1__BE-Lon__co2__cosp__all__*.png'),
    'ex_4_mean_cosp_all_classes': _latest_example_output('4_mean_cosp_all_classes__wd1__BE-Lon__co2__cosp__all__*.png'),
    'ex_5_cf_vs_ws_st': _latest_example_output('5_CF_vs_ws__st__BE-Lon__co2__cosp__all__*.png'),
    'ex_5_cf_vs_ws_unst': _latest_example_output('5_CF_vs_ws__unst__BE-Lon__co2__cosp__all__*.png'),
}

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
