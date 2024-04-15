import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PSutils'
copyright = '2024, Jixing Zhong'
author = 'Jixing Zhong'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'numpydoc',
    'sphinx.ext.autodoc',
    "sphinx.ext.intersphinx",
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    "sphinx.ext.extlinks",
    'sphinx_autodoc_typehints',  # needs to be after napoleon
    # "sphinx.ext.linkcode",
    # "sphinx_design",
    # "sphinx_search.extension",
    'sphinx_rtd_theme'
]

# autosummary_generate = True
napoleon_google_docstring = False
napoleon_numpy_docstring = True
# napoleon_include_init_with_doc = False
# napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]

numpydoc_class_members_toctree = False

templates_path = ['_templates']
exclude_patterns = []

myst_url_schemes = ("http", "https")

# autodoc_typehints_format = "short"
intersphinx_mapping = dict(
    anndata=("https://anndata.readthedocs.io/en/stable/", None),
    matplotlib=("https://matplotlib.org/stable/", None),
    numpy=("https://numpy.org/doc/stable/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy/", None),
    seaborn=("https://seaborn.pydata.org/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    skimage=('https://scikit-image.org/docs/stable/', None),
    scanpy=('https://scanpy.readthedocs.io/en/stable/', None),
    pydeseq2=('https://pydeseq2.readthedocs.io/en/latest/', None),
    SpatialDE=('https://pmbio.github.io/SpatialDE/index.html', None),

)



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
# html_static_path = ['_static']
html_static_path = []
