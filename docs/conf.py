"""Sphinx configuration for fossil_sim."""

import shutil
from pathlib import Path

from fossil_sim import __version__

DOCS = Path(__file__).parent
ROOT = DOCS.parent

project = "fossil_sim"
author = "Luca Pagano"
copyright = "2026, Luca Pagano"
release = __version__
version = __version__

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
]

myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist"]
myst_heading_anchors = 3
nb_execution_mode = (
    "off"  # notebooks are committed executed (scripts/build_notebooks.py)
)

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "author_year"

autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "astropy": ("https://docs.astropy.org/en/stable", None),
}

exclude_patterns = ["_build"]
html_theme = "pydata_sphinx_theme"
html_title = f"fossil_sim {version}"
html_theme_options = {
    "github_url": "https://github.com/paganol/fossil_sim",
    "use_edit_page_button": False,
}

# Notebooks live in notebooks/ at the repository root; copy them into the
# documentation source tree (git-ignored) so that Sphinx can render them.
shutil.copytree(
    ROOT / "notebooks", DOCS / "tutorials" / "notebooks", dirs_exist_ok=True
)
