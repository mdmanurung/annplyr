from __future__ import annotations

import shutil
import sys
from datetime import datetime
from importlib.metadata import PackageNotFoundError, metadata
from pathlib import Path

from sphinxcontrib import katex

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "extensions"))

try:
    info = metadata("annplyr")
    project = info["Name"]
    author = info["Author"]
    version = info["Version"]
    urls = dict(project_url.split(", ") for project_url in info.get_all("Project-URL"))
    repository_url = urls["Source"]
except PackageNotFoundError:
    project = "annplyr"
    author = "annplyr developers"
    version = "0.3.0"
    repository_url = "https://github.com/mdmanurung/annplyr"

copyright = f"{datetime.now():%Y}, {author}."
release = version

bibtex_bibfiles = ["references.bib"]
templates_path = ["_templates"]
nitpicky = True
needs_sphinx = "4.0"

html_context = {
    "display_github": True,
    "github_user": "mdmanurung",
    "github_repo": "annplyr",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.katex",
    "sphinx_autodoc_typehints",
    "sphinx_design",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinxext.opengraph",
    *[path.stem for path in (HERE / "extensions").glob("*.py")],
]

autosummary_generate = True
autodoc_member_order = "groupwise"
default_role = "literal"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True
napoleon_custom_sections = [("Ownership", "returns_style")]
myst_heading_anchors = 6
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True
typehints_defaults = "braces"
always_use_bars_union = True

doctest_global_setup = """
import numpy as np
import pandas as pd
from anndata import AnnData
import annplyr as ap

adata = AnnData(
    X=np.array(
        [
            [2.0, 0.0, 1.0],
            [0.0, 1.0, 3.0],
            [4.0, 2.0, 0.0],
            [1.0, 3.0, 2.0],
        ],
        dtype=np.float32,
    ),
    obs=pd.DataFrame(
        {
            "batch": pd.Categorical(["A", "B", "A", "B"], categories=["A", "B", "unused"]),
            "donor_id": ["donor_1", "donor_2", "donor_1", "donor_2"],
            "cell_type": ["B cell", "T cell", "B cell", "T cell"],
            "n_counts": pd.array([1200, 900, 1500, 1100], dtype="Int64"),
            "pct_counts_mt": np.array([4.0, 12.0, 3.0, 7.5], dtype=np.float32),
        },
        index=["cell_0", "cell_1", "cell_2", "cell_3"],
    ),
    var=pd.DataFrame(index=["MS4A1", "CD79A", "NKG7"]),
)
adata.layers["counts"] = adata.X.copy()
adata.raw = adata
"""

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "anndata": ("https://anndata.scverse.org/en/stable/", None),
    "narwhals": ("https://narwhals-dev.github.io/narwhals/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "development/v0.3-executor-handoff.md",
    "development/v0.3-revised-plan.md",
    "development/v0.3-task-list.md",
]

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_title = project
html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
    "path_to_docs": "docs/",
    "navigation_with_keys": False,
}

pygments_style = "default"
katex_prerender = shutil.which(katex.NODEJS_BINARY) is not None
nitpick_ignore = [
    ("py:class", "annplyr._expr.AnnplyrSelector"),
    ("py:class", "annplyr._expr.Desc"),
]
