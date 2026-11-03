from __future__ import annotations

import shutil
import sys
from datetime import datetime
from importlib.metadata import PackageNotFoundError, metadata
from pathlib import Path

from sphinxcontrib import katex

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "extensions"))
sys.path.insert(0, str(HERE))

try:
    info = metadata("annplyr")
    project = info["Name"]
    author = info["Author"] or info["Author-email"].split(" <")[0]
    version = info["Version"]
    urls = dict(project_url.split(", ") for project_url in info.get_all("Project-URL"))
    repository_url = urls["Source"]
except PackageNotFoundError:
    project = "annplyr"
    author = "Mikhael Manurung"
    version = "0.5.0"
    repository_url = "https://github.com/mdmanurung/annplyr"

copyright = f"{datetime.now():%Y}, {author}."
release = version

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
    "sphinx.ext.napoleon",
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
nb_execution_mode = "cache"
nb_execution_raise_on_error = True
nb_execution_timeout = 900
nb_merge_streams = True
typehints_defaults = "braces"
always_use_bars_union = True

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
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
    # Maintainer notes stay in the repo; scripts and tests read some of them.
    "development/**",
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
    "icon_links": [
        {
            "name": "Changelog",
            "url": f"{repository_url}/blob/main/CHANGELOG.md",
            "icon": "fa-solid fa-clock-rotate-left",
            "type": "fontawesome",
        },
    ],
}

pygments_style = "default"
katex_prerender = shutil.which(katex.NODEJS_BINARY) is not None
nitpick_ignore = [
    ("py:class", "annplyr._expr.AnnplyrSelector"),
    ("py:class", "annplyr._expr.Desc"),
]
