# Sphinx configuration for the omataxonomy documentation.
from importlib import metadata

project = "omataxonomy"
author = "Adrian Altenhoff"
copyright = "Adrian Altenhoff, DessimozLab"
release = metadata.version("omataxonomy")
version = ".".join(release.split(".")[:2])

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

myst_enable_extensions = ["colon_fence"]
source_suffix = {".md": "markdown"}
exclude_patterns = ["_build", "examples"]

autodoc_member_order = "bysource"
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

html_theme = "furo"
html_title = f"omataxonomy {release}"
