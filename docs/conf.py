# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import os, sys, shutil, subprocess
import re
from pathlib import Path

DIR = Path(__file__).parent.resolve()


# -- Project information -----------------------------------------------------

project = "SDTFE"
copyright = "2022, Steve Rangel, et al."
author = "Steve Rangel, et al."


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc.typehints",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "myst_parser",
]

autosectionlabel_prefix_document = True

autodoc_typehints = "description"
add_module_names = False

autosummary_generate = False
napoleon_numpy_docstring = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "README.md"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_theme_options = {"prev_next_buttons_location": None}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


def prepare(app):
    with open(DIR.parent / "README.md") as f:
        contents = f.read()

    # Filter out section titles for index.rst for LaTeX
    if app.builder.name == "latex":
        contents = re.sub(r"^(.*)\n[-~]{3,}$", r"**\1**", contents, flags=re.MULTILINE)

    with open(DIR / "README.md", "w") as f:
        f.write(contents)


def clean_up(app, exception):
    (DIR / "README.md").unlink()


def setup(app):
    app.add_css_file("css/custom.css")
    # Copy the readme in
    app.connect("builder-inited", prepare)

    # Clean up the generated readme
    app.connect("build-finished", clean_up)
