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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

#project = 'lammpstutorials'
copyright = 'All source code is available under the GNU General Public License v3.0'
author = 'Simon Gravelle'

# The full version, including alpha/beta/rc tags
#release = '0.1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinxcontrib.googleanalytics']


# google analytics
googleanalytics_id = 'G-W1WGEC5GQ8'


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo'
html_title = "    "

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = ["custom.css"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'tango' # 'murphy' # 

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
#html_css_files = [
#    'css/custom.css',
#]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# note: do not add .ipynb when nbspinx is enabled, 
# otherwise you get the "missing title" error
source_suffix = ['.rst']

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
#html_theme_options = {}
html_logo = "../../figures/logo.png"
#html_theme_options = {
#    "light_logo": "../static/logo_MAICOS_light.png",
 #   "dark_logo": "../static/logo_MAICOS_dark.png",
#}

# hide extra stuff
html_show_copyright = False
html_show_sphinx = False
#html_show_search_summary = False
html_short_title = "True"
