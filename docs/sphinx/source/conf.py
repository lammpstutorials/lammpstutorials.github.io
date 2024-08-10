# Configuration file for the Sphinx documentation builder.

copyright = 'All source code is available under the GNU General Public License v3.0'
author = 'See the AUTHORS.rst'

extensions = ['sphinx_togglebutton',
              'sphinx_favicon',
              'sphinxcontrib.bibtex']

templates_path = ['_templates']

exclude_patterns = []

html_theme = 'furo'
html_title = "    "

html_static_path = ['_static']
html_css_files = ["custom.css"]

pygments_style = 'tango'

source_suffix = ['.rst']

html_logo = "_static/logo.png"

favicons = [
    {"href": "favicon-32x32.png"},
]

html_show_copyright = False
html_show_sphinx = False
html_short_title = "True"

bibtex_bibfiles = ['journal-article.bib']
