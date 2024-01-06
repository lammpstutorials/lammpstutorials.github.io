# Configuration file for the Sphinx documentation builder.

copyright = 'All source code is available under the GNU General Public License v3.0'
author = 'Simon Gravelle'

extensions = ['sphinxcontrib.googleanalytics',
              'sphinx_togglebutton',
              'sphinx_favicon'] # 'rst2pdf.pdfbuilder'

googleanalytics_id = 'G-W1WGEC5GQ8'

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

#html_context = {
#    "lammpstutorials.github.io ": True, # Integrate GitHub
#    "github_user": "lammpstutorials", # Username
#    "github_repo": "lammpstutorials.github.io ", # Repo name
#    "github_version": "version2.0", # Version
#}
