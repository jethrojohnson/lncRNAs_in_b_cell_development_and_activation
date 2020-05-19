# -*- coding: utf-8 -*-
#
# Test documentation build configuration file, created by
# sphinx-quickstart on Mon Mar 23 15:27:57 2009.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys, os, itertools, re

################################################################
## CGAT pipelines - values to change

# path were documentation source resides
# Use environment variable SPHINX_DOCSDIR. If unset, the default
# is a level up in the src directory.
docsdir=os.environ.get( "SPHINX_DOCSDIR", "/ifs/devel/jethro/proj010/pipeline_docs/pipeline_proj010_lncrna" )

if not os.path.exists(docsdir):
    raise ValueError( "documentation directory '%s' not found" % docsdir )

# General information about the project.
project = u'CGAT - LncRNA Analysis Pipeline'
copyright = u'2011, Andreas Heger'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '1.0'
# The full version, including alpha/beta/rc tags.
release = '1'

################################################################
################################################################
################################################################
## The pipeline assumes that sphinxreport is called within the
## working directory. If the report is in a separate build directory,
## change the paths below.
## 
## directory with export directory from pipeline
## This should be a directory in the build directory - you can
## link from here to a directory outside the build tree, though.
exportdir = os.path.abspath('export')

datadir = os.path.abspath(".")

################################################################
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path = [os.path.abspath('.'), 
            os.path.abspath('%s/trackers' % docsdir)]  + sys.path

sphinxreport_cachedir=os.path.abspath("_cache")

# add warnings into the document
sphinxreport_show_errors = True

# database backend. Possible values are mysql, psql and sqlite
sphinxreport_sql_backend="sqlite:///%s/csvdb" % datadir

sphinxreport_urls=("data", "rst", "code")

# http port to use for serve.py
http_port=8080

# derived sets. This is a dictionary associating
# a slice with a track
trackers_derived_slices = {} 

# default slices to plot
trackers_default_slices = ["all"]

# images to plot
# Note that adding pdfs slows down the pipeline for complicated
# plots and should only be created as a last step and when
# latex is to be generated
renderer_images = None

# static images to create for each plot                                                                                                                                                                                                      
# a tuple of (id, format, dpi).
sphinxreport_images=( ( "hires", "hires.png", 200),
                      ( "eps", "eps", 50 ), )

# -- General configuration -----------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.doctest', 'sphinx.ext.coverage', 
              'sphinx.ext.pngmath', 'sphinx.ext.ifconfig',
              'SphinxReport.only_directives', 
              'SphinxReport.report_directive', 
              'sphinx.ext.inheritance_diagram',
              'SphinxReport.errors_directive',
              'SphinxReport.roles' ]

# Add any paths that contain templates here, relative to this directory.
templates_path = [os.path.relpath( '%s/_templates' % docsdir )]

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8'

# The master toctree document.
master_doc = 'contents'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of documents that shouldn't be included in the build.
#unused_docs = []

# List of directories, relative to source directory, that shouldn't be searched
# for source files.
exclude_trees = ['_build']

# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  Major themes that come with
# Sphinx are currently 'default' and 'sphinxdoc'.
html_theme = 'cgat'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation
#html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path=[os.path.join(os.path.dirname(docsdir), "themes")]

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = os.path.join( docsdir, "_templates", "cgat_logo.png" )

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
html_additional_pages = {'index': 'index.html', 'gallery':'gallery.html'}

# If false, no module index is generated.
#html_use_modindex = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
#html_show_sourcelink = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# If nonempty, this is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = ''

# Output file base name for HTML help builder.
htmlhelp_basename = 'Testdoc'


# -- Options for LaTeX output --------------------------------------------------

# The paper size ('letter' or 'a4').
#latex_paper_size = 'letter'

# The font size ('10pt', '11pt' or '12pt').
#latex_font_size = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('contents', 'Test.tex', ur'Test Documentation',
   ur'Andreas Heger', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# Additional stuff for the LaTeX preamble.
latex_preamble = """
   \usepackage{amsmath}
   \usepackage{amsfonts}
   \usepackage{amssymb}
   \usepackage{txfonts}
"""

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_use_modindex = True
