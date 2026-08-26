project = "PackMem2"
copyright = "2026, Patrick Fuchs"
author = "Maya Zygadlo, Pierre Poulain, Romain Gautier, Patrick Fuchs"

import packmem2
# To uncomment when we have a release
#release = packmem2.__version__

extensions = ["sphinx.ext.mathjax",
              "sphinx.ext.autodoc",
              "sphinx.ext.napoleon",
              "nbsphinx",
              "myst_parser"
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

html_theme = "sphinx_book_theme"
# To uncomment when we have a logo.
#html_logo = "img/buildH_logo_small.png"
