# -*- coding: utf-8 -*-
"""Sphinx configuration."""
import os
import sys

import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("../tests/"))
sys.path.insert(0, os.path.abspath("../scripts/"))

project = "BIC086"
author = "Ryan B Patterson-Cross"
copyright = "2021, IMS-MRL Bioinformatics and Biostatistic Core"
extensions = ["sphinx_rtd_theme", "sphinx.ext.autodoc", "sphinx.ext.napoleon"]

napoleon_google_docstrings = False
napoleon_numpy_docstrings = True
napoleon_use_param = False

html_theme = "sphinx_rtd_theme"
