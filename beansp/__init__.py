# -*- coding: utf-8 -*-

"""Top-level package for BEANSp."""

__author__ = """Adelle Goodwin and Duncan Galloway"""
__email__ = 'adelle.goodwin@curtin.edu.au'
__version__ = '2.51.4'

# this allows "from beans import *"
# __all__ = ["beans"]

# Simplify imports from the package.
#   enables use "from beansp import Beans"
#   rather than "from beansp.beans import Beans"
from .beans import Beans
# from .beans import *
