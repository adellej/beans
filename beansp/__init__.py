# -*- coding: utf-8 -*-

"""Top-level package for BEANSp."""

__author__ = 'Adelle Goodwin'
__email__ = 'adelle.goodwin@monash.edu'
__version__ = '0.9.2'

# this allows "from beans import *"
# __all__ = ["beans"]

# Simplify imports from the package.
#   enables use "from beansp import Beans"
#   rather than "from beansp.beans import BEans"
from .beans import Beans
# from .beans import *
