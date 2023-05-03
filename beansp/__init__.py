# -*- coding: utf-8 -*-

"""Top-level package for BEANSp."""

__author__ = """Adelle Goodwin"""
__email__ = 'adelle.goodwin@monash.edu'
### this just does not work at all - importing version from a dedicated file,
### that is not __init__.py
# from .version import __version__
__version__ = '0.9.2'

# this allows "from beans import *"
# __all__ = ["beans"]

# from .beans import *
from .beans import Beans
# from beansp.beans import Beans
