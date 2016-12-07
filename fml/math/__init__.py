"""
FML kernels
=====
Provides
  1. cho_solve
"""

from .cho_solve import cho_solve, cho_invert
from .distance import manhattan_distance, l2_distance, p_distance

__all__ = ['cho_solve', "distance"]
