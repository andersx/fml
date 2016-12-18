"""
FML kernels
=====
Provides
  1. gaussian_kernel
  2. laplacian_kernel
"""

from kernels import laplacian_kernel, gaussian_kernel, get_alpha, get_prediction, get_alpha_from_distance, get_prediction_from_distance, manhattan_distance, get_atomic_kernels_arad, get_atomic_kernels_gaussian
from farad_kernels import fget_alpha_arad
from farad_kernels import fget_kernel_arad
from farad_kernels import fget_kernels_arad
__all__ = ['kernels']
