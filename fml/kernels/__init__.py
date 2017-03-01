# MIT License
#
# Copyright (c) 2016 Anders Steen Christensen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
FML kernels
=====
Provides
  1. gaussian_kernel
  2. laplacian_kernel
"""

from kernels import laplacian_kernel, gaussian_kernel, get_alpha, get_prediction, get_alpha_from_distance, get_prediction_from_distance, manhattan_distance, get_atomic_kernels_arad, get_atomic_kernels_gaussian, get_atomic_kernels_laplacian, get_atomic_symmetric_kernels_arad
from farad_kernels import fget_alpha_arad
from farad_kernels import fget_kernel_arad
from farad_kernels import fget_kernels_arad
from farad_kernels import fget_atomic_kernels_arad
from farad_kernels import fget_atomic_distance_arad
from force_kernels import force_kernel
from force_kernels import pyforce_kernel
__all__ = ['kernels']
