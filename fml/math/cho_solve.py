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

from fcho_solve import fcho_solve
from numpy import zeros

def cho_solve(A, y):
    """ Solves [A x = y] for x using a Cholesky decomposition
        via calls to LAPACK dpotrf and dpotrs in the F2PY module.

        Arguments:
        ==============
        A -- the A-matrix (symmetric and positive definite).
        y -- the right-hand side of the equation (vector).

        Returns:
        ==============
        x -- the vector for with the equation has been solved.
    """

    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix')

    if len(y.shape) != 1 or y.shape[0] != A.shape[1]:
        raise ValueError('expected matrix and vector of same stride size')

    n = A.shape[0]

    x = zeros((n))
    ffcho_solve(A,y,x)

    return x
