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

from numpy import empty

from fdistance import fmanhattan_distance
from fdistance import fl2_distance
from fdistance import fp_distance_integer, fp_distance_double

def manhattan_distance(A, B):

    if len(A.shape) != 2 or len(B.shape) != 2:
        raise ValueError('expected matrices of dimension=2')

    if B.shape[0] != A.shape[0]:
        raise ValueError('expected matrices containing vectors of same size')

    na = A.shape[1]
    nb = B.shape[1]

    D = empty((na, nb), order='F')

    fmanhattan_distance(A, B, D)

    return D

def l2_distance(A, B):

    if len(A.shape) != 2 or len(B.shape) != 2:
        raise ValueError('expected matrices of dimension=2')

    if B.shape[0] != A.shape[0]:
        raise ValueError('expected matrices containing vectors of same size')

    na = A.shape[1]
    nb = B.shape[1]

    D = empty((na, nb), order='F')

    fl2_distance(A, B, D)

    return D

def p_distance(A, B, p=2):

    if len(A.shape) != 2 or len(B.shape) != 2:
        raise ValueError('expected matrices of dimension=2')

    if B.shape[0] != A.shape[0]:
        raise ValueError('expected matrices containing vectors of same size')

    na = A.shape[1]
    nb = B.shape[1]

    D = empty((na, nb), order='F')


    if (type(p) == type(1)):
        if (p == 2):
            fl2_distance(A, B, D)
        else:
            fp_distance_integer(A, B, D, p)

    elif (type(p) == type(1.0)):
        fp_distance_double(A, B, D, p)
    else:
        raise ValueError('expected exponent of integer or float type')

    return D
