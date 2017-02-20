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
from farad import atomic_arad_l2_distance_all

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


def get_l2_distance_arad(X1, X2, Z1, Z2, \
        width=0.2, cut_distance=6.0, r_width=1.0, c_width=0.5):
    """ Calculates the Gaussian distance matrix D for atomic ARAD for two
        sets of molecules

        K is calculated using an OpenMP parallel Fortran routine.

        Arguments:
        ==============
        X1 -- np.array of ARAD descriptors for molecules in set 1.
        X2 -- np.array of ARAD descriptors for molecules in set 2.
        Z1 -- List of lists of nuclear charges for molecules in set 1.
        Z2 -- List of lists of nuclear charges for molecules in set 2.

        Keyword arguments:
        width --
        cut_distance --
        r_width --
        c_width --

        Returns:
        ==============
        D -- The distance matrices for each sigma (4D-array, Nmol1 x Nmol2 x Natom1 x Natoms2)
    """

    print X1.shape
    print X2.shape

    amax = X1.shape[1]

    assert X1.shape[3] == amax, "ERROR: Check ARAD decriptor sizes! code = 1"
    assert X2.shape[1] == amax, "ERROR: Check ARAD decriptor sizes! code = 2"
    assert X2.shape[3] == amax, "ERROR: Check ARAD decriptor sizes! code = 3"

    nm1 = len(Z1)
    nm2 = len(Z2)

    assert X1.shape[0] == nm1,  "ERROR: Check ARAD decriptor sizes! code = 4"
    assert X2.shape[0] == nm2,  "ERROR: Check ARAD decriptor sizes! code = 5"

    N1 = []
    for Z in Z1:
        N1.append(len(Z))

    N2 = []
    for Z in Z2:
        N2.append(len(Z))

    N1 = np.array(N1,dtype=np.int32)
    N2 = np.array(N2,dtype=np.int32)

    nsigmas = len(sigmas)
    
    c1 = []
    for charges in Z1:
        c1.append(np.array([PTP[int(q)] for q in charges], dtype=np.int32))

    Z1_arad = np.zeros((nm1,amax,2))

    for i in range(nm1):
        for j, z in enumerate(c1[i]):
            Z1_arad[i,j] = z

    c2 = []
    for charges in Z2:
        c2.append(np.array([PTP[int(q)] for q in charges], dtype=np.int32))

    Z2_arad = np.zeros((nm2,amax,2))

    for i in range(nm2):
        for j, z in enumerate(c2[i]):
            Z2_arad[i,j] = z

    return atomic_arad_l2_distance_all(X1, X2, Z1_arad, Z2_arad, N1, N2, \
                nm1, nm2, width, cut_distance, r_width, c_width)
