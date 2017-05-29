#!/usr/bin/env python2

import sys
# sys.path.append('../fml/')

import os
import numpy as np
import fml
from fml.kernels import get_atomic_kernels_aras

import random

t_width = np.pi / 4.0 # 0.7853981633974483
d_width = 0.2
cut_distance = 6.0
r_width = 1.0
c_width = 0.5

def aras_kernels(mols1, mols2):

    sigmas = [0.1 * 2**i for i in range(20)]

    amax = mols1[0].aras_descriptor.shape[0]

    nm1 = len(mols1)
    nm2 = len(mols2)

    X1 = np.array([mol.aras_descriptor for mol in mols1]).reshape((nm1,amax,3+amax*2,amax))
    X2 = np.array([mol.aras_descriptor for mol in mols2]).reshape((nm2,amax,3+amax*2,amax))

    Z1 = [mol.nuclear_charges for mol in mols1]
    Z2 = [mol.nuclear_charges for mol in mols2]

    print X1.shape
    print X2.shape

    K = get_atomic_kernels_aras(X1, X2, Z1, Z2, sigmas, \
    t_width=np.pi/2.0, d_width=0.2, cut_distance=5.0, r_width=1.0, c_width=1.0,
    order=1, scale_angular=1.0/10.0)

    return K



if __name__ == "__main__":
    
    mols = []
    path = "xyz/"
    filenames = os.listdir(path)

    np.set_printoptions(linewidth=99999999999999999)
    print "Generating ARAS descriptors from FML interface ..."
    # for filename in sorted(filenames[:2]):
    # for filename in sorted(["06.xyz"]):
    for filename in sorted(filenames):

        mol = fml.Molecule()
        mol.read_xyz(path + filename)
        print mol.atomtypes
        print filename
        mol.generate_aras_descriptor(size=11)
        mols.append(mol)

    print mols[0].aras_descriptor[0]
    K = aras_kernels(mols, mols)

    print K[8]


