#!/usr/bin/env python2

import sys
sys.path.append('../fml/')

import os
import numpy as np

import fml

from fml.math.distance import get_l2_distance_arad
from fml.kernels import get_atomic_kernels_arad

if __name__ == "__main__":


    # Get list of xyzfilenames
    path = "../tests/xyz/"
    filenames = os.listdir(path)

    mols = []

    # Make list of fml.Molecule()
    for filename in filenames:

        # Initialize molecule
        mol = fml.Molecule()

        # Read coordinates etc from xyz file
        mol.read_xyz(path + filename)

        # Generate ARAD descriptor, size = max size of molecule
        mol.generate_arad_descriptor(size=20)

        # Add mols to list
        mols.append(mol)


    x1 = []
    z1 = []

    # Make list of molecules
    for mol in mols[:5]:

        x1.append(mol.arad_descriptor)
        z1.append(mol.nuclear_charges)

    x2 = []
    z2 = []

    # Make another list of molecules
    for mol in mols[5:]:

        x2.append(mol.arad_descriptor)
        z2.append(mol.nuclear_charges)

    # Convert descriptors to np.array
    x1 = np.array(x1)
    x2 = np.array(x2)

    # Get Distance Matrix from FML
    D = get_l2_distance_arad(x1,x2,z1,z2)

    print D.shape

    # Get Gaussian kernel matrix for ARAD descriptor (for a list of sigmas)
    sigmas = [0.1, 1.0, 10.0, 100.0]
    K = get_atomic_kernels_arad(x1, x2, z1, z2, sigmas)

    print K.shape
