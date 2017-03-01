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

import numpy as np
import copy

from arad import ARAD

from data import NUCLEAR_CHARGE

from representations import fgenerate_coulomb_matrix
from representations import fgenerate_local_coulomb_matrix
from representations import fgenerate_atomic_coulomb_matrix

HOF_DFTB3 = dict()
HOF_DFTB3["H"] = -172.3145
HOF_DFTB3["C"] = -906.4342
HOF_DFTB3["N"] = -1327.2991
HOF_DFTB3["O"] = -1936.6161
HOF_DFTB3["S"] = -1453.3907


class Molecule:

    def __init__(self):
        self.natoms = -1
        self.energy = float("nan")
        self.molid = -1
        self.dftb3_energy = float("nan")
        self.dftb3_hof = float("nan")

        self.atomtypes = []
        self.nuclear_charges = []
        self.coordinates = []
        self.active_atoms = []

        # Container for misc properties
        self.properties = []
        self.properties2 = []

    def generate_coulomb_matrix(self, size=23):
        self.coulomb_matrix = fgenerate_coulomb_matrix(self.nuclear_charges, \
                self.coordinates, self.natoms, size)

    def generate_local_coulomb_matrix(self, calc="all",size=23):
        self.local_coulomb_matrix = fgenerate_local_coulomb_matrix( \
                self.nuclear_charges, self.coordinates, self.natoms, size)

    def generate_atomic_coulomb_matrix(self, calc="all",size=23):
        self.atomic_coulomb_matrix = fgenerate_atomic_coulomb_matrix( \
                self.nuclear_charges, self.coordinates, self.natoms, size)

    def generate_arad_descriptor(self, size=23):
        arad_object = ARAD(maxMolSize=size,maxAts=size)
        self.arad_descriptor = arad_object.describe(np.array(self.coordinates), \
                np.array(self.nuclear_charges))

        assert (self.arad_descriptor).shape[0] == size, "ERROR: Check ARAD descriptor size!"
        assert (self.arad_descriptor).shape[2] == size, "ERROR: Check ARAD descriptor size!"


    def read_xyz(self, filename):

        f = open(filename, "r")
        lines = f.readlines()
        f.close()

        self.natoms = int(lines[0])

        for line in lines[2:]:
            tokens = line.split()

            if len(tokens) < 4:
                break

            self.atomtypes.append(tokens[0])
            self.nuclear_charges.append(NUCLEAR_CHARGE[tokens[0]])

            x = float(tokens[1])
            y = float(tokens[2])
            z = float(tokens[3])

            self.coordinates.append(np.array([x, y, z]))

        self.coordinates = np.array(self.coordinates)



def get_lines(filename):

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    return lines

def parse_molecules(filename):

    lines = get_lines(filename)

    mols = []

    mol = Molecule()

    for line in lines:

        tokens = line.split()

        if len(tokens) == 1:

            if mol.natoms > 0:
                mols.append(mol)

            mol = Molecule()
            mol.natoms = int(tokens[0])

        if len(tokens) == 2:
            mol.molid = int(tokens[0])
            mol.energy = float(tokens[1])
            mol.dftb3_energy = parse_dft3_energy(mol.molid)


        if len(tokens) == 7:

            atom_type = tokens[0]
            mol.atomtypes.append(atom_type)
            mol.nuclear_charges.append(NUCLEAR_CHARGE[atom_type])
            x = float(tokens[4])
            y = float(tokens[5])
            z = float(tokens[6])

            mol.coordinates.append(np.array([x, y, z]))

            mol.dftb3_hof = 0.0
            mol.dftb3_hof += mol.dftb3_energy

            for atom in ["H", "C", "N", "O", "S"]:

                n = mol.atomtypes.count(atom)
                mol.dftb3_hof -= n * HOF_DFTB3[atom]

    # for mol in mols:
    #     print mol.molid, mol.energy, mol.dftb3_hof

    return mols

def parse_dft3_energy(molid):

    filename = "../logfiles/" + str(molid) + ".log"
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    energy = float("nan")
    for line in lines:
        if "Total Energy" in line:
            tokens = line.split()
            energy = float(tokens[2]) * 627.51

    return energy
