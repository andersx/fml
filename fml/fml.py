import numpy as np
from representations import generate_coulomb_matrix
from representations import generate_local_coulomb_matrix

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
        self.coordinates = []

    def generate_coulomb_matrix(self):
        self.coulomb_matrix = generate_coulomb_matrix(self.atomtypes, self.coordinates)

    def generate_local_coulomb_matrix(self):
        self.local_coulomb_matrix = generate_local_coulomb_matrix(self.atomtypes, \
            self.coordinates, calc="all")

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

            mol.atomtypes.append(tokens[0])
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
