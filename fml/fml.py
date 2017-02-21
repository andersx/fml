import numpy as np
import copy
from representations import fgenerate_coulomb_matrix
from representations import fgenerate_local_coulomb_matrix
from representations import fgenerate_atomic_coulomb_matrix

HOF_DFTB3 = dict()
HOF_DFTB3["H"] = -172.3145
HOF_DFTB3["C"] = -906.4342
HOF_DFTB3["N"] = -1327.2991
HOF_DFTB3["O"] = -1936.6161
HOF_DFTB3["S"] = -1453.3907

NUCLEAR_CHARGE = dict()
NUCLEAR_CHARGE["H"] = 1.0
NUCLEAR_CHARGE["C"] = 6.0
NUCLEAR_CHARGE["N"] = 7.0
NUCLEAR_CHARGE["O"] = 8.0
NUCLEAR_CHARGE["F"] = 9.0
NUCLEAR_CHARGE["Si"] = 14.0
NUCLEAR_CHARGE["S"] = 16.0
NUCLEAR_CHARGE["Cl"] = 17.0
NUCLEAR_CHARGE["Ge"] = 32.0

class ARAD(object):

    def __init__(self,maxMolSize = 30,maxAts = 30,cut = 4., debug=False):
        self.tag = 'coords'
        self.debug = debug
        self.maxMolSize = maxMolSize
        self.maxAts = maxAts
        self.cut = cut
        np.set_printoptions(threshold= 'nan')

        self.PTP = {\
            1  :[1,1] ,2:  [1,8]#Row1

            ,3  :[2,1] ,4:  [2,2]#Row2\
            ,5  :[2,3] ,6:  [2,4] ,7  :[2,5] ,8  :[2,6] ,9  :[2,7] ,10 :[2,8]\

            ,11 :[3,1] ,12: [3,2]#Row3\
            ,13 :[3,3] ,14: [3,4] ,15 :[3,5] ,16 :[3,6] ,17 :[3,7] ,18 :[3,8]\

            ,19 :[4,1] ,20: [4,2]#Row4\
            ,31 :[4,3] ,32: [4,4] ,33 :[4,5] ,34 :[4,6] ,35 :[4,7] ,36 :[4,8]\
            ,21 :[4,9] ,22: [4,10],23 :[4,11],24 :[4,12],25 :[4,13],26 :[4,14],27 :[4,15],28 :[4,16],29 :[4,17],30 :[4,18]\

            ,37 :[5,1] ,38: [5,2]#Row5\
            ,49 :[5,3] ,50: [5,4] ,51 :[5,5] ,52 :[5,6] ,53 :[5,7] ,54 :[5,8]\
            ,39 :[5,9] ,40: [5,10],41 :[5,11],42 :[5,12],43 :[5,13],44 :[5,14],45 :[5,15],46 :[5,16],47 :[5,17],48 :[5,18]\

            ,55 :[6,1] ,56: [6,2]#Row6\
            ,81 :[6,3] ,82: [6,4] ,83 :[6,5] ,84 :[6,6] ,85 :[6,7] ,86 :[6,8]
                ,72: [6,10],73 :[6,11],74 :[6,12],75 :[6,13],76 :[6,14],77 :[6,15],78 :[6,16],79 :[6,17],80 :[6,18]\
            ,57 :[6,19],58: [6,20],59 :[6,21],60 :[6,22],61 :[6,23],62 :[6,24],63 :[6,25],64 :[6,26],65 :[6,27],66 :[6,28],67 :[6,29],68 :[6,30],69 :[6,31],70 :[6,32],71 :[6,33]\

            ,87 :[7,1] ,88: [7,2]#Row7\
            ,113:[7,3] ,114:[7,4] ,115:[7,5] ,116:[7,6] ,117:[7,7] ,118:[7,8]\
                ,104:[7,10],105:[7,11],106:[7,12],107:[7,13],108:[7,14],109:[7,15],110:[7,16],111:[7,17],112:[7,18]\
            ,89 :[7,19],90: [7,20],91 :[7,21],92 :[7,22],93 :[7,23],94 :[7,24],95 :[7,25],96 :[7,26],97 :[7,27],98 :[7,28],99 :[7,29],100:[7,30],101:[7,31],101:[7,32],102:[7,14],103:[7,33]}

    def describe(self,coords,ocupationList,cell = None):
        L = len(coords)
        coords = np.asarray(coords)
        ocupationList = np.asarray(ocupationList)
        M =  np.zeros((self.maxMolSize,5,self.maxAts))

        if cell is not None:
            coords = dot(coords,cell)
            #print cell
            nExtend = (floor(self.cut/linalg.norm(cell,2,axis = 0)) + 1).astype(int)
            #print nExtend
            for i in range(-nExtend[0],nExtend[0] + 1):
                for j in range(-nExtend[1],nExtend[1] + 1):
                    for k in range(-nExtend[2],nExtend[2] + 1):

                        #print i, j, k
                        #print i*cell[0,:] + j*cell[1,:] + k*cell[2,:]
                        if i == -nExtend[0] and j  == -nExtend[1] and k  == -nExtend[2]:
                            coordsExt = coords + i*cell[0,:] + j*cell[1,:] + k*cell[2,:]
                            ocupationListExt = copy.copy(ocupationList)
                        else:
                            ocupationListExt = append(ocupationListExt,ocupationList)
                            coordsExt = append(coordsExt,coords + i*cell[0,:] + j*cell[1,:] + k*cell[2,:],axis = 0)

        else:
            coordsExt = np.copy(coords)
            ocupationListExt = np.copy(ocupationList)

        M[:,0,:] = 1E+100

        for i in range(L):
            #print '+'
            #Calculate Distance
            cD = - coords[i] + coordsExt[:]

            ocExt =  np.asarray([self.PTP[o] for o in  ocupationListExt])


            #Obtaining angles
            angs = np.sum(cD[:,np.newaxis] * cD[np.newaxis,:], axis = 2)
            D1 = np.sqrt(np.sum(cD**2, axis = 1))
            D2 = D1[:,np.newaxis]*D1[np.newaxis,:]  + 10**-12
            #print D2.shape
            angs = np.arccos(angs/D2)

            angs = np.nan_to_num(angs)


            #Obtaining cos and sine terms
            cosAngs = np.cos(angs) * (1. - np.sin(np.pi * D1[np.newaxis,:]/(2. * self.cut))) #cos(pi * D1[newaxis,:]/(2 * self.cut))
            sinAngs = np.sin(angs) * (1. - np.sin(np.pi * D1[np.newaxis,:]/(2. * self.cut))) #cos(pi * D1[newaxis,:]/(2 * self.cut))



            args = np.argsort(D1)

            D1 = D1[args]

            ocExt = np.asarray([ocExt[l] for l in args])


            #cosAngs = asarray([[cosAngs[k,l] for l in args] for k in args])
            #sinAngs = asarray([[sinAngs[k,l] for l in args] for k in args])


            cosAngs = cosAngs[args,:]
            cosAngs = cosAngs[:,args]
            sinAngs = sinAngs[args,:]
            sinAngs = sinAngs[:,args]

            args = np.where(D1 < self.cut)[0]
            #print args

            D1 = D1[args]
            #print D1#[:50]

            ocExt = np.asarray([ocExt[l] for l in args])

            #cosAngs = asarray([[cosAngs[k,l] for l in args] for k in args])
            #sinAngs = asarray([[sinAngs[k,l] for l in args] for k in args])

            cosAngs = cosAngs[args,:]
            cosAngs = cosAngs[:,args]
            sinAngs = sinAngs[args,:]
            sinAngs = sinAngs[:,args]


            #print all(cosAngs == cosAngsVali)
            #print all(sinAngs == sinAngsVali)
            D1 = D1[1:]
            ocExt = ocExt[1:]
            cosAngs = cosAngs[1:,1:]
            sinAngs = sinAngs[1:,1:]

            '''
            print D1.shape
            print ocExt.shape
            print cosAngs.shape
            print sinAngs.shape
            '''

            M[i,0,: len(D1)] = D1
            M[i,1,: len(D1)] = ocExt[:,0]
            M[i,2,: len(D1)] = ocExt[:,1]
            M[i,3,: len(D1)] = np.sum(cosAngs,axis = 1)
            M[i,4,: len(D1)] = np.sum(sinAngs,axis = 1)


        return M


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
