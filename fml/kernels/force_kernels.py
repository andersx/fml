# MIT License
#
# Copyright (c) 2017 Anders Steen Christensen, Felix Faber
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
import fml
from fforce_kernels import fsingle_force_kernel
from fforce_kernels import fcovariant_force_kernel

PTP = {\
         1  :[1,1] ,2:  [1,8]#Row1
       
        ,3  :[2,1] ,4:  [2,2]#Row2\
        #           C          N          O
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


ATOMTYPE = {
    "H" : 1,
    "C" : 6,
    "N" : 7,
    "O" : 8,
    "S" : 16
}



def periodic_distance(ai, aj, r_width, c_width):

    zi = PTP[ATOMTYPE[ai]]
    zj = PTP[ATOMTYPE[aj]]

    # Row-distance
    dr = np.exp(- abs(zi[0] - zj[0])/r_width)
    
    # Column-distance
    dc = np.exp(- abs(zi[1] - zj[1])/c_width)

    return dr * dc


def periodic_distance_table(r_width, c_width):

    n = len(PTP)
    n = 10

    D = np.zeros((n,n))

    for ai in range(n):
        for aj in range(n):

            zi = PTP[ai+1]
            zj = PTP[aj+1]
            
            # Row-distance
            dr = np.exp(- abs(zi[0] - zj[0])/r_width)
    
            # Column-distance
            dc = np.exp(- abs(zi[1] - zj[1])/c_width)

            D[ai,aj] = dr * dc

    return D


def pyforce_kernel(moli, molj, sigma_space, r_width=1.0, c_width=0.5):

    K = np.zeros((moli.natoms, molj.natoms, 3, 3))

    for i, xi in enumerate(moli.coordinates):
        for j, xj in enumerate(molj.coordinates):

            for ii, xii in enumerate(moli.coordinates):
                for jj, xjj in enumerate(molj.coordinates):

                    if (i == ii) or (j == jj):
                        continue
                    
                    ri = xi - xii
                    rj = xj - xjj

                    ri_norm = np.linalg.norm(ri)
                    rj_norm = np.linalg.norm(rj)

                    alpha_ij = (ri_norm**2 + rj_norm**2) / (4 * sigma_space**2)
                    gamma_ij = (ri_norm * rj_norm) / (2 * sigma_space**2)

                    phi_ij = np.exp(-alpha_ij) / (gamma_ij**2) * \
                            (gamma_ij * np.cosh(gamma_ij) - np.sinh(gamma_ij))

                    K[i,j] += (np.outer(ri, rj) * phi_ij) * periodic_distance(moli.atomtypes[ii], molj.atomtypes[jj], r_width, c_width)


                    # print phi_ij, ri_norm, rj_norm
                    # K[i,j] += (phi_ij)

            # print i, j, periodic_distance(moli.atomtypes[i], molj.atomtypes[j], r_width, c_width)

            K[i,j] *= periodic_distance(moli.atomtypes[i], molj.atomtypes[j], r_width, c_width)

    L = 1.0 / (2.0 * np.sqrt(np.pi * sigma_space**2))**3
    K *= L

    return K

def single_force_kernel(moli, molj, sigma_space, r_width=1.0, c_width=0.5):

    xi = moli.coordinates
    xj = molj.coordinates

    xi = np.swapaxes(xi,0,1)
    xj = np.swapaxes(xj,0,1)

    qi = np.array([PTP[ATOMTYPE[ai]] for ai in moli.atomtypes], dtype=np.int32)
    qj = np.array([PTP[ATOMTYPE[aj]] for aj in molj.atomtypes], dtype=np.int32)

    K = fsingle_force_kernel(xi, xj, moli.natoms, molj.natoms, qi, qj, sigma_space, r_width, c_width)

    K = np.swapaxes(K,0,3)
    K = np.swapaxes(K,1,2)

    return K


def force_kernel_all(moli, molj, sigma_space, r_width=1.0, c_width=0.5):
    """ Calculates the covariant kernel matrix.

        Arguments:
        ==============
        moli -- list of fml.Molecules()
        molj -- list of fml.Molecules()
        sigma -- The vaule of sigma for which to calculate the Kernel matrices.

        Keyword Argugments:
        ==============
        r_width = sigma for the row-width in the periodic table.
        c_width = sigma for the column-width in the periodic table.

        Returns:
        ==============
        K -- The kernel matrices (dimensions Nmoli x Nmoli x amax x amax x 3 x 3)
             where Nmoli is the number of molecules in moli, amax is the number of atoms
             in the largest molecule, etc.
    """

    amax = np.amax([mol.natoms for mol in moli] + [mol.natoms for mol in molj])

    nmi = len(moli)
    nmj = len(molj)

    ci = np.zeros((nmi, amax, 3))
    for i, mol in enumerate(moli):
        ci[i,:mol.natoms,:3] += np.array(mol.coordinates)[:mol.natoms,:3]

    cj = np.zeros((nmj, amax, 3))
    for i, mol in enumerate(molj):
        cj[i,:mol.natoms,:3] += np.array(mol.coordinates)[:mol.natoms,:3]

    qi = np.zeros((nmi, amax, 2), dtype=np.int32)
    for i, mol in enumerate(moli):
        qi[i,:mol.natoms,:2] += np.array([PTP[ATOMTYPE[ai]] for ai in mol.atomtypes], dtype=np.int32)[:mol.natoms,:2]
    
    qj = np.zeros((nmj, amax, 2), dtype=np.int32)
    for i, mol in enumerate(molj):
        qj[i,:mol.natoms,:2] += np.array([PTP[ATOMTYPE[ai]] for ai in mol.atomtypes], dtype=np.int32)[:mol.natoms,:2]


    ni = np.array([mol.natoms for mol in moli], dtype=np.int32)
    nj = np.array([mol.natoms for mol in molj], dtype=np.int32)


    K = fforce_kernel(ci, cj, ni, nj, qi, qj, nmi, nmj, amax, sigma_space, r_width, c_width)

    K = np.swapaxes(K,0,5)
    K = np.swapaxes(K,1,4)
    K = np.swapaxes(K,2,3)

    return K

def covariant_force_kernel(moli, molj, sigma_space, r_width=1.0, c_width=0.5):
    """ Calculates the covariant kernel matrix.

        Arguments:
        ==============
        moli -- list of fml.Molecules()
        molj -- list of fml.Molecules()
        sigma -- The vaule of sigma for which to calculate the Kernel matrices.

        Keyword Argugments:
        ==============
        r_width = sigma for the row-width in the periodic table.
        c_width = sigma for the column-width in the periodic table.

        Returns:
        ==============
        K -- The kernel matrices (dimensions Nmoli x Nmoli x amax x amax x 3 x 3)
             where Nmoli is the number of molecules in moli, amax is the number of atoms
             in the largest molecule, etc.
    """

    amax = np.amax([mol.natoms for mol in moli] + [mol.natoms for mol in molj])

    nacti = np.array([len(mol.active_atoms) for mol in moli] , dtype=np.int32)
    nactj = np.array([len(mol.active_atoms) for mol in molj] , dtype=np.int32)
    
    actmax = np.amax([len(mol.active_atoms) for mol in moli] + \
                     [len(mol.active_atoms) for mol in molj])

    acti = np.zeros((len(moli), actmax), dtype=np.int32)
    actj = np.zeros((len(molj), actmax), dtype=np.int32)

    for i, mol in enumerate(moli):
        for j, k in enumerate(mol.active_atoms):
            acti[i,j] = k + 1

    for i, mol in enumerate(molj):
        for j, k in enumerate(mol.active_atoms):
            actj[i,j] = k + 1

    nmi = len(moli)
    nmj = len(molj)

    # print nacti
    # print acti
    # print actmax

    # Make coordinate array
    ci = np.zeros((nmi, amax, 3))
    for i, mol in enumerate(moli):
        ci[i,:mol.natoms,:3] += np.array(mol.coordinates)[:mol.natoms,:3]

    # Make coordinate array
    cj = np.zeros((nmj, amax, 3))
    for i, mol in enumerate(molj):
        cj[i,:mol.natoms,:3] += np.array(mol.coordinates)[:mol.natoms,:3]

    # Make array of position in the periodic table
    qi = np.zeros((nmi, amax, 2), dtype=np.int32)
    for i, mol in enumerate(moli):
        qi[i,:mol.natoms,:2] += np.array([PTP[ATOMTYPE[ai]] for ai in mol.atomtypes], dtype=np.int32)[:mol.natoms,:2]
    
    # Make array of position in the periodic table
    qj = np.zeros((nmj, amax, 2), dtype=np.int32)
    for i, mol in enumerate(molj):
        qj[i,:mol.natoms,:2] += np.array([PTP[ATOMTYPE[ai]] for ai in mol.atomtypes], dtype=np.int32)[:mol.natoms,:2]

    
    # Make array of position in the periodic table
    zi = np.zeros((nmi, amax), dtype=np.int32)
    for i, mol in enumerate(moli):
        zi[i,:mol.natoms] += np.array([int(ATOMTYPE[ai]) for ai in mol.atomtypes], dtype=np.int32)[:mol.natoms]
    
    # Make array of position in the periodic table
    zj = np.zeros((nmj, amax), dtype=np.int32)
    for i, mol in enumerate(molj):
        zj[i,:mol.natoms] += np.array([int(ATOMTYPE[ai]) for ai in mol.atomtypes], dtype=np.int32)[:mol.natoms]


    # Make array of molecule sizes
    ni = np.array([mol.natoms for mol in moli], dtype=np.int32)

    # Make array of molecule sizes
    nj = np.array([mol.natoms for mol in molj], dtype=np.int32)


    pd = periodic_distance_table(r_width, c_width)

    K = fcovariant_force_kernel(ci, cj, ni, nj, qi, qj, zi, zj, nmi, nmj, amax,
            acti, actj, nacti, nactj, actmax, sigma_space, pd, r_width, c_width)

    K = np.swapaxes(K,0,5)
    K = np.swapaxes(K,1,4)
    K = np.swapaxes(K,2,3)

    return K

if __name__ == "__main__":

    import sys

    moli = fml.Molecule()
    molj = fml.Molecule()

    moli.read_xyz(sys.argv[1])
    molj.read_xyz(sys.argv[2])

    K = pyforce_kernel(moli, molj, 0.5)


    # # L = fsingle_force_kernel(moli.coordinates, molj.coordinates, moli.natoms, molj.natoms, 0.5, 0.6)

    L = single_force_kernel(moli, molj, 0.5)
    print L
    # print "K"
    # print K
   

    mols = [moli, molj]
    lol = force_kernel(mols, mols, 0.5)
    # print "lol"
    L2 = lol[0,1]
    print L.shape
    print L2.shape
    print K - L2[:3,:,:,:]
    print np.amax(np.abs(K - L2[:3,:,:,:]))
