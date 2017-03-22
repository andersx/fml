#!/usr/bin/env python2

import sys
sys.path.append('../fml/')

import os
import numpy as np
import fml

import random

t_width = np.pi / 4.0 # 0.7853981633974483
d_width = 0.2
cut_distance = 6.0
r_width = 1.0
c_width = 0.5

PTP = {\
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


def periodic_distance(a, b):

    ra = PTP[int(a)][0]
    rb = PTP[int(b)][0]
    ca = PTP[int(a)][1]
    cb = PTP[int(b)][1]

    return (r_width**2 + (ra - rb)**2) * (c_width**2 + (ca - cb)**2)

def ksi(x):
    return 1.0 - np.sin(np.pi * x/(2.0 * cut_distance))

def aras_scalar(atom1, atom2, q1, q2, n1, n2, qall1, qall2):

    ksi1 = ksi(atom1[0,:n1])
    ksi2 = ksi(atom2[0,:n2])

    a = 1.0 / (np.sum(ksi1) * np.sum(ksi2)) \
            * np.sqrt(np.pi) * d_width * r_width**4 * c_width**4 \
            / periodic_distance(q1, q2)

    b = 0.0

    for i in range(n1):
        c = 0.0

        for j in range(n2):
            d = 0.0

            for k in range(n1):
                e = 0.0

                for l in range(n2):
                    e += ksi2[l] * np.exp(-(atom1[i,k] - atom2[k,l])**2 / (4.0 * t_width**2))

                d += e * ksi1[k]

            c += d * np.exp(-(atom1[0,i] - atom2[0,j])**2 / (4.0 * d_width**2)) * ksi2[j] \
                    / periodic_distance(qall1[i], qall2[j])

        b += c * ksi1[i]


    return a * b

def aras_distance(mol1, mol2, max_size=30):


    d = np.zeros((mol1.natoms, mol2.natoms))
    aa = np.zeros((mol1.natoms))
    bb = np.zeros((mol2.natoms))

    for i in range(mol1.natoms):
        atom1 = mol1.aras_descriptor[i]
        aa[i] = aras_scalar(atom1, atom1, 
            mol1.nuclear_charges[i],
            mol1.nuclear_charges[i], 
            mol1.natoms, mol1.natoms,
            mol1.nuclear_charges, mol1.nuclear_charges)

    for i in range(mol2.natoms):
        atom2 = mol2.aras_descriptor[i]
        bb[i] = aras_scalar(atom2, atom2, 
            mol2.nuclear_charges[i],
            mol2.nuclear_charges[i], 
            mol2.natoms, mol2.natoms,
            mol2.nuclear_charges, mol2.nuclear_charges)

    for i in range(mol1.natoms):
        atom1 = mol1.aras_descriptor[i]
        for j in range(mol2.natoms):
            atom2 = mol1.aras_descriptor[j]

            ab = aras_scalar(atom1, atom2, 
                mol1.nuclear_charges[i],
                mol2.nuclear_charges[j],
                mol1.natoms, mol2.natoms,
                mol1.nuclear_charges, mol2.nuclear_charges)
            d[i,j] = aa[i] + bb[j] - 2.0 * ab 

    return d



def get_kernel(mols1, mols2, sigma, max_size=30):


    K = np.zeros((len(mols1), len(mols2)))

    for i, mol1 in enumerate(mols1):
        for j, mol2 in enumerate(mols2):

            print i, j
            d = aras_distance(mol1, mol2)

            d *= -0.5 / (sigma**2)

            np.exp(d, d)

            K[i,j] = np.sum(d)

    return K



if __name__ == "__main__":
    
    mols = []
    path = "100mols/"
    filenames = os.listdir(path)


    print "Generating ARAS descriptors from FML interface ..."
    for filename in filenames:

        mol = fml.Molecule()
        mol.read_xyz(path + filename)
        mol.generate_aras_descriptor(size=30)
        mols.append(mol)


    random.shuffle(mols)

    train = mols[:40]
    test = mols[-40:]

    d = aras_distance(mols[0], mols[1])

    K = get_kernel(train, train, 100.0)
    Ks = get_kernel(mols, test, 100.0)
    
    np.save("K.npy", K)
    np.save("Ks.npy", Ks)

