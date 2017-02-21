#!/usr/bin/env python2

import sys
sys.path.append('../fml/')

import os
import numpy as np
import fml

from fml.math.distance import get_l2_distance_arad
from fml.kernels import get_atomic_kernels_arad

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
#@jit#(nopython = True, nogil = True)   


#@jit#(nopython = True, nogil = True)   


def condense(X1,X2,Z1,Z2,width = 0.2, cutDist = 6.,RWidth = 1 ,CWidth = 0.5,num_cores = 1, printTime = False,molecular = False):#RWidth =1./sqrt(2),Cwidth = 1./sqrt(2)):#,RWidth =1.,Cwidth = 1.):
    def createZINP(Z,l):
        Z_inp = -1*np.ones((len(Z),l,2)).astype(int)
        for i in range(len(Z)):
            Z_inp[i,:len(Z[i])] = np.asarray([PTP[z] for z in Z[i]])
           
        return Z_inp
   
    RWidth = float(RWidth)
    Cwidth = float(CWidth)
    Z1_inp = createZINP(Z1,len(X1[0]))
    Z2_inp = createZINP(Z2,len(X2[0]))
   
    from time import time

   

    start = time()
    DistMatrix = run(X1,X2,Z1_inp,Z2_inp,width,cutDist,RWidth,Cwidth,num_cores,molecular = molecular)
    if printTime:
        print 'Time Paralell: ', time() - start
    #print DistMatrix[0,0]
    #print DistMatrix.shape
    return DistMatrix


class KeyboardInterruptError(Exception): pass

def run(X1,X2,Z1,Z2,width,cutDist,RWidth,Cwidth,Ncores = 2,NIters = 20, molecular = True ):
    NIters = min(len(X2),NIters)
    from multiprocessing import Pool
    try:
       
        p = Pool(Ncores)
       
        w_list  = np.repeat(width,len(X1))
        cD_list = np.repeat(cutDist,len(X1))
        RW_list = np.repeat(RWidth,len(X1))
        CW_list = np.repeat(Cwidth,len(X1))
        const_list =np.repeat(np.NaN,len(X1))
       
        w_list_2  = np.repeat(width,len(X2))
        cD_list_2 = np.repeat(cutDist,len(X2))
        RW_list_2 = np.repeat(RWidth,len(X2))
        CW_list_2 = np.repeat(Cwidth,len(X2))
       
       
       
   
        args1 = [X1,Z1,w_list,cD_list,RW_list,CW_list]
        args2 = [X2,Z2,w_list_2,cD_list_2,RW_list_2,CW_list_2]   
       
        if molecular:
            pass_nonmNoDist
            D1 = np.asarray(p.map(pass_nonmNoDist,zip(*args1)))
            D2 = np.asarray(p.map(pass_nonmNoDist,zip(*args2)))
            D1 = unravelKerMatrix(D1)
            D2 = unravelKerMatrix(D2)
       
        else:
            D1 = np.asarray(p.map(pass_calcNorm,zip(*args1)))
            D2 = np.asarray(p.map(pass_calcNorm,zip(*args2)))
       
        i_2 = 0
        x2 = np.repeat(np.asarray([X2[i_2:min(len(X2),i_2 + NIters)]]),len(X1),axis = 0)
        z2 = np.repeat(np.asarray([Z2[i_2:min(len(X2),i_2 + NIters)]]),len(X1),axis = 0)
        d2 = np.repeat(np.asarray([D2[i_2:min(len(X2),i_2 + NIters)]]),len(X1),axis = 0)
        d2 = np.repeat(np.asarray([D2[i_2:min(len(X2),i_2 + NIters)]]),len(X1),axis = 0)
        if molecular:
            args = [X1,x2,Z1,z2,w_list,cD_list,RW_list,CW_list]
            results = np.asarray(p.map(pass_calcProdnoDist,zip(*args)))   
        else:
            args = [X1,x2,Z1,z2,D1,d2,w_list,cD_list,RW_list,CW_list,const_list]
            results = np.asarray(p.map(pass_calcProd,zip(*args)))

        for i_2 in range(NIters,len(X2),NIters):
            x2 = np.repeat(np.asarray([X2[i_2:min(len(X2),i_2 + NIters)]]),len(X1),axis = 0)
            z2 = np.repeat(np.asarray([Z2[i_2:min(len(X2),i_2 + NIters)]]),len(X1),axis = 0)
            d2 = np.repeat(np.asarray([D2[i_2:min(len(X2),i_2 + NIters)]]),len(X1),axis = 0)
           
            if molecular:
                args = [X1,x2,Z1,z2,w_list,cD_list,RW_list,CW_list]
                results = append(results,
                    np.asarray(p.map(pass_calcProdnoDist,zip(*args))),
                    axis = 1)   
            else:
                args = [X1,x2,Z1,z2,D1,d2,w_list,cD_list,RW_list,CW_list,const_list]
                results = append(results,
                    np.asarray(p.map(pass_calcProd,zip(*args))),
                    axis = 1)
        p.close()
        p.join()
            #print results.shape
       
        if molecular:
            results = unravelKerMatrix(results)
            #print results[0:5,0:5]
            #print D1[0:5]
            #print D2[0:5]           
            #print results.shape,D1[:,newaxis].shape,D2[newaxis].shape

            D = D1[:,newaxis] +  D2[newaxis]
            results =  D - 2 * results
            #print results[0:5,0:5]
       
        return results   
    except KeyboardInterrupt:
        print 'got ^C'
        p.terminate()
        p.join()
        raise KeyboardInterrupt
   
    except Exception, inst:
        p.terminate()
        p.join()
        print 'got exception: %r' % (inst,)
        print type(inst)     # the exception instance
        print inst.args      # arguments stored in .args
        print inst           # __str__ allows args to be printed directly
        raise inst
   
def pass_calcNorm(args):
    try:
        return calcNorm(*args)
    except KeyboardInterrupt:
        raise KeyboardInterruptError()

def calcNorm(X,Z,width,cutDist,RWidth,Cwidth):
   
    DistMatrix  = np.zeros(len(X))
    for i in  range(len(X)):
        if Z[i,0] == -1:
            break
        dd = _MDist(X[i],X[i],width,cutDist,RWidth,Cwidth)
        DistMatrix[i] = dd

    return DistMatrix   

#@jit
def pass_calcProd(args):
    try:   
        return calcProd(*args)
    except KeyboardInterrupt:
        raise KeyboardInterruptError()
#@jit#(nopython = True, nogil = True)               
def calcProd(X1,X2,Z1,Z2,D1,D2,width,cutDist,RWidth,Cwidth,const):
    #maxRDist = 8*RWidth
    #maxCDist = 8*Cwidth
   

    DistMatrix  = const * np.ones((len(X2), len(X1), len(X2[0])))
    for i in range(len(X2)):
        for j_1 in  range(len(X1)):

            if Z1[j_1,0] == -1:
                break

            for j_2 in range(len(X2[0])):

                if Z2[i,j_2,0] == -1:
                    break

                RDist = abs(Z1[j_1,0] - Z2[i,j_2,0])
                CDist = abs(Z1[j_1,1] - Z2[i,j_2,1])
               
                #if RDist < maxRDist and CDist < maxCDist:
                dd = _MDist(X1[j_1],X2[i,j_2],width,cutDist,RWidth,Cwidth)
                dd = dd *_StochDist(RDist,CDist,RWidth,Cwidth)
               
               
                dd = D1[j_1] + D2[i,j_2] - 2 * dd
               
                if dd < 0:
                    assert dd > -1E-13, 'Something is wrong in calcProd: ' + str(dd)       
                    DistMatrix[i,j_1,j_2] = 0
               
               
                DistMatrix[i,j_1,j_2] = dd
               
    return DistMatrix

def pass_calcProdnoDist(args):
    try:
        return calcProd_noDist(*args)
    except KeyboardInterrupt:
        raise KeyboardInterruptError()

def calcProd_noDist(X1,X2,Z1,Z2,width,cutDist,RWidth,Cwidth):
   
    DistMatrix  = np.zeros((len(X2), len(X1), len(X2[0])))
    for i in range(len(X2)):
        for j_1 in  range(len(X1)):

            if Z1[j_1,0] == -1:
                break

            for j_2 in range(len(X2[0])):

                if Z2[i,j_2,0] == -1:
                    break

                RDist = abs(Z1[j_1,0] - Z2[i,j_2,0])
                CDist = abs(Z1[j_1,1] - Z2[i,j_2,1])
               
                dd = _MDist(X1[j_1],X2[i,j_2],width,cutDist,RWidth,Cwidth)
                dd = dd *_StochDist(RDist,CDist,RWidth,Cwidth)
                DistMatrix[i,j_1,j_2] = dd
               
    return DistMatrix

def pass_nonmNoDist(args):
    try:
        return nonmNoDist(*args)
    except KeyboardInterrupt:
        raise KeyboardInterruptError()

def nonmNoDist(X1,Z1,width,cutDist,RWidth,Cwidth):
   
    DistMatrix  = np.zeros((len(X1), len(X1)))
    for j_1 in  range(len(X1)):
        if Z1[j_1,0] == -1:
            break
        for j_2 in range(len(X1)):
            if Z1[j_2,0] == -1:
                break
            RDist = abs(Z1[j_1,0] - Z1[j_2,0])
            CDist = abs(Z1[j_1,1] - Z1[j_2,1])
            dd = _MDist(X1[j_1],X1[j_2],width,cutDist,RWidth,Cwidth)
            dd = dd *_StochDist(RDist,CDist,RWidth,Cwidth)
            DistMatrix[j_1,j_2] = dd
               
    return DistMatrix

def _MDist(X1,X2,width,cutDist,RWidth,Cwidth):
    maxGausDist = 8*width
    #maxRDist = 8*RWidth
    #maxCDist = 8*Cwidth   
    AAdist = 0

    for m_1 in range(len(X1[0])):
       
        if X1[0,m_1] > cutDist:
            break

        for m_2 in range(len(X2[0])):
            if X2[0,m_2] > cutDist:
                break
            if  abs(X2[0,m_2] - X1[0,m_1]) < maxGausDist:
                RDist = abs(X1[1,m_1] - X2[1,m_2])
                CDist = abs(X1[2,m_1] - X2[2,m_2])
                #if RDist < maxRDist and CDist < maxCDist:
                d = _dist(X1[0,m_1], X2[0,m_2],width,cutDist)
                d = d *_StochDist(RDist,CDist,RWidth,Cwidth)
                AAdist += d * (1 + X1[3,m_1]*X2[3,m_2] + X1[4,m_1]*X2[4,m_2])
    return AAdist   

def _dist(x1,x2,w1,w2):
    return np.exp(-((x1-x2)**2)/(4*w1**2))*(1 - np.sin(np.pi * x1/(2 * w2)))*(1 - np.sin(np.pi * x2/(2 * w2))) #*(math.cos(pi * x1/(2 * w2))*math.cos(pi * x2/(2 * w2)))
                      
                      
def _StochDist(D1,D2,W1,W2):
    rdist =  W1**2/(W1**2 + D1**2) * W2**2/(W2**2 + D2**2)
    # print "rdist", rdist
    return rdist
    # return W1**2/(W1**2 + D1**2) * W2**2/(W2**2 + D2**2)
    #math.exp(-(D1**2)/(4*W1**2) -(D2**2)/(4*W2**2))
   
def unravelKerMatrix(KM):
   
    KM_unravled = np.nansum(KM,axis = (-2,-1))
    return KM_unravled



if __name__ == "__main__":
    
    
    mols = []
    path = "tests/xyz/"
    filenames = os.listdir(path)


    print "Generating ARAD descriptors from FML interface ..."
    for filename in filenames:

        mol = fml.Molecule()
        mol.read_xyz(path + filename)
        mols.append(mol)
        mol.generate_arad_descriptor(size=30)

    
    arad = fml.ARAD()


    x1 = []
    z1 = []

    print "Generating ARAD descriptors from reference implementation ..."
    for mol in mols[:5]:

        z1.append(np.array(mol.nuclear_charges))
        x1.append(arad.describe(np.array(mol.coordinates),np.array(mol.nuclear_charges)))


    x2 = []
    z2 = []

    for mol in mols[5:]:

        z2.append(np.array(mol.nuclear_charges))
        x2.append(arad.describe(np.array(mol.coordinates),np.array(mol.nuclear_charges)))

    print "Calculating ARAD Distance matrix from reference code ..."
    Da = condense(x1,x2,z1,z2, molecular=False)

    np.set_printoptions(linewidth=10000000)

    print "Calculating ARAD Distance matrix from FML code ... "
    Db = get_l2_distance_arad(np.array(x1),np.array(x2),z1,z2)

    print "Example distance matrix (reference)"
    print Da[3,2][:10,:10]
    print Da.shape
    print "Example distance matrix (FML ARAD)"
    print Db[3,2][:10,:10]
    print Db.shape


    print "Max abs element distance difference:", np.nanmax(np.abs(Da - Db))
    print "Max abs element distance difference index:", np.nanargmax(np.abs(Da - Db))


    sigmas = [0.1, 1.0, 10.0, 100.0]

    Ka = np.zeros((len(sigmas),len(x1),len(x2)))

    print "Calculating kernel from reference distance matrix"
    for s, sigma in enumerate(sigmas):

        for i in range(len(x1)):
            for j in range(len(x2)):

                Ka[s,i,j] = np.nansum(np.exp(Da[i,j] * (-1.0 / sigma**2)))

    print "Calculating kernel from FML distance matrix"
    Kb = get_atomic_kernels_arad(np.array(x1), np.array(x2), z1, z2, sigmas)

    print "Example kernel matrix (reference)"
    print Ka[-2]


    print "Example kernel matrix (ARAD)"
    print Kb[-2]

    print "Max abs element kernel difference:", np.nanmax(np.abs(Ka - Kb))
    print "Max abs element kernel difference index:", np.nanargmax(np.abs(Ka - Kb))
