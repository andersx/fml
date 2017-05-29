# MIT License
#
# Copyright (c) 2017 Felix Faber and Anders Steen Christensen
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

from numpy import *
import numpy as np
from copy import deepcopy
from pydispersion import getc6, getc6_pbc

class ARAS(object):

    def getAngle(self,sp,norms):
        angles = np.zeros(sp.shape)
        mask1 = np.logical_and(np.abs(sp - norms) > self.epsilon ,np.abs(norms) > self.epsilon)
        angles[mask1] = np.arccos(np.clip(sp[mask1]/norms[mask1], -1.0, 1.0))
        return angles

    def __init__(self,maxMolSize = 30,maxAts = 30,cut = 5., debug=False):
        self.tag = 'coords'
        self.debug = debug
        self.maxMolSize = maxMolSize
        self.maxAts = maxAts
        self.cut = cut
        self.epsilon = 10.* np.finfo(float).eps

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
        coords = asarray(coords)
        ocupationList = asarray(ocupationList)
        M =  zeros((self.maxMolSize,3+self.maxAts*2,self.maxAts))

        if cell is not None:
            coords = dot(coords,cell)
            nExtend = (floor(self.cut/linalg.norm(cell,2,axis = 0)) + 1).astype(int)
            #print nExtend
            for i in range(-nExtend[0],nExtend[0] + 1):
                for j in range(-nExtend[1],nExtend[1] + 1):
                    for k in range(-nExtend[2],nExtend[2] + 1):

                        #print i, j, k
                        #print i*cell[0,:] + j*cell[1,:] + k*cell[2,:]
                        if i == -nExtend[0] and j  == -nExtend[1] and k  == -nExtend[2]:
                            coordsExt = coords + i*cell[0,:] + j*cell[1,:] + k*cell[2,:]
                            ocupationListExt = copy(ocupationList)
                        else:
                            ocupationListExt = append(ocupationListExt,ocupationList)
                            coordsExt = append(coordsExt,coords + i*cell[0,:] + j*cell[1,:] + k*cell[2,:],axis = 0)

        else:
            coordsExt = copy(coords)
            ocupationListExt = copy(ocupationList)

        M[:,0,:] = 1E+100

        C6 = None
        
        if cell is not None:
            C6 = getc6_pbc(coords, len(ocupationList), ocupationList, cell)
        else:
            C6 = getc6(coords, len(ocupationList), ocupationList)

        C9 = np.sqrt(C6[np.newaxis]*C6[:,np.newaxis]*C6[:,:,np.newaxis])

        # print "C6", C6.shape
        # print C6
        # print "C9", C9.shape
        # print C9
        # print "done"

        for i in range(L):
            c6 = deepcopy(C6[i])
            c9 = deepcopy(C9[i])
            #print '+'
            #Calculate Distance
            cD = - coords[i] + coordsExt[:]

            ocExt =  asarray(ocupationListExt)

            #Obtaining angles
            sp = sum(cD[:,newaxis] * cD[newaxis,:], axis = 2)
            D1 = sqrt(sum(cD**2, axis = 1))
            D2 = D1[:,newaxis]*D1[newaxis,:]
            angs = self.getAngle(sp, D2)
            if any(np.isnan(angs)):
                print "WARNING: Possible error in angle calculation!"
                print "scalar product: ", sp[np.isnan(angs)]
                print "norm: ", D2[np.isnan(angs)]
                print sp[np.isnan(angs)]/D2[np.isnan(angs)]
                print "ang: ", angs[np.isnan(angs)]

            args = argsort(D1)
            # print "args", args
            # print "1", c6
            c6 = c6[args]
            # print "2", c6
            c9 = c9[args,:]
            c9 = c9[:,args]
            D1 = D1[args]
            ocExt = asarray([ocExt[l] for l in args])
            angs = angs[args,:]
            angs = angs[:,args]

            args = where(D1 < self.cut)[0]
            D1 = D1[args]
            ocExt = asarray([ocExt[l] for l in args])
            angs = angs[args,:]
            angs = angs[:,args]
            c6 = c6[args]
            c9 = c9[args,:]
            c9 = c9[:,args]
            M[i,0,: len(D1)] = D1
            M[i,1,: len(D1)] = ocExt[:]
            M[i,3: len(D1) + 3,: len(D1)] = angs
            M[i,2,: len(D1)] = c6[:]
            M[i,3+ self.maxAts: self.maxAts + len(D1) + 3,: len(D1)] = c9[:]
        return M
