#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:21:31 2021

@author: sengelstad6
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 15:11:54 2021

@author: sengelstad6
"""

import os, shutil

import pyCAPS

import matplotlib.pyplot as plt

from tacs.pytacs import pyTACS

from tacs import functions

import numpy as np

class tacsaim:
    @classmethod
    def __init__(self,csmFile,debug=False):
        self.csmFile = csmFile
        
        csmFileDir = os.path.join(os.getcwd(),"CSM",csmFile)
        
        self.debug = debug
        
        #mkdir(self.problemName + str(self.iProb))
        self.myProblem = pyCAPS.Problem('myCAPS', capsFile=csmFileDir, outLevel=0)
        
        self.geom = self.myProblem.geometry
        self.deskeys = self.geom.despmtr.keys()
        print(self.deskeys)
        
        self.egads = self.myProblem.analysis.create(aim="egadsTessAIM")
        
        self.tacs = self.myProblem.analysis.create(aim = "tacsAIM",
                                         name = "tacs")
    @classmethod
    def updateMesh(self,x):
        x = np.array(x)
        #print(x)
        ind = 0
        for key in self.deskeys:
            self.geom.despmtr[key].value = x[ind]
            ind += 1
        
        #geom.view()
        #myProblem.geometry.view()
        
        self.egads.input.Edge_Point_Min = 15
        self.egads.input.Edge_Point_Max = 20
        
        self.egads.input.Mesh_Elements = "Quad"
        
        self.egads.input.Tess_Params = [.25,.01,15]
        
        self.NumberOfNode = self.egads.output.NumberOfNode
        
        self.tacs.input.File_Format = "Large"
        self.tacs.input.Mesh_File_Format = "Large"
        
        # Link the mesh
        self.tacs.input["Mesh"].link(self.egads.output["Surface_Mesh"])
        
        # Set analysis type
        self.tacs.input.Analysis_Type = "Static"
        
        madeupium    = {"materialType" : "isotropic",
                        "youngModulus" : 72.0E9 ,
                        "poissonRatio": 0.33,
                        "density" : 2.8E3,
                        "tensionAllow" :  20.0e7}
        
        self.tacs.input.Material = {"Madeupium": madeupium}
        
        # Set properties
        shell  = {"propertyType" : "Shell",
                  "membraneThickness" : 0.006,
                  "material"        : "madeupium",
                  "bendingInertiaRatio" : 1.0, # Default
                  "shearMembraneRatio"  : 5.0/6.0} # Default
        #need to add the strength here
        
        self.tacs.input.Property = {"plate": shell}
        # Set constraints
        constraint = {"groupName" : "edge",
                      "dofConstraint" : 123456}
        
        self.tacs.input.Constraint = {"edgeConstraint": constraint}
        
        # Set load
        load = {"groupName" : "plate",
                "loadType" : "Pressure",
                "pressureForce" : 2.e6}
        
        # Set loads
        self.tacs.input.Load = {"appliedPressure": load }
        
        self.tacs.input.Design_Variable = {"plateLength" : {},
                                           "plateWidth" : {}}
    @classmethod
    def printSens(self,debugInd=-1):
        self.tacs.preAnalysis()
        structOptions = {'writeSolution': True, }
        
        datFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name + '.dat')
        
        # Load BDF file
        FEASolver = pyTACS(datFile, options=structOptions)
        # Set up TACS Assembler
        FEASolver.initialize()
        #add functions
        evalFuncs = ['wing_mass', 'ks_vmfailure']
        
        # Read in forces from BDF and create tacs struct problems
        SPs = FEASolver.createTACSProbsFromBDF()
        for caseID in SPs:
               SPs[caseID].addFunction('wing_mass', functions.StructuralMass)
               SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, KSWeight=100.0)
        #print("ran pytacs")
        # Solve each structural problem and write solutions
        funcs = {}; funcsSens = {}
        for caseID in SPs:
            SPs[caseID].solve()
            SPs[caseID].evalFunctions(funcs,evalFuncs=evalFuncs)
            SPs[caseID].evalFunctionsSens(funcsSens,evalFuncs=evalFuncs)
            SPs[caseID].writeSolution(outputDir=os.path.dirname(__file__))
            coords = SPs[caseID].getNodes()
            
        self.funcKeys = funcs.keys()
        nfunc = len(self.funcKeys)
        
        #reorder the nodes
        for key in self.funcKeys:
            dfdX = funcsSens[key]['Xpts']
            
            #reorder nodes, start at
            nnodes = int(len(dfdX)/3)
            coords = coords.reshape(nnodes,3)
            coords_nastran = np.zeros((nnodes,3))
            globalNodes = np.linspace(1.0,nnodes,nnodes)
            TACSnodeMap = FEASolver.meshLoader.getLocalNodeIDsFromGlobal(globalNodes) #error in nodemap, less nodes
            #print("TACS Node Map:", TACSnodeMap)
            
            #print("dfdX shape: ",dfdX.shape)
            dfdX = dfdX.reshape(nnodes,3)
            
            dfdX_bdf = np.zeros((nnodes,3))
            #print("Length TACS node map: ",len(TACSnodeMap))
            #print("nnodes in mesh: ",nnodes)
            if (not(self.debug)):
                for bdfind in range(nnodes):
                    tacsNode = TACSnodeMap[bdfind]
                    dfdX_bdf[bdfind,:] = dfdX[tacsNode,:]
                    coords_nastran[bdfind,:] = coords[tacsNode,:]
            else: #only one comp of TACS gradient is 1, used to check CAPS sensitivity
                node = int(debugInd / 3.0)
                direc = np.mod(debugInd,3)
                dfdX_bdf[node,direc] = 1.0
            #print(coords_nastran)
            #put back into funcSens dict
            funcsSens[key]['Xpts'] = dfdX_bdf
            
        filename = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".sens")
        with open(filename, "w") as f:
            f.write("{} {}\n".format(nfunc,self.NumberOfNode))
            for key in self.funcKeys:
                cSens = funcsSens[key]['Xpts']
                f.write(key + "\n")
                f.write("{}\n".format(funcs[key]))
                for nodeind in range(self.NumberOfNode): # d(Func1)/d(xyz)
                    f.write("{} {} {}\n".format(cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
                    
        self.tacs.postAnalysis()
    @classmethod
    def dict2vec(self,dictionary):
        keys = dictionary.keys()
        vector = np.zeros((len(keys)))
        index = 0
        for key in keys:
            vector[index] = dictionary[key]
            index += 1
        return vector
    @classmethod
    def dict2matrix(self,dictionary):
        keys1 = list(dictionary.keys())
        keys2 = dictionary[keys1[0]].keys()
        
        matrix = np.zeros((len(keys1),len(keys2)))
        row = 0
        for key1 in keys1:
            col = 0
            for key2 in keys2:
                matrix[row,col] = dictionary[key1][key2]
                col += 1
            row += 1
        return matrix
    @classmethod
    def gradient(self,x):
        self.updateMesh(x)
        if (self.debug):
            nodes = [10,19,28,39,68]
            direcs = [0,1,2]
            node = nodes[4]
            direc = direcs[1]
            self.printSens(3*(node-1)+direc)
        else:
            self.printSens()
        
        desKeys = self.tacs.input.Design_Variable.keys()
        sens = {}
        for key in self.funcKeys:
            sens[key] = {}
            for desKey in desKeys:
                sens[key][desKey] = self.tacs.dynout[key].deriv(desKey)
        sensMatrix = self.dict2matrix(sens)        
        return sensMatrix
    @classmethod
    def moveBDF(self,ind):
        bdfFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".bdf")
        newLoc = os.path.join(os.getcwd(),"BDF",self.csmFile+str(ind)+".bdf")
        shutil.move(bdfFile,newLoc)
    @classmethod
    def moveSens(self,ind):
        sensFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".sens")
        newLoc = os.path.join(os.getcwd(),"BDF",self.csmFile+str(ind)+".sens")
        shutil.move(sensFile,newLoc)
    @classmethod
    def func(self,x):
        self.updateMesh(x)
        if (self.debug):
            self.printSens(0)
        else:
            self.printSens()
            
        func = {}
        for key in self.funcKeys:
            func[key] = self.tacs.dynout[key].value
        funcVec = self.dict2vec(func)
        return funcVec
    @classmethod
    def func2(self,filename):
        self.printSens2(filename)
        
        func = {}
        for key in self.funcKeys:
            func[key] = self.tacs.dynout[key].value
        funcVec = self.dict2vec(func)
        return funcVec
    @classmethod
    def finiteDifference(self, x, function, gradient, h=1e-4):
        x = np.array(x)        
        p = np.random.uniform(size=x.shape)
        p = p / np.linalg.norm(p)
        
        print(p)
        
        func1 = function(x-p*h)
        func2 = function(x+p*h)
        
        fdGrad = (func2-func1)/2/h
        mygrad = gradient(x)
        print(mygrad,p)
        directDeriv = np.matmul(mygrad, p)
        
        print("Functions: ",func1,func2)
        
        print('Finite Difference Gradient')
        print(fdGrad)
        print('Chain Rule Gradient')
        print(directDeriv)
        
        error = (directDeriv - fdGrad)/fdGrad
        print('Gradient Error')
        print(error)
        return np.linalg.norm(error)/2**0.5
    @classmethod
    def finiteDiffPlot(self,x):
        orderOfMag = -1*np.linspace(1.0,8.0,30)
        stepSizes = 10**(orderOfMag)
        errors = np.zeros((30))
        
        for i in range(len(stepSizes)):
            h = stepSizes[i]
            errors[i] = self.finiteDifference(x,tacsaim.func, tacsaim.gradient, h)
        
        fig, ax = plt.subplots()
        ax.loglog(stepSizes,errors,'k-',linewidth=3)
        ax.set_xlabel('Step Sizes')
        ax.set_ylabel('Error in Obj Gradient')
        plt.show()
    @classmethod
    def bdfChange(self,run):
        h = 1.0e-5
        L = 0.5; W = 1.0
        if (run == 1):
            x = np.array( [L,W] )
            myprob.func(x)
        elif (run == 2): 
            x = np.array( [L+h,W] )
            myprob.func(x)
        myprob.moveBDF(run)
        myprob.moveSens(run)
    @classmethod
    def checkTACSgradient(self):
        print("case 1: ",self.func2("panel.csm1"))
        print("case 2: ",self.func2("panel.csm2"))
        return
##### MAIN #######
myprob = tacsaim("panel.csm",debug=True)
#to check the sensitivities I was manually 
#changing node and direc in gradient function
print(myprob.gradient([0.5,1.0]))