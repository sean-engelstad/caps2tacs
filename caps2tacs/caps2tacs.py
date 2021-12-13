#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:21:31 2021

@author: sengelstad6
"""

import os
import numpy as np

import pyCAPS

from tacs.pytacs import pyTACS
from tacs import functions

class Caps2Tacs:
    
    def __init__(self, csmFile, setupFunction):
        #get the CSM file and initialize the aims
        
        self.csmFile = csmFile
            
        #initialize things
        self.initAim()
        
        #give it a function to setup the loads, BCs, mesh in CAPS
        #myFunction(egadsAim,tacsAim)    
        self.setupCapsProblem(setupFunction)
        
    def initAim(self):
        #initialize the egads and tacs aims
        
        self.caps = pyCAPS.Problem('myCAPS', capsFile=self.csmFile, outLevel=0)
        
        self.geom = self.caps.geometry
        self.deskeys = self.geom.despmtr.keys()
        desDict = {}
        for key in self.deskeys:
            desDict[key] = {}
            
        self.nvar = len(self.deskeys)
        
        self.egads = self.caps.analysis.create(aim="egadsTessAIM")
        
        #self.aim is the tacs AIM
        self.aim = self.caps.analysis.create(aim = "tacsAIM",
                                         name = "tacs")
        self.aim.input.Design_Variable = desDict
        
    def setupCapsProblem(self, setupFunction):
        #run the setup function to setup egads and tacs aims with 
        #BCs, material properties, loads, etc.
        setupFunction(self.egads,self.aim)
        
        self.NumberOfNode = self.egads.output.NumberOfNode
        
    def updateDesign(self,x):
        x = np.array(x)
        #print(x)
        ind = 0
        for key in self.deskeys:
            self.geom.despmtr[key].value = x[ind]
            ind += 1
            
    def viewGeometry(self,x):
        #view the mesh
        self.updateDesign(x)
        self.aim.geometry.view()        
        
    def fullSolve(self,x=None):
        #full solve computes the function values and gradients from caps to tacs 
        #for each design vector x
    
        #update design with design vector x, if x is none use current
        if (not(x is None)):
            self.updateDesign(x)
            
        #run tacs/pytacs and reorder the sensitivity
        self.runTACS()
        self.reorderNodes()
        self.getFunctionNames()
        
        #print the sensitivity and compute func and grad from caps aim
        self.printSensitivity()
        self.getGradients()
        self.getFunctions()
            
    def runTACS(self, evalFuncs=None):
        #run tacs aim preanalysis to generate BDF and DAT file
        self.aim.preAnalysis()
        
        #read the data file into pytacs
        datFile = os.path.join(self.aim.analysisDir, self.aim.input.Proj_Name + '.dat')
        self.FEASolver = pyTACS(datFile)
        
        # Set up TACS Assembler
        self.FEASolver.initialize()
        
        #add functions
        if (evalFuncs is None):
            evalFuncs = ['wing_mass', 'ks_vmfailure']
        
        # Read in forces from BDF and create tacs struct problems
        SPs = self.FEASolver.createTACSProbsFromBDF()
        for caseID in SPs:
           SPs[caseID].addFunction('wing_mass', functions.StructuralMass)
           SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)

        # Solve each structural problem and write solutions
        self.pytacsFunc = {}; self.pytacsSens = {}
        for caseID in SPs:
            SPs[caseID].solve()
            SPs[caseID].evalFunctions(self.pytacsFunc,evalFuncs=evalFuncs)
            SPs[caseID].evalFunctionsSens(self.pytacsSens,evalFuncs=evalFuncs)
            SPs[caseID].writeSolution(outputDir=os.path.dirname(__file__))
            
        #get the keys for each function
        self.funcKeys = self.pytacsFunc.keys()
        self.nfunc = len(self.funcKeys)
        
    def reorderNodes(self):
        self.meshSens = {}
        #reorder the nodes
        for key in self.funcKeys:
            dfdX = self.pytacsSens[key]['Xpts']
            
            #dfdX is 3*nnodes want it to be nnodes x 3
            nnodes = int(len(dfdX)/3)
            
            #get the nodemap from pytacs
            bdfNodes = np.linspace(1.0,nnodes,nnodes)
            tacsNodeMap = self.FEASolver.meshLoader.getLocalNodeIDsFromGlobal(bdfNodes,nastranOrdering=True)
            
            #             
            dfdX = dfdX.reshape((nnodes,3))
            dfdX_bdf = np.zeros((nnodes,3))
            
            for bdfind in range(nnodes):
                    #node 1 in bdf, corresponds to 0 bdfind, and plug in bdfind-1 to get tacsInd cause cyclically off by 1
                    tacsInd = tacsNodeMap[bdfind] #tacsNodeMap is cyclically off by 1
                    dfdX_bdf[bdfind,:] = dfdX[tacsInd,:]
            #print(coords_nastran)
            #put back into funcSens dict
            self.meshSens[key] = dfdX_bdf
            
    def printSensitivity(self):
        #print the sensitivity file
        
        printNodes = False #setting to print nodes as first entry
        
        #where to print .sens file
        filename = os.path.join(self.aim.analysisDir, self.aim.input.Proj_Name+".sens")
        
        #open the file
        with open(filename, "w") as f:
            
            #write (nfunctions, nnodes) in first lien
            f.write("{} {}\n".format(self.nfunc,self.NumberOfNode))
            
            #for each function mass, stress, etc.
            for key in self.funcKeys:
                
                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                cSens = self.meshSens[key]
                
                #write the value of the function
                f.write(key + "\n")
                f.write("{}\n".format(self.pytacsFunc[key]))
                
                #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
                for nodeind in range(self.NumberOfNode): # d(Func1)/d(xyz)
                    if (printNodes):
                        f.write("{} {} {}\n".format(nodeind, cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
                    else:
                        f.write("{} {} {}\n".format(cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
        
        #run aim postanalysis
        self.aim.postAnalysis()
        
    def dict2vec(self,dictionary):
        keys = dictionary.keys()
        vector = np.zeros((len(keys)))
        index = 0
        for key in keys:
            vector[index] = dictionary[key]
            index += 1
        return vector
    
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
    
    def gradientDict2vec(self, dictionary):
        keys1 = list(dictionary.keys())
        keys2 = dictionary[keys1[0]].keys()
        
        newdict = {}
        for key1 in keys1:
            vec = np.zeros((len(keys2)))
            index = 0
            for key2 in keys2:
                vec[index] = dictionary[key1][key2]
                index += 1
            newdict[key1] = vec
        return newdict
    
    def getGradients(self, asDict=True):
        #get the chain rule gradient tacs->caps, from caps tacs aim
        
        desKeys = self.aim.input.Design_Variable.keys()
        sens = {}
        for key in self.funcKeys:
            sens[key] = {}
            for desKey in desKeys:
                sens[key][desKey] = self.aim.dynout[key].deriv(desKey)
        self.gradientDict = sens
        self.gradients = self.gradientDict2vec(sens)
            
        if (asDict):
            return self.gradientDict
        else:
            return self.gradients
    def getFunctions(self, asDict=True):
        #get the function values from the caps tacs aim
            
        func = {}
        for key in self.funcKeys:
            func[key] = self.aim.dynout[key].value
            
        self.functionDict = func
        self.functions = self.dict2vec(func)
    
        if (asDict):
            return self.functionDict
        else:
            return self.functions
        
    def printValues(self):
        print("\n"); print("\n")
        for key in self.funcKeys:
            print(key + " = ",self.functionDict[key])
            print("gradient = ",self.gradients[key])
            print("\n")
    
############################################################
    #Output Functions section
        
    def getFunctionNames(self):
        self.massStr = ""
        self.stressStr = ""
        for key in self.funcKeys:
            #print("Key: ",key)
            if ("mass" in key):
                self.massStr = key
            if ("failure" in key):
                self.stressStr = key
    
    def reSolve(self,x=None):
        if (not(x is None)):
            self.fullSolve(x)
    
    def mass(self, x=None):
        self.reSolve(x)
        return self.functionDict[self.massStr]
    
    def vm_stress(self,x=None):
        self.reSolve(x)
        return self.functionDict[self.stressStr]
    
    def mass_grad(self,x=None):
        self.reSolve(x)
        return self.gradients[self.massStr]
    
    def vm_stress_grad(self,x=None):
        self.reSolve(x)
        return self.gradients[self.stressStr]
    
    def checkGradients(self, x=None, h=1e-4):
        #Finite Difference Check
        functions = [self.mass,self.vm_stress]
        gradients = [self.mass_grad, self.vm_stress_grad]
        names = ["mass","ks_vm_stress"]
        
        if (x is None):
            x = np.ones((self.nvar))
        
        errors = []
        for i in range(2):
            myfunc = functions[i]
            mygrad = gradients[i]
            name = names[i]
            
            x = np.array(x)        
            p = np.random.uniform(size=x.shape)
            p = p / np.linalg.norm(p)
            
            #print(p)
            
            func1 = myfunc(x-p*h)
            func2 = myfunc(x+p*h)
            
            fdGrad = (func2-func1)/2/h
            cgrad = mygrad(x)
            #print(mygrad,p)
            directDeriv = np.matmul(cgrad, p)
            
            error = (directDeriv - fdGrad)/fdGrad
            errors.append(error)
            
        #print the functions
        print("\n"); print("\n")
        for i in range(2):
            names = name[i]; error = errors[i]
            print(name + ' FD gradient error',error)
        print("\n")
    def printDesignVariables(self,x):
        print("")
        print("Design Variables::")
        
        desDict = self.aim.input.Design_Variable
        ind = 0
        for key in desDict.keys():
            print(key, " = ", x[ind])
            ind += 1