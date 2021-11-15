#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 15:11:54 2021

@author: sengelstad6
"""

import os

import pyCAPS

from tacs.pytacs import pyTACS

from tacs import functions

import time

class CAPS2TACS:
    @classmethod
    def __init__(self,csmFile):
        self.csmFile = csmFile
        
        filename = os.path.join("stiffPanel2.csm")

        #mkdir(self.problemName + str(self.iProb))
        self.myProblem = pyCAPS.Problem('myCAPS', capsFile=filename, outLevel=0)
        
        self.geom = self.myProblem.geometry
        self.deskeys = self.geom.despmtr.keys()
        
        self.egads = self.myProblem.analysis.create(aim="egadsTessAIM")
        
        self.tacs = self.myProblem.analysis.create(aim = "tacsAIM",
                                         name = "tacs")
        self.debug = False
    @classmethod
    def cleanup(self):
        return
    @classmethod
    def generateMesh(self,x):
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
        
        self.tacs.input.Property = {"plate": shell,
                               "stiffener": shell}
        
        # Set constraints
        constraint = {"groupName" : "plateEdge",
                      "dofConstraint" : 123456}
        
        self.tacs.input.Constraint = {"edgeConstraint": constraint}
        
        # Set load
        load = {"groupName" : "plate",
                "loadType" : "Pressure",
                "pressureForce" : 2.e6}
        
        # Set loads
        self.tacs.input.Load = {"appliedPressure": load }
        
        self.tacs.input.Design_Variable = {"plateLength" : {},
                                      "plateWidth"  : {},
                                      "stiffHeight" : {}}
        
        # Run Small format
        self.tacs.preAnalysis()
    @classmethod
    def printSens(self):
        structOptions = {'writeSolution': True, }
        
        datFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name + '.dat')
        
        # Load BDF file
        FEASolver = pyTACS(datFile, options=structOptions)
        # Set up TACS Assembler
        FEASolver.createTACSAssembler()
        #add functions
        FEASolver.addFunction('wing_mass', functions.StructuralMass)
        FEASolver.addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5,
                      KSWeight=100.0)
        evalFuncs = ['wing_mass', 'ks_vmfailure']
        
        # Read in forces from BDF and create tacs struct problems
        SPs = FEASolver.createTACSProbsFromBDF()
        #print("ran pytacs")
        # Solve each structural problem and write solutions
        funcs = {}; funcsSens = {}
        for caseID in SPs:
            FEASolver(SPs[caseID])
            FEASolver.evalFunctions(SPs[caseID], funcs,evalFuncs=evalFuncs)
            FEASolver.evalFunctionsSens(SPs[caseID], funcsSens,evalFuncs=evalFuncs)
            FEASolver.writeSolution(outputDir=os.path.dirname(__file__))
        
        self.funcKeys = funcs.keys()
        nfunc = len(self.funcKeys)
        
        #dfdX = funcsSens['load_set_001_wing_mass']['Xpts']
        
        
        #reorder nodes, start at
        #TACSnodeMap = FEASolver._getGlobalToLocalNodeIDDict() #error in nodemap, less nodes
        #print(TACSnodeMap)
        
        #nnodes = int(len(dfdX)/3)
        #dfdX = dfdX.reshape(nnodes,3)
        
        #dfdX_bdf = np.zeros((nnodes,3))
        #print(len(TACSnodeMap))
        #for i in range(nnodes):
        #    #i is bdf node, tacsNode is globalNode
        #    tacsNode = TACSnodeMap[i]
        #    dfdX_bdf[i,:] = dfdX[tacsNode,:]
        
        filename = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".sens")
        with open(filename, "w") as f:
            f.write("{} {}\n".format(nfunc,self.NumberOfNode))
            for key in self.funcKeys:
                cSens = funcsSens[key]['Xpts']
                cSens = cSens.reshape(self.NumberOfNode,3)
                f.write(key + "\n")
                f.write("{}\n".format(funcs[key]))
                for nodeind in range(self.NumberOfNode): # d(Func1)/d(xyz)
                    f.write("{} {} {}\n".format(cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
    @classmethod
    def exactGrad(self,x):
        self.generateMesh(x)
        self.printSens()
        self.tacs.postAnalysis()
        
        desKeys = self.tacs.input.Design_Variable.keys()
        sens = {}
        for key in self.funcKeys:
            sens[key] = {}
            for desKey in desKeys:
                sens[key][desKey] = self.tacs.dynout[key].deriv(desKey)
        self.sens = sens
        return sens
    
    @classmethod
    def func(self,x):
        self.generateMesh(x)
        
        self.printSens()
        
        time.sleep(5)
        print('have slept')
        #self.geom.view()
        
        self.tacs.postAnalysis()
        
        func = {}
        for key in self.funcKeys:
            func[key] = self.tacs.dynout[key].value
        return func
    @classmethod
    def finiteDiff(self,x):
        self.genBDF(x) #use this to get desvar names
        sampleFunc = self.computeFunc(x) #use this to get func names
        h = 1.0e-5
        fdSens = {}
        for funcKey in sampleFunc.keys():
            fdSens[funcKey] = {}
            vind = 0
            for desvarKey in self.deskeys:
                dx = [0.0] * len(x)
                dx[vind] += h
                print(x,x+dx)
                func0 = self.func(x-dx)[funcKey]
                funcf = self.func(x+dx)[funcKey]
                #central difference finite diff
                fdSens[funcKey][desvarKey] = (funcf-func0)/h/2
                vind += 1
        return fdSens
    def compareGrad(self,x):
        return
myprob = CAPS2TACS("stiffPanel2.csm")
#print(myprob.exactGrad([1.0,1.0,1.0]))
h = 1.0e-5
print(myprob.func([1.0 ,1.0,1.0]))
myprob.debug = True
print(myprob.func([1.0 ,1.0 + h,1.0]))