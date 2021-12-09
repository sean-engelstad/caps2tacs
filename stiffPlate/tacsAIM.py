#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:21:31 2021

@author: sengelstad6
"""

import os, shutil

import pyCAPS

import matplotlib.pyplot as plt

from tacs.pytacs import pyTACS

from tacs import functions

import numpy as np

from paropt import ParOpt

from mpi4py import MPI

class TacsAim(ParOpt.Problem):
    def __init__(self,csmFile,debug=False):
        self.csmFile = csmFile
        
        csmFileDir = os.path.join(os.getcwd(),"CSM",csmFile)
        
        self.debug = debug
        
        #mkdir(self.problemName + str(self.iProb))
        self.myProblem = pyCAPS.Problem('myCAPS', capsFile=csmFileDir, outLevel=0)
        
        self.geom = self.myProblem.geometry
        self.deskeys = self.geom.despmtr.keys()
        #print(self.deskeys)
        
        self.egads = self.myProblem.analysis.create(aim="egadsTessAIM")
        
        self.tacs = self.myProblem.analysis.create(aim = "tacsAIM",
                                         name = "tacs")
        self.ct = 0
        
        # Set up the topology optimization problem
        self.nvar = 2 #num desvar
        nelems = 256
        ncon = 1
        nblock = 1
        super(TacsAim, self).__init__(MPI.COMM_SELF, self.nvar, 1, nblock)
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
        
        loadval = 1e7
        loadpernode = loadval / self.NumberOfNode 
        
        # Set load
        load = {"groupName" : "plate",
                "loadType" : "Pressure",
                "pressureForce" : loadpernode}
        
        # Set loads
        self.tacs.input.Load = {"appliedPressure": load }
        
        self.tacs.input.Design_Variable = {"plateLength" : {},
                                           "plateWidth" : {}}
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
        coordsDict = {}
        
        #reorder the nodes
        for key in self.funcKeys:
            dfdX = funcsSens[key]['Xpts']
            
            #reorder nodes, start at
            nnodes = int(len(dfdX)/3)
            coords = coords.reshape(nnodes,3)
            coords_nastran = np.zeros((nnodes,3))
            globalNodes = np.linspace(1.0,nnodes,nnodes)
            tacsNodeMap = FEASolver.meshLoader.getLocalNodeIDsFromGlobal(globalNodes,nastranOrdering=True)
            #print(tacsNodeMap)
            
            #print("dfdX shape: ",dfdX.shape)
            dfdX = dfdX.reshape((nnodes,3))
            
            dfdX_bdf = np.zeros((nnodes,3))
            #print("Length TACS node map: ",len(TACSnodeMap))
            #print("nnodes in mesh: ",nnodes)
            if (not(self.debug)):
                for bdfind in range(nnodes):
                    #node 1 in bdf, corresponds to 0 bdfind, and plug in bdfind-1 to get tacsInd cause cyclically off by 1
                    tacsInd = tacsNodeMap[bdfind] #tacsNodeMap is cyclically off by 1
                    dfdX_bdf[bdfind,:] = dfdX[tacsInd,:]
                    coords_nastran[bdfind,:] = coords[tacsInd,:]
            else: #only one comp of TACS gradient is 1, used to check CAPS sensitivity
                node = int(debugInd / 3.0)
                direc = np.mod(debugInd,3)
                dfdX_bdf[node,direc] = 1.0
            #print(coords_nastran)
            #put back into funcSens dict
            funcsSens[key]['Xpts'] = dfdX_bdf
            
        printNodes = False
        filename = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".sens")
        with open(filename, "w") as f:
            f.write("{} {}\n".format(nfunc,self.NumberOfNode))
            for key in self.funcKeys:
                cSens = funcsSens[key]['Xpts']
                f.write(key + "\n")
                f.write("{}\n".format(funcs[key]))
                for nodeind in range(self.NumberOfNode): # d(Func1)/d(xyz)
                    if (printNodes):
                        f.write("{} {} {}\n".format(nodeind, cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
                    else:
                        f.write("{} {} {}\n".format(cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
                    
        filename2 = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".coords")
        with open(filename2, "w") as f:
            f.write("Nastran Mesh Coordinates:\n")
            coords = coords_nastran
            f.write(key + "\n")
            for nodeind in range(self.NumberOfNode): # d(Func1)/d(xyz)
                f.write("{} {} {}\n".format(coords[nodeind,0], coords[nodeind,1], coords[nodeind,2]))
                    
        self.tacs.postAnalysis()
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
    def gradient(self,x):
        self.updateMesh(x)
        if (self.debug):
            nodes = [10,19,28,39,68]
            direcs = [0,1,2]
            node = nodes[4]
            direc = direcs[0]
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
    def moveBDF(self,ind):
        bdfFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".bdf")
        newLoc = os.path.join(os.getcwd(),"BDF",self.csmFile+str(ind)+".bdf")
        shutil.move(bdfFile,newLoc)
    def moveSens(self,ind):
        sensFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".sens")
        newLoc = os.path.join(os.getcwd(),"BDF",self.csmFile+str(ind)+".sens")
        shutil.move(sensFile,newLoc)
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
    def mass(self, x):
        return self.func(x)[0]
    def vm_stress(self, x):
        return self.func(x)[1]/1e-7
    def mass_grad(self, x):
        grad = self.gradient(x)
        return grad[0,:]
    def vm_stress_grad(self, x):
        grad = self.gradient(x)/1e-7
        return grad[1,:]
    def checkGradients(self, x=None, function=None, gradient=None, h=1e-4):
        if (function is None):
            function = self.func
        if (gradient is None):
            gradient = self.gradient
        if (x is None):
            x = np.ones(2)
            #x = np.ones((self.nvar))
        
        x = np.array(x)        
        p = np.random.uniform(size=x.shape)
        p = p / np.linalg.norm(p)
        
        #print(p)
        
        func1 = function(x-p*h)
        func2 = function(x+p*h)
        
        fdGrad = (func2-func1)/2/h
        mygrad = gradient(x)
        #print(mygrad,p)
        directDeriv = np.matmul(mygrad, p)
        
        #print("Functions: ",func1,func2)
        
        print('Finite Difference Gradient')
        print(fdGrad)
        print('Chain Rule Gradient')
        print(directDeriv)
        
        error = (directDeriv - fdGrad)/fdGrad
        print('Gradient Error')
        print(error)
        return abs(error) #error vector for each function
    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        lb[:] = 1e-3
        ub[:] = 3.0
        x[:] = 0.95
        return

    def evalObjCon(self, x):
        """
        Return the objective, constraint and fail flag
        """
        self.ct += 1
        print(self.ct)
        maxStress = 10.0
        fail = 0
        obj = self.mass(x[:]) #mass
        con = [maxStress - self.vm_stress(x[:])] #vmstress

        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """

        fail = 0
        g[:] = self.mass_grad(x[:])
        A[0][:] = -self.vm_stress_grad(x[:])

        return fail
    
##### MAIN #######
tacsproj = TacsAim("panel.csm",debug=False)

tacsproj.checkGradients()

options = {
    'algorithm': 'tr',
    'tr_init_size': 0.05,
    'tr_min_size': 1e-6,
    'tr_max_size': 10.0,
    'tr_eta': 0.25,
    'tr_infeas_tol': 1e-6,
    'tr_l1_tol': 1e-3,
    'tr_linfty_tol': 0.0,
    'tr_adaptive_gamma_update': True,
    'tr_max_iterations': 1000,
    'max_major_iters': 100,
    'penalty_gamma': 1e3,
    'qn_subspace_size': 10,
    'qn_type': 'bfgs',
    'abs_res_tol': 1e-8,
    'starting_point_strategy': 'affine_step',
    'barrier_strategy': 'mehrotra_predictor_corrector',
    'use_line_search': False}

options = {
    'algorithm': 'mma'}

# Set up the optimizer
opt = ParOpt.Optimizer(tacsproj, options)

#Set a new starting point
opt.optimize()
x, z, zw, zl, zu = opt.getOptimizedPoint()
print(x[:])
print(tacsproj.evalObjCon(x[:]))