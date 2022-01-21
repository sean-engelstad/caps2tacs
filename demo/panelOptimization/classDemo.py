#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:21:31 2021

@author: sengelstad6
"""

import os, shutil

import pyCAPS

from tacs.pytacs import pyTACS

from tacs import functions

import numpy as np

from paropt import ParOpt

from mpi4py import MPI

class Caps2Tacs(ParOpt.Problem):
    def __init__(self,csmFile,debug=False):
        self.csmFile = csmFile
        
        csmFileDir = csmFile
        
        self.debug = debug
        
        self.myProblem = pyCAPS.Problem('myCAPS', capsFile=csmFileDir, outLevel=0)
        
        self.geom = self.myProblem.geometry
        self.deskeys = self.geom.despmtr.keys()
        #print(self.deskeys)
        
        self.egads = self.myProblem.analysis.create(aim="egadsTessAIM")
        
        self.tacs = self.myProblem.analysis.create(aim = "tacsAIM",
                                         name = "tacs")
        self.ct = 0
        
        self.setBounds()
        
        self.changeLoad = False
        
        # Set up the topology optimization problem
        self.nvar = 2 #num desvar
        nelems = 256
        ncon = 2
        nblock = 1
        super(Caps2Tacs, self).__init__(MPI.COMM_SELF, self.nvar, ncon, nblock)
    def loadMag(self,x):
        return 1/(x[0]*x[1])**(1.5)
    def loadMagGrad(self,x):
        g = np.zeros((2))
        coeff = -1.5/(x[0]*x[1])**(2.5)
        g[0] = coeff * x[1]
        g[1] = coeff * x[0]
        return g
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
        
        if (self.changeLoad):
        	pressload = 1.0e5 * self.loadMag(x)
        else:
        	pressload = 1.0e5
        	
        # Set load
        load = {"groupName" : "plate",
                "loadType" : "Pressure",
                "pressureForce" : pressload}
        
        # Set loads
        self.tacs.input.Load = {"appliedPressure": load }
        
        self.tacs.input.Design_Variable = {"plateLength" : {},
                                           "plateWidth" : {}}
    def printSens(self,debugInd=-1):
        self.tacs.preAnalysis()
        #doesn't take struct Options right now, wan't to remove print
        structOptions = {'writeConnectivity': False, 
                         'writeNodes': False,
                         'writeDisplacements': False,
                         'writeStrains': False,
                         'writeStresses': False,
                         'writeExtras': False,
                         'printIterations': False}
        
        datFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name + '.dat')
        
        # Load BDF file
        FEASolver = pyTACS(datFile)
        # Set up TACS Assembler
        FEASolver.initialize()
        #add functions
        evalFuncs = ['wing_mass', 'ks_vmfailure']
        
        # Read in forces from BDF and create tacs struct problems
        SPs = FEASolver.createTACSProbsFromBDF()
        for caseID in SPs:
               SPs[caseID].addFunction('wing_mass', functions.StructuralMass)
               SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)
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
            
        printNodes = True
        filename = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".sens")
        with open(filename, "w") as f:
            f.write("{} {}\n".format(nfunc,self.NumberOfNode))
            for key in self.funcKeys:
                cSens = funcsSens[key]['Xpts']
                f.write(key + "\n")
                f.write("{}\n".format(funcs[key]))
                for nodeind in range(self.NumberOfNode): # d(Func1)/d(xyz)
                    if (printNodes):
                        f.write("{} {} {} {}\n".format(nodeind+1, cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
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
        #store gradient to be pulled separately
        self.grad = sensMatrix        
        return self.grad
    def moveBDF(self,ind):
        bdfFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".bdf")
        newLoc = os.path.join(os.getcwd(),"BDF",self.csmFile+str(ind)+".bdf")
        shutil.move(bdfFile,newLoc)
    def moveSens(self,ind):
        sensFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".sens")
        newLoc = os.path.join(os.getcwd(),"BDF",self.csmFile+str(ind)+".sens")
        shutil.move(sensFile,newLoc)
    def function(self,x,getVec=True):
        self.updateMesh(x)
        if (self.debug):
            self.printSens(0)
        else:
            self.printSens()
            
        func = {}
        for key in self.funcKeys:
            func[key] = self.tacs.dynout[key].value
        if (getVec):
            self.func = self.dict2vec(func)
        else:
            self.func = func
        return self.func
    def setBounds(self,maxStress=1):
        self.maxStress = maxStress
    def mass(self,x=None):
        if (x is None):
            m = self.func[1]
        else:
            m = self.function(x)[1]
        return m
    def vm_stress(self,x=None):
        if (x is None):
            stress = self.func[0]
        else:
            stress = self.function(x)[0]
        return stress
    def mass_grad(self,x=None):
        if (x is None):
            mgrad = self.grad[1,:]
        else:
            mgrad = self.gradient(x)[1,:]
        return mgrad
    def vm_stress_grad(self, x=None):
        if (x is None):
            stress_grad = self.grad[0,:]
            stress = self.func[0]
        else:
            stress_grad = self.gradient(x)[0,:]
            stress = self.function(x)[0]
        
        if (self.changeLoad):
            #chain rule for changing load
            sgrad = stress_grad * self.loadMag(x) + stress * self.loadMagGrad(x)
        else:
            sgrad = stress_grad
        return sgrad
    
    def checkGradients(self, x=None, functions=None, gradients=None, names = None,h=1e-4):
        if (functions is None):
            functions = [self.mass,self.vm_stress]
        if (gradients is None):
            gradients = [self.mass_grad, self.vm_stress_grad]
        if (names is None):
            names = ["mass","ks_vm_stress"]
        if (x is None):
            x = np.ones(2)
            #x = np.ones((self.nvar))
        
        self.function(x); self.gradient(x)
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
        for i in range(2):
            name = names[i]; error = errors[i]
            print(name + ' FD gradient error',error)
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
        #minimize mass, with maxstress limit, 
        #call function
        self.function(x)
        
        self.ct += 1
        print(self.ct)
        fail = 0
        obj = self.mass(x[:]) #mass
        stress = self.vm_stress(x[:])
        con = [2.0 - stress, stress-0.5] #vmstress

        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """
        #compute gradients
        self.gradient(x)
        
        fail = 0
        g[:] = self.mass_grad(x[:])
        stress_grad = self.vm_stress_grad(x[:])
        A[0][:] = -stress_grad
        A[1][:] = stress_grad

        return fail
