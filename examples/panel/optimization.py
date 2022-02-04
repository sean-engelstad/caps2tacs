# -*- coding: utf-8 -*-
import os
from caps2tacs import Caps2Tacs
from paropt import ParOpt
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions
import numpy as np
import matplotlib.pyplot as plt

def capsFunction(egadsAim,tacsAim):
    #setup function for panel.csm
    
	#Egads Aim section, for mesh
    egadsAim.input.Edge_Point_Min = 15
    egadsAim.input.Edge_Point_Max = 20
    
    egadsAim.input.Mesh_Elements = "Quad"
    
    egadsAim.input.Tess_Params = [.25,.01,15]
    
    #increase the precision in the BDF file
    tacsAim.input.File_Format = "Large"
    tacsAim.input.Mesh_File_Format = "Large"
    
    # Link the mesh
    tacsAim.input["Mesh"].link(egadsAim.output["Surface_Mesh"])
    
    # Set analysis type
    tacsAim.input.Analysis_Type = "Static"
    
    #materials section
    madeupium    = {"materialType" : "isotropic",
                    "youngModulus" : 72.0E9 ,
                    "poissonRatio": 0.33,
                    "density" : 2.8E3,
                    "tensionAllow" :  20.0e7}
    
    tacsAim.input.Material = {"Madeupium": madeupium}
    
    # Material properties section
    shell  = {"propertyType" : "Shell",
              "membraneThickness" : 0.006,
              "material"        : "madeupium",
              "bendingInertiaRatio" : 1.0, # Default
              "shearMembraneRatio"  : 5.0/6.0} # Default
    
    tacsAim.input.Property = {"plate": shell}
    
    # constraint section
    constraint = {"groupName" : "edge",
                  "dofConstraint" : 123456}
    
    tacsAim.input.Constraint = {"edgeConstraint": constraint}
    
    #loads section
    pressload = 1.0e5
    
    # Set load
    load = {"groupName" : "plate",
            "loadType" : "Pressure",
            "pressureForce" : pressload}
    
    # Set loads
    tacsAim.input.Load = {"appliedPressure": load }
    
    #design variables
    tacsAim.input.Design_Variable = {"plateLength" : {},
    								"plateWidth" : {}}
    
def pytacsFunction(obj, datFile):
    #locally defined method that runs pytacs with your desired functions
    #and prints the function values and sensitivity, for each data file & bdf pair

    #initialize pytacs with that data file
    obj.FEASolver = pyTACS(datFile)
        
    # Set up TACS Assembler
    obj.FEASolver.initialize()
    
    #choose the functions to evaluate
    evalFuncs = ['wing_mass', 'ks_vmfailure']
    
    #read the bdf & dat file into pytacs FEAsolver
    #SPs represents "StructuralProblems"
    SPs = obj.FEASolver.createTACSProbsFromBDF()
    
    # Read in forces from BDF and create tacs struct problems
    for caseID in SPs:
       SPs[caseID].addFunction('wing_mass', functions.StructuralMass)
       SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)

    # Solve each structural problem and write solutions
    func = {}; sens = {}
    for caseID in SPs:
        SPs[caseID].solve()
        SPs[caseID].evalFunctions(func,evalFuncs=evalFuncs)
        SPs[caseID].evalFunctionsSens(sens,evalFuncs=evalFuncs)
        SPs[caseID].writeSolution(outputDir=os.path.dirname(__file__))

    #store function and sensitivity values    
    obj.func = func
    obj.sens = sens

#ParOpt Optimization Class
class Optimization(ParOpt.Problem):
    def __init__(self):
        desvarList = ["plateLength", "plateWidth"]
        self.problem = Caps2Tacs("panel.csm", capsFunction, pytacsFunction, desvarList)
        
        nvar = 2 #number of design var
        #ncon = 2
        ncon = 0
        nblock = 1
        super(Optimization, self).__init__(MPI.COMM_SELF, nvar, ncon, nblock)

        self.objs = []
        self.cons = []
    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        lb[:] = 1e-3
        ub[:] = 3.0
        x[:] = 0.95
        return
    def getNames(self):
        massStr = ""
        stressStr = ""
        for key in self.problem.funcKeys:
            #print("Key: ",key)
            if ("mass" in key):
                massStr = key
            if ("failure" in key):
                stressStr = key
        return massStr, stressStr
    def evalObjCon(self, x):
        """
        Return the objective, constraint and fail flag
        """
        #run the solver
        self.problem.solveStructuralProblem(x[:])
        
        massKey, stressKey = self.getNames()

        fail = 0
        obj = self.problem.func[stressKey] #mass
        self.objs.append(obj)

        #maxConstr = 2.0 - self.problem.func[stressKey]
        #minConstr = self.problem.func[stressKey] - 0.5
        #con = [maxConstr,minConstr] #vmstress
        con = []

        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """
        #run the solver
        self.problem.solveStructuralProblem(x[:])
        
        massKey, stressKey = self.getNames()
        
        fail = 0
        g[:] = self.problem.grad[stressKey]
        
        #stress_grad = self.problem.grad[stressKey]
        #A[0][:] = -stress_grad
        #A[1][:] = stress_grad
        
        return fail
    def printObjCon(self, x):
        fail, obj, con = self.evalObjCon(x)
        print("\n")
        print("Objective Value: ",obj)
        print("Constraint Value(s): ",con)
    def plotObjectiveProgress(self):
        niter = len(self.objs)
        iterations = np.arange(0,niter)
        plt.plot(iterations, self.objs, 'k-')
        plt.xlabel('Iteration #')
        plt.ylabel('stress obj')
        plt.show()

## Optimization problem defined here ##

#run options: check, run
option = "check"
    
myOpt = Optimization()


if (option == "check"):
    def mass(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey = myOpt.getNames()
        return myOpt.problem.func[massKey]
    def stress(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey = myOpt.getNames()
        return myOpt.problem.func[stressKey]
    def massGrad(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey = myOpt.getNames()
        return myOpt.problem.grad[massKey]
    def stressGrad(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey = myOpt.getNames()
        return myOpt.problem.grad[stressKey]
    
    x = [1.0, 1.0 ]
    funcs = [mass,stress]
    gradients = [massGrad,stressGrad]
    names = ["mass","stress"]
    myOpt.problem.checkGradients(x,funcs,gradients,names,h=1e-4)
    
elif (option == "run"):
   
    options = {
        'algorithm': 'mma'}
    
    # Set up the optimizer
    opt = ParOpt.Optimizer(myOpt, options)
    
    #Set a new starting point
    opt.optimize()
    x, z, zw, zl, zu = opt.getOptimizedPoint()
    
    #print optimized information
    myOpt.printObjCon(x)
    print("\n")
    print("Final Design: ")
    myOpt.problem.printDesignVariables(x[:])

    myOpt.plotObjectiveProgress()