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
              "membraneThickness" : 0.01,
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
    
    #geom design variables
    DVdict = {}
    geomDVs = ["plateLength","plateWidth"]
    for geomDV in geomDVs:
        DVdict[geomDV] = {}

    #thick design variables
    def makeThicknessDV(capsGroup, thickness):
      desvar    = {"groupName" : capsGroup,
              "initialValue" : thickness,
              "lowerBound" : thickness*0.5,
              "upperBound" : thickness*1.5,
              "maxDelta"   : thickness*0.1,
              "fieldName" : "T"}
      return desvar

    thickDV = makeThicknessDV("plate",0.01)
    DVdict["thick"] = thickDV
        
    tacsAim.input.Design_Variable = DVdict
    
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
        desvarList = ["plateLength", "plateWidth","thick"]
        self.problem = Caps2Tacs("panel.csm", capsFunction, pytacsFunction, desvarList)
        
        nvar = 3 #number of design var
        ncon = 1
        nblock = 1
        super(Optimization, self).__init__(MPI.COMM_SELF, nvar, ncon, nblock)

        self.objs = []
        self.cons = []
        self.start = True

    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        lbs = [0.1, 0.1, 0.01]
        ubs = [1.0, 1.0, 0.1]
        xs = [0.5, 0.5, 0.05]
        for i in range(3):
            lb[i] = lbs[i]
            ub[i] = ubs[i]
            x[i] = xs[i]
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
        con = np.zeros(1, dtype=ParOpt.dtype)

        #run the solver
        self.problem.solveStructuralProblem(x[:])
        
        massKey, stressKey = self.getNames()

        #setup initial max stress
        if (self.start):
            self.maxStress = 2 * self.problem.func[stressKey]
            self.start = False

        fail = 0
        obj = self.problem.func[massKey] #mass
        self.objs.append(obj)

        con[0] = 1 - self.problem.func[stressKey] / self.maxStress

        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """
        #run the solver
        self.problem.solveStructuralProblem(x[:])
        
        massKey, stressKey = self.getNames()
        
        fail = 0
        g[:] = self.problem.grad[massKey]
        print(g[:])
        A[0][:] = -1 * self.problem.grad[stressKey] / self.maxStress
        print(A[0][:])

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

#run options: check, run, eval
option = "run"
    
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
    
    x = [1.0, 1.0, 0.01 ]
    funcs = [mass,stress]
    gradients = [massGrad,stressGrad]
    names = ["mass","stress"]
    myOpt.problem.checkGradients(x,funcs,gradients,names,h=1e-4)
    
elif (option == "run"):
   
    # options = {
    # 'algorithm': 'mma',
    # 'mma_init_asymptote_offset': 0.5,
    # 'mma_min_asymptote_offset': 0.01,
    # 'mma_max_iterations': 20}
    filename = "paropt.out"
    options = {
        'algorithm': 'ip',
        'abs_res_tol': 1e-4,
        'starting_point_strategy': 'affine_step',
        'barrier_strategy': 'monotone',
        'start_affine_multiplier_min': 0.01,
        'penalty_gamma': 1000.0,
        'qn_subspace_size': 10,
        'qn_type': 'bfgs',
        'output_file': filename}
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

elif (option == "eval"):

    D1 = [1.0,1.0,0.05]
    D2 = [1.0, 1.0, 0.04]
    D3 = [1.0,1.0, 0.03]

    myOpt.problem.solveStructuralProblem(D1[:])
    myOpt.problem.printResults()
    myOpt.problem.solveStructuralProblem(D2[:])
    myOpt.problem.printResults()
    myOpt.problem.solveStructuralProblem(D3[:])
    myOpt.problem.printResults()