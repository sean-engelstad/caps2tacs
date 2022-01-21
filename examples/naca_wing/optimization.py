# -*- coding: utf-8 -*-
import os
from caps2tacs import Caps2Tacs
from paropt import ParOpt
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions

def capsFunction(egadsAim,tacsAim):
    #setup function for panel.csm
    
    #Egads Aim section, for mesh
    egadsAim.input.Edge_Point_Min = 5
    egadsAim.input.Edge_Point_Max = 10
    
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
    
    tacsAim.input.Material = {"madeupium": madeupium}
    
    # Material properties section
    ribshell  = {"propertyType" : "Shell",
              "membraneThickness" : 0.02,
              "material"        : "madeupium",
              "bendingInertiaRatio" : 1.0, # Default
              "shearMembraneRatio"  : 5.0/6.0} # Default
    
    sparshell = {"propertyType" : "Shell",
              "membraneThickness" : 0.05,
              "material"        : "madeupium",
              "bendingInertiaRatio" : 1.0, # Default
              "shearMembraneRatio"  : 5.0/6.0} # Default
    
    tacsAim.input.Property = {"wing": ribshell}
    
    # constraint section
    constraints = {"groupName" : "wingRoot",
                  "dofConstraint" : 123456}
    
    tacsAim.input.Constraint = {"fixRoot": constraints}
    
    # Set load
    liftload = {"groupName"         : "bottomWing",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 1.0e2,
        "directionVector"   : [0.0, 1.0, 0.0]}
    
    # Set loads
    tacsAim.input.Load = {"lift": liftload }

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
        self.problem = Caps2Tacs("nacaFD_wing.csm", capsFunction, pytacsFunction, printNodes=True)
        
        self.nvar = 6 #number of design var
        ncon = 2
        nblock = 1
        super(Optimization, self).__init__(MPI.COMM_SELF, self.nvar, ncon, nblock)
        self.setBounds()
    def setBounds(self,maxStress=1.0,minStress=0.0):
        self.maxStress = maxStress
        self.minStress = minStress
    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        #area, aspectRatio, taper, twistAngle, leadingEdgeSweep, dihedral
        lowerBounds =   [20.0, 3.0,  0.3,  1.0,  3.0,  1.0 ]
        initialValues = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0 ]
        upperBounds =   [100.0,10.0, 10.0, 10.0, 50.0, 20.0]
        for i in range(self.nvar):
            lb[i] = lowerBounds[i]
            ub[i] = upperBounds[i]
            x[i] = initialValues[i]
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
        self.problem.solveStructuralProblem(x[:])
        
        massKey, stressKey = self.getNames()

        fail = 0
        obj = self.problem.func[massKey] #mass

        maxConstr = self.maxStress - self.problem.func[stressKey]
        minConstr = self.problem.func[stressKey] - self.minStress
        con = [maxConstr,minConstr] #vmstress

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
        
        stress_grad = self.problem.grad[stressKey]
        A[0][:] = -stress_grad
        A[1][:] = stress_grad
        
        return fail
    def printObjCon(self, x):
        fail, obj, con = self.evalObjCon(x)
        print("\n")
        print("Objective Value: ",obj)
        print("Constraint Value(s): ",con)

## Optimization problem defined here ##

#run options: check, run
option = "check"
#option = "run"

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
    
    x = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0 ]
    funcs = [mass,stress]
    gradients = [massGrad,stressGrad]
    names = ["mass","stress"]
    myOpt.problem.checkGradients(x,funcs,gradients,names,h=1e-4)
elif (option == "run"):
    myOpt.setBounds(maxStress=2.0,minStress=0.5)
    
    options = {
    'algorithm': 'mma',
    'mma_init_asymptote_offset': 0.5,
    'mma_min_asymptote_offset': 0.01,
    'mma_bound_relax': 1e-4,
    'mma_max_iterations': 100}
    
    # Set up the optimizer
    opt = ParOpt.Optimizer(myOpt, options)
    
    #Set a new starting point
    opt.optimize()
    x, z, zw, zl, zu = opt.getOptimizedPoint()
    
    #print optimized information
    myOpt.printObjCon(x)
    print("\n")
    print("Final design: ")
    myOpt.problem.printDesignVariables(x[:])
