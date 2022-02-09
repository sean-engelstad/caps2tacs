# -*- coding: utf-8 -*-
import os
from caps2tacs import *
from paropt import ParOpt
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions
import numpy as np
import matplotlib.pyplot as plt

def capsFunction(egadsAim,tacsAim):
    #setup function for panel.csm
    
    #Egads Aim section, for mesh
    egadsAim.input.Edge_Point_Min = 20
    egadsAim.input.Edge_Point_Max = 25
    
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
    OMLshell = {"propertyType" : "Shell",
              "membraneThickness" : 0.01,
              "material"        : "madeupium",
              "bendingInertiaRatio" : 1.0, # Default
              "shearMembraneRatio"  : 5.0/6.0} # Default
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
    
    tacsAim.input.Property = {"rib": ribshell,
    "spar" : sparshell,
    "OML" : OMLshell}
    
    # constraint section
    constraint1 = {"groupName" : "wingRoot",
                  "dofConstraint" : 123456}

    tacsAim.input.Constraint = {"fixRoot": constraint1}
    
    # Set load
    liftload = {"groupName"         : "bottomWing",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 1.0e2,
        "directionVector"   : [0.0, 1.0, 0.0]}
    
    # Set loads
    tacsAim.input.Load = {"lift": liftload }
                                     
    #geom design variables
    DVdict = {}
    DVRdict = {}
    geomDVs = ["area","aspect", "taper", "twist", "lesweep", "dihedral"]
    capsGroups = ["rib", "spar", "OML"]
    thickDVs = ["thick1", "thick2", "thick3"]
    thick0 = 0.02
    for geomDV in geomDVs:
        DVdict[geomDV] = {}

    #thick design variables
    # def makeThicknessDV(capsGroup, thickness):
    #   desvar    = {"groupName" : capsGroup,
    #           "initialValue" : thickness,
    #           "lowerBound" : thickness*0.5,
    #           "upperBound" : thickness*1.5,
    #           "maxDelta"   : thickness*0.1,
    #           "fieldName" : "T"}
    #           #"variableType" : "Property"}
    #   return desvar
    # def makeThicknessDVR(capsGroup):
    #     DVR = {"variableType": "Property",
    #     "fieldName" : "T",
    #     "constantCoeff" : 0.0,
    #     "groupName" : capsGroup,
    #     "linearCoeff" : 1.0}
    #     return DVR

    useDVR = True
    for i in range(3):
        if (useDVR):
            DVRdict[thickDVs[i]] = makeThicknessDVR(capsGroups[i])
        else:
            DVdict[thickDVs[i]] = makeThicknessDV(capsGroups[i],thick0)
        
    
    tacsAim.input.Design_Variable = DVdict   
    if (useDVR): tacsAim.input.Design_Variable_Relation = DVRdict
        

def pytacsFunction(obj, datFile):
    useMPI = True
    if (useMPI):
        comm = MPI.COMM_WORLD


    #locally defined method that runs pytacs with your desired functions
    #and prints the function values and sensitivity, for each data file & bdf pair

    #initialize pytacs with that data file
    obj.FEASolver = pyTACS(datFile)
        
    # Set up TACS Assembler
    obj.FEASolver.initialize()
    
    #choose the functions to evaluate
    evalFuncs = ['wing_mass', 'ks_vmfailure','compliance']

    # vec = obj.FEASolver.assembler.createVec()
    # vec.getArray()[:] = 1.0
    # obj.FEASolver.assembler.applyBCs(vec)
    # obj.FEASolver.assembler.setVariables(vec)
    # obj.FEASolver.outputViewer.writeToFile('bcs.f5')
    
    #read the bdf & dat file into pytacs FEAsolver
    #SPs represents "StructuralProblems"
    SPs = obj.FEASolver.createTACSProbsFromBDF()
    
    # Read in forces from BDF and create tacs struct problems
    for caseID in SPs:
       SPs[caseID].addFunction('wing_mass', functions.StructuralMass)
       SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)
       SPs[caseID].addFunction('compliance', functions.Compliance)
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
        desvarList = ["area","aspect","taper","twist","lesweep","dihedral","thick1", "thick2", "thick3"]
        self.problem = Caps2Tacs("naca_thick.csm", capsFunction, pytacsFunction, desvarList)
        
        self.nvar = 9 #number of design var
        ncon = 1 #number of constraint
        nblock = 1
        super(Optimization, self).__init__(MPI.COMM_SELF, self.nvar, ncon, nblock)

        self.objs = []
        self.start = True
    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        #area, aspectRatio, taper, twistAngle, leadingEdgeSweep, dihedral, ribT, sparT, OMLT
        lowerBounds =   [20.0, 3.0,  0.3,  1.0,  3.0,  1.0, 0.01, 0.01, 0.01 ]
        initialValues = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0, 0.03, 0.05, 0.02]
        upperBounds =   [100.0,10.0, 1.0, 10.0, 50.0, 20.0, 0.1, 0.1, 0.1]
        for i in range(self.nvar):
            lb[i] = lowerBounds[i]
            ub[i] = upperBounds[i]
            x[i] = initialValues[i]
        return
    def getNames(self):
        massStr = ""
        stressStr = ""
        compStr = ""
        for key in self.problem.funcKeys:
            #print("Key: ",key)
            if ("mass" in key):
                massStr = key
            if ("failure" in key):
                stressStr = key
            if ("compliance" in key):
                compStr = key
        return massStr, stressStr, compStr
    def evalObjCon(self, x):
        """
        Return the objective, constraint and fail flag
        """
        self.problem.solveStructuralProblem(x[:])
        
        massKey, stressKey, compKey = self.getNames()
        
        if (self.start):
        	self.start = False
        	self.maxStress = self.problem.func[stressKey]

        fail = 0
        obj = self.problem.func[massKey] #mass
	
        con = np.zeros(1, dtype=ParOpt.dtype)
        con[0] = 1 - self.problem.func[stressKey] / self.maxStress

        self.objs.append(obj)

        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """
        #run the solver
        self.problem.solveStructuralProblem(x[:])
        
        massKey, stressKey, compKey = self.getNames()
        
        fail = 0
        g[:] = self.problem.grad[compKey]
        A[0][:] = -1 * self.problem.grad[stressKey]
        
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
        plt.ylabel('mass obj')
        plt.show()

## Optimization problem defined here ##

#run options: check, run, eval
option = "eval"

myOpt = Optimization()


if (option == "check"):
    def mass(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey, compKey = myOpt.getNames()
        return myOpt.problem.func[massKey]
    def stress(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey, compKey = myOpt.getNames()
        return myOpt.problem.func[stressKey]
    def compliance(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey, compKey = myOpt.getNames()
        return myOpt.problem.func[compKey]
    def massGrad(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey, compKey = myOpt.getNames()
        return myOpt.problem.grad[massKey]
    def stressGrad(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey, compKey = myOpt.getNames()
        return myOpt.problem.grad[stressKey]
    def complianceGrad(x):
        myOpt.problem.solveStructuralProblem(x[:])
        massKey, stressKey, compKey = myOpt.getNames()
        return myOpt.problem.grad[compKey]
    
    D = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0, 0.1, 0.1, 0.1 ]
    funcs = [mass,stress,compliance]
    gradients = [massGrad,stressGrad,complianceGrad]
    names = ["mass","stress","compl"]
    myOpt.problem.checkGradients(D,funcs,gradients,names, h=1e-4)

elif (option == "run"):
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
    D, z, zw, zl, zu = opt.getOptimizedPoint()
    
    #print optimized information
    myOpt.printObjCon(D)
    print("\n")
    print("Final design: ")
    myOpt.problem.printDesignVariables(D[:])

    myOpt.plotObjectiveProgress()

elif (option == "eval"):
    D = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0, 0.1 ,0.1, 0.1 ]
    p = np.random.uniform(size=6)
    p = p / np.linalg.norm(p)
    h = 1.0e-5
    print(p)
    print(p*h)
    myOpt.problem.solveStructuralProblem(D)
