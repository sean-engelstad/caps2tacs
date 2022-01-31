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
    
    #put in the design variables
    tacsAim.input.Design_Variable = {"area" : {},
                                     "aspect" : {},
                                     "taper" : {},
                                     "twist" : {},
                                     "lesweep" : {},
                                     "dihedral" : {},
                                     "thickRatio" : {}}

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
    
    checkTACSgradient = False
    if (checkTACSgradient):
        n = 100
        oMag = np.linspace(-3,-10,n)
        hvec = 10**oMag
        massFD = np.zeros((n))
        stressFD = np.zeros((n))
        hind = 0
        for dh in hvec:
            if (hind == 0):
                funcs = {}
            funcs_new = {}
            #design variable check
            for caseID in SPs:            
                problem = SPs[caseID]
                
                if (hind == 0):
                    #gradient check section
                    x_orig = problem.getDesignVars()
                    
                    # Reset design variables
                    problem.setDesignVars(x_orig)
                    
                    # Perform a fd/cs sensisitivity check on nodal coordinate sensitivity
                    xpts_orig = problem.getNodes()
                    
                    #solve original
                    problem.setNodes(xpts_orig)
                    problem.solve()
                    problem.evalFunctions(funcs)
                    
                    nnodes = obj.FEASolver.getNumOwnedNodes()
                
                # Get number of nodes owned by this proc
                
                xpts_pert = np.random.rand(3 * nnodes)
                xpts_pert = xpts_pert / np.linalg.norm(xpts_pert)
                
                xpts_new = xpts_orig + xpts_pert * dh
                
                
                
                # Solve np.mean(g)new problem (with slight shift)
                problem.setNodes(xpts_new)
                problem.solve()
                problem.evalFunctions(funcs_new)
    
                # Loop through each function and compare sensitivities
                for funcName in funcs:
                    print("old", funcs)
                    print("new", funcs_new)
                    dfunc_approx = np.real((funcs_new[funcName] - funcs[funcName]) / dh)
                    # Project sensitivity against perturbation vector
                    dfunc_exact_local = np.real(np.dot(sens[funcName]['Xpts'], xpts_pert))
                    # The sens vector is distributed across multiple processors,
                    # accumulate sensitivity contribution from each processor to get total sensitivity
                    dfunc_exact = comm.allreduce(dfunc_exact_local)
                    if comm.rank == 0:
                        print('Func name:      ', funcName)
                        print('FD (XptSens):      ', dfunc_approx)
                        print('Result (XptSens):  ', dfunc_exact)
                        print('Rel err (XptSens): ', (dfunc_exact - dfunc_approx) / dfunc_exact)
                        rel_error = (dfunc_exact - dfunc_approx) / dfunc_exact
                        abs_error = np.abs(rel_error)
                        if ("mass" in funcName):
                            mass_error = abs_error
                        else:
                            stress_error = abs_error
            massFD[hind] = mass_error
            stressFD[hind] = stress_error
            hind += 1
        
        #plot result
        plt.subplot(2,1,1)
        plt.loglog(hvec,massFD,'k-')
        plt.xlabel('step size h')
        plt.ylabel('mass FD error')
        plt.subplot(2,1,2)
        plt.loglog(hvec,stressFD,'k-')
        plt.xlabel('step size h')
        plt.ylabel('stress FD error')
        plt.show()
        exit()


#ParOpt Optimization Class
class Optimization(ParOpt.Problem):
    def __init__(self):
        #design variable order used for parOpt, and this will build the design dict
        desvarList = ["area","aspect","taper","twist","lesweep","dihedral","thickRatio"]
        
        self.problem = Caps2Tacs("rect_wing.csm", capsFunction, pytacsFunction, desvarList)
        
        self.nvar = 7 #number of design var
        ncon = 1
        nblock = 1
        super(Optimization, self).__init__(MPI.COMM_SELF, self.nvar, ncon, nblock)
        
        self.objs = []
    def getVarsAndBounds(self, x, lb, ub):
        #"""Get the variable values and bounds"""
        #area, aspectRatio, taper, twistAngle, leadingEdgeSweep, dihedral, thicknessRatio
        lowerBounds =   [20.0, 3.0,  0.3,  1.0,  3.0,  1.0, 0.05 ]
        initialValues = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0, 0.12 ]
        upperBounds =   [100.0,10.0, 1.0, 10.0, 50.0, 10.0, 0.20]
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

        maxStress = 0.35
        maxConstr = maxStress - self.problem.func[stressKey]
        con = [maxConstr] #vmstress
        
        self.objs.append(obj)

        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """desDict
        Return the objective, constraint and fail flag
        """        
        #run the solver, x[:] opens the parOpt object to an array
        self.problem.solveStructuralProblem(x[:])
        
        massKey, stressKey = self.getNames()
        
        fail = 0
        g[:] = self.problem.grad[massKey]
        
        A[0][:] = -self.problem.grad[stressKey]
        
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

#run options: check, run
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
    
    x = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0, 0.12 ]
    funcs = [mass,stress]
    gradients = [massGrad,stressGrad]
    names = ["mass","stress"]
    myOpt.problem.checkGradients(x,funcs,gradients,names,h=1e-6)
    
elif (option == "run"):    
    options = {
    'algorithm': 'mma',
    'mma_init_asymptote_offset': 0.5,
    'mma_min_asymptote_offset': 0.01,
    'mma_max_iterations': 10}
    
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
    
    myOpt.plotObjectiveProgress()
