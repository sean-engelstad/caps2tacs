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
    egadsAim.input.Edge_Point_Min = 10
    egadsAim.input.Edge_Point_Max = 15
    
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
    
    #design variables
    tacsAim.input.Design_Variable = {"area" : {},
                                     "aspect" : {},
                                     "taper" : {},
                                     "twist" : {},
                                     "lesweep" : {},
                                     "dihedral" : {}}

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

    #exit(0)
    checkTACSgradient = False
    if (checkTACSgradient):
        option = 1
        if (option == 1):
            n = 20
            oMag = np.linspace(-3,-10,n)
            hvec = 10**oMag
            massFD = np.zeros((n))
            stressFD = np.zeros((n))
            compFD = np.zeros((n))
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
                            elif ("failure" in funcName):
                                stress_error = abs_error
                            elif ("compl" in funcName):
                                comp_error = abs_error
                massFD[hind] = mass_error
                stressFD[hind] = stress_error
                compFD[hind] = comp_error
                hind += 1
            
            #plot result
            plt.subplot(3,1,1)
            plt.loglog(hvec,massFD,'k-')
            plt.xlabel('step size h')
            plt.ylabel('mass FD error')
            plt.subplot(3,1,2)
            plt.loglog(hvec,stressFD,'k-')
            plt.xlabel('step size h')
            plt.ylabel('stress FD error')
            plt.subplot(3,1,3)
            plt.loglog(hvec,compFD,'k-')
            plt.xlabel('step size h')
            plt.ylabel('compliance FD error')
            plt.show()
            exit()
        elif (option == 2):
            dh = 1e-6

            for caseID in SPs:            
                problem = SPs[caseID]

            funcs = {}
            #gradient check option 2 - component wise gradients of tacs
            x_orig = problem.getDesignVars()
            
            # Reset design variables
            problem.setDesignVars(x_orig)
            
            # Perform a fd/cs sensisitivity check on nodal coordinate sensitivity
            xpts_orig = problem.getNodes()
            
            #solve original func
            problem.setNodes(xpts_orig)
            problem.solve()
            problem.evalFunctions(funcs)
            
            nnodes = obj.FEASolver.getNumOwnedNodes()

            fdGrad = {}
            for funcName in funcs:
                fdGrad[funcName] = np.zeros((3*nnodes))

            for direc in range(3*nnodes):
                ei = np.zeros((3*nnodes))
                ei[direc] = 1
                
                # Get number of nodes owned by this proc
                
                xpts_new = xpts_orig + ei * dh
                
                funcs_new = {}

                # Solve np.mean(g)new problem (with slight shift)
                problem.setNodes(xpts_new)
                problem.solve()
                problem.evalFunctions(funcs_new)

                for funcName in funcs:
                    fdGrad[funcName][direc] = (funcs_new[funcName] - funcs[funcName])/dh



            for funcName in funcs:
                #error in gradients
                diffGrad = fdGrad[funcName] - sens[funcName]['Xpts']
                relError = np.linalg.norm(diffGrad)/np.linalg.norm(fdGrad[funcName])
                sensMag = np.linalg.norm(sens[funcName]['Xpts'])
                fdMag = np.linalg.norm(fdGrad[funcName])
                print(funcName + " fd Grad mag ", fdMag)
                print(funcName + " chainRule Grad mag ", sensMag)
                print(funcName + " total ei fd error ", error)
        #exit()

#ParOpt Optimization Class
class Optimization(ParOpt.Problem):
    def __init__(self):
        desvarList = ["area","aspect","taper","twist","lesweep","dihedral"]
        self.problem = Caps2Tacs("nacaFD_wing.csm", capsFunction, pytacsFunction, desvarList)
        
        self.nvar = 6 #number of design var
        ncon = 0 #number of constraint
        nblock = 1
        super(Optimization, self).__init__(MPI.COMM_SELF, self.nvar, ncon, nblock)

        self.objs = []
    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        #area, aspectRatio, taper, twistAngle, leadingEdgeSweep, dihedral
        lowerBounds =   [20.0, 3.0,  0.3,  1.0,  3.0,  1.0 ]
        initialValues = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0 ]
        upperBounds =   [100.0,10.0, 1.0, 10.0, 50.0, 20.0]
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

        fail = 0
        obj = self.problem.func[compKey] #mass

        #maxConstr = 2000.0 - self.problem.func[stressKey]
        con = [] #vmstress

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
        
        #A[0][:] = -1 * self.problem.grad[stressKey]
        
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
        plt.ylabel('compliance obj')
        plt.show()

## Optimization problem defined here ##

#run options: check, run, eval
option = "check"

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
    
    D = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0 ]
    funcs = [mass,stress,compliance]
    gradients = [massGrad,stressGrad,complianceGrad]
    names = ["mass","stress","compl"]
    myOpt.problem.checkGradients(D,funcs,gradients,names, h=1e-4)

elif (option == "run"):
    options = {
    'algorithm': 'mma',
    'mma_init_asymptote_offset': 0.5,
    'mma_min_asymptote_offset': 0.01,
    'mma_max_iterations': 20}
    
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
    D = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0 ]
    p = np.random.uniform(size=6)
    p = p / np.linalg.norm(p)
    h = 1.0e-5
    print(p)
    print(p*h)
    myOpt.problem.solveStructuralProblem(D+p*h)