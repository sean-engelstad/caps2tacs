# -*- coding: utf-8 -*-
import os
from caps2tacs import Caps2Tacs
from mpi4py import MPI
from tacs.pytacs import pyTACS
from scipy.optimize import minimize
from tacs import functions
import numpy as np
import matplotlib.pyplot as plt

def capsFunction(egadsAim,tacsAim):
    #setup function for naca_small.csm
    
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

    #return the capsGroups you want to have thickDVs in that order
    #[thick1, thick2, thick3]
    capsDVgroups = ["rib", "spar", "OML"]

    return capsDVgroups
        

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
       #SPs[caseID].addFunction('compliance', functions.Compliance)
    # Solve each structural problem and write solutions
    func = {}; sens = {}
    for caseID in SPs:
        SPs[caseID].solve()
        print("finished pytacs solve")
        SPs[caseID].evalFunctions(func,evalFuncs=evalFuncs)
        print("finished pytacs funcs")
        SPs[caseID].evalFunctionsSens(sens,evalFuncs=evalFuncs)
        print("finished pytacs sens")
        SPs[caseID].writeSolution(outputDir=os.path.dirname(__file__))
        print("finished pytacs file")

    #store function and sensitivity values    
    obj.func = func
    obj.sens = sens

desvarList = ["area","aspect","taper","ctwist","lesweep","dihedral","thick1", "thick2", "thick3"]
problem = Caps2Tacs("naca_small.csm", capsFunction, pytacsFunction, desvarList)

def getNames():
    #get the function names for mass, stress, and compliance
    massStr = ""
    stressStr = ""
    compStr = ""
    for key in problem.funcKeys:
        #print("Key: ",key)
        if ("mass" in key):
            massStr = key
        if ("failure" in key):
            stressStr = key
        if ("compliance" in key):
            compStr = key
    return massStr, stressStr, compStr

def mass(x):
    #solve structural problem with design variable x
    problem.solveStructuralProblem(x[:])

    #get the names of each function
    massKey, stressKey, compKey = getNames()

    #get the 
    mass = problem.func[massKey]

    return mass

def massGrad(x):
    #solve structural problem with design variable x
    problem.solveStructuralProblem(x[:])

    #get the names of each function
    massKey, stressKey, compKey = getNames()

    #get mass gradient
    massGrad = problem.grad[massKey]

    return massGrad

#initial design variable
x0 = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0, 0.03, 0.05, 0.02]
bnds = (9*[0.01],9*[100])

#minimize the mass with scipy minimize
res = minimize(mass, x0, method="BFGS", jac=massGrad, bounds=bnds,options={'disp': True})
print(res)
