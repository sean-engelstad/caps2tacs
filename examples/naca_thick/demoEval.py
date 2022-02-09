# -*- coding: utf-8 -*-
import os
from caps2tacs import Caps2Tacs
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions
import numpy as np
import matplotlib.pyplot as plt

def makeThicknessDV(capsGroup, thickness):
    desvar    = {"groupName" : capsGroup,
          "initialValue" : thickness,
          "lowerBound" : thickness*0.5,
          "upperBound" : thickness*1.5,
          "maxDelta"   : thickness*0.1,
          "fieldName" : "T"}
    return desvar

def makeThicknessDVR(capsGroup):
    DVR = {"variableType": "Property",
    "fieldName" : "T",
    "constantCoeff" : 0.0,
    "groupName" : capsGroup,
    "linearCoeff" : 1.0}
    return DVR

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

    useDVR = True
    for i in range(3):
        #make Design_Variable_Relations
        if (useDVR):
            DVRdict[thickDVs[i]] = makeThicknessDVR(capsGroups[i])
        
        #also make Design_Variables
        DVdict[thickDVs[i]] = makeThicknessDV(capsGroups[i],thick0)
        
    
    tacsAim.input.Design_Variable = DVdict   
    if useDVR: tacsAim.input.Design_Variable_Relation = DVRdict

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

D = [40.0, 6.0,  0.5,  5.0,  30.0, 5.0, 0.1 ,0.1, 0.1 ]
desvarList = ["area","aspect","taper","twist","lesweep","dihedral","thick1", "thick2", "thick3"]
problem = Caps2Tacs("naca_thick.csm", capsFunction, pytacsFunction, desvarList)
problem.solveStructuralProblem(D)
