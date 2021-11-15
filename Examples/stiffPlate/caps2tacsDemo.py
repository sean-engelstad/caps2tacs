#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 16:52:03 2021

@author: sengelstad6
"""

import os

import pyCAPS

from tacs.pytacs import pyTACS

from tacs import functions
#probably want to make a function df/dX(desvarX) which runs TACS each time, then put that into parOpt

#filename = os.path.join("stiffPanel.csm")
filename = os.path.join("stiffPanel2.csm")

#mkdir(self.problemName + str(self.iProb))
myProblem = pyCAPS.Problem('myCAPS', capsFile=filename, outLevel=0)

#myProblem.geometry.view()

mesh = myProblem.analysis.create(aim="egadsTessAIM")

tacsAnalysis = myProblem.analysis.create(aim = "tacsAIM",
                                 name = "tacs")

mesh.input.Edge_Point_Min = 15
mesh.input.Edge_Point_Max = 20

mesh.input.Mesh_Elements = "Quad"

mesh.input.Tess_Params = [.25,.01,15]

NumberOfNode = mesh.output.NumberOfNode

# Link the mesh
tacsAnalysis.input["Mesh"].link(mesh.output["Surface_Mesh"])

# Set analysis type
tacsAnalysis.input.Analysis_Type = "Static"

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3,
                "tensionAllow" :  20.0e7}

tacsAnalysis.input.Material = {"Madeupium": madeupium}

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default
#need to add the strength here

tacsAnalysis.input.Property = {"plate": shell,
                       "stiffener": shell}

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

tacsAnalysis.input.Constraint = {"edgeConstraint": constraint}

# Set load
load = {"groupName" : "plate",
        "loadType" : "Pressure",
        "pressureForce" : 2.e6}

# Set loads
tacsAnalysis.input.Load = {"appliedPressure": load }

tacsAnalysis.input.Design_Variable = {"plateLength" : {},
                              "plateWidth"  : {},
                              "stiffHeight" : {}}

# Run Small format
tacsAnalysis.preAnalysis()

#tacs.geometry.view()



#pytacs
structOptions = {'writeSolution': True, }

caps_dir = os.path.join(os.getcwd(),'myCAPS/Scratch/tacs')
datFile = os.path.join(caps_dir, 'nastran_CAPS.dat')

# Load BDF file
FEASolver = pyTACS(datFile, options=structOptions)
# Set up TACS Assembler
FEASolver.createTACSAssembler()
#add functions
FEASolver.addFunction('wing_mass', functions.StructuralMass)
FEASolver.addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5,
              KSWeight=100.0)
evalFuncs = ['wing_mass', 'ks_vmfailure']

# Read in forces from BDF and create tacs struct problems
SPs = FEASolver.createTACSProbsFromBDF()
#print("ran pytacs")
# Solve each structural problem and write solutions
funcs = {}; funcsSens = {}
for caseID in SPs:
    FEASolver(SPs[caseID])
    FEASolver.evalFunctions(SPs[caseID], funcs,evalFuncs=evalFuncs)
    FEASolver.evalFunctionsSens(SPs[caseID], funcsSens,evalFuncs=evalFuncs)
    FEASolver.writeSolution(outputDir=os.path.dirname(__file__))

funcKeys = funcs.keys()
nfunc = len(funcKeys)

dfdX = funcsSens['load_set_001_wing_mass']['Xpts']


#reorder nodes, start at
TACSnodeMap = FEASolver._getGlobalToLocalNodeIDDict() #error in nodemap, less nodes
#print(TACSnodeMap)

#nnodes = int(len(dfdX)/3)
#dfdX = dfdX.reshape(nnodes,3)

#dfdX_bdf = np.zeros((nnodes,3))
#print(len(TACSnodeMap))
#for i in range(nnodes):
#    #i is bdf node, tacsNode is globalNode
#    tacsNode = TACSnodeMap[i]
#    dfdX_bdf[i,:] = dfdX[tacsNode,:]

filename = os.path.join(tacsAnalysis.analysisDir, tacsAnalysis.input.Proj_Name+".sens")
with open(filename, "w") as f:
    f.write("{} {}\n".format(nfunc,NumberOfNode))
    for key in funcKeys:
        cSens = funcsSens[key]['Xpts']
        cSens = cSens.reshape(NumberOfNode,3)
        f.write(key + "\n")
        f.write("{}\n".format(funcs[key]))
        for nodeind in range(NumberOfNode): # d(Func1)/d(xyz)
            f.write("{} {} {}\n".format(cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))

tacsAnalysis.postAnalysis()

#print(dir(tacsAnalysis))

desKeys = tacsAnalysis.input.Design_Variable.keys()
ndesVar = len(desKeys)
func = {}
sens = {}
for key in funcKeys:
    func[key] = tacsAnalysis.dynout[key].value
    sens[key] = {}
    for desKey in desKeys:
        sens[key][desKey] = tacsAnalysis.dynout[key].deriv(desKey)
#could convert func and sens to array/vec instead of nested dictionary
print("Functions: ",func)
print("Sensitivities: ",sens)