# -*- coding: utf-8 -*-
import os
import pyCAPS
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions
import numpy as np
import matplotlib.pyplot as plt

##-----------Initialize PyCAPS----------##

#initialize pycaps
caps = pyCAPS.Problem('myCAPS', capsFile="naca_small.csm", outLevel=0)

#initialize egads Aim
egadsAim = caps.analysis.create(aim="egadsTessAIM")

#initialize tacs Aim
tacsAim = caps.analysis.create(aim = "tacsAIM", name = "tacs")

##-----------PyCAPS setup function-------##
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

##-----------setup DVs and DVRs-----------##
#setup the Design_Variable and Design_Variable_Relations

#list of design variables, with thick1-thickN for thickness DVs
desvars = ["area","aspect","taper","ctwist","lesweep","dihedral","thick1", "thick2", "thick3"]
nvar = len(desvars)

#where the thickDVs have thick1 is for "rib", thick2 for "spar" etc.
capsGroups = ["rib","spar","OML"]
thickness = [0.01, 0.02, 0.03]

def makeThicknessDV(capsGroup, thickness):
    #thick DV dictionary for Design_Variable Dict
    desvar    = {"groupName" : capsGroup,
            "initialValue" : thickness,
            "lowerBound" : thickness*0.5,
            "upperBound" : thickness*1.5,
            "maxDelta"   : thickness*0.1}
    return desvar
    
def makeThicknessDVR(DVname):
    #thick DV dictionary for Design_Variable_Relation Dict
    DVR = {"variableType": "Property",
    "fieldName" : "T",
    "constantCoeff" : 0.0,
    "groupName" : DVname,
    "linearCoeff" : 1.0}
    return DVR

#make initial DV and DVRdict
DVdict = {}
DVRdict = {}

#add thickDVs and geomDVs to caps
thickCt = 0
for desvar in desvars:
    if ("thick" in desvar):
        #add thickDV entry into DV_Relations and DV Dicts
        DVRdict[desvar] = makeThicknessDVR(desvar)
        DVdict[desvar] = makeThicknessDV(capsGroups[thickCt],thickness[thickCt])
        thickCt += 1
    else: #geomDV, add empty entry into DV dicts
        DVdict[desvar] = {}

#input DVdict and DVRdict into tacsAim
tacsAim.input.Design_Variable = DVdict
tacsAim.input.Design_Variable_Relation = DVRdict

##---------PYTACS to run TACS-------------##
#first run the preanalysis to prepare to run tacs through pytacs
tacsAim.preAnalysis()

#setup MPI COMM for pytacs
comm = MPI.COMM_WORLD

#data file
datFile = os.path.join(tacsAim.analysisDir, tacsAim.input.Proj_Name + '.dat')

#initialize pytacs with that data file
FEASolver = pyTACS(datFile)
    
# Set up TACS Assembler
FEASolver.initialize()

#choose the functions to evaluate
evalFuncs = ['wing_mass', 'ks_vmfailure']

#read the bdf & dat file into pytacs FEAsolver
#SPs represents "StructuralProblems"
SPs = FEASolver.createTACSProbsFromBDF()

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
    #print("finished pytacs funcs")
    SPs[caseID].evalFunctionsSens(sens,evalFuncs=evalFuncs)
    #print("finished pytacs sens")
    SPs[caseID].writeSolution(outputDir=os.path.dirname(__file__))
    #print("finished pytacs file")


##-----------------Reorder the nodes-------------------##
#get bookkeeping info on func such as funcKeys, nfunctions, nnodes
funcKeys = func.keys()
nfunc = len(funcKeys)
funcKeysList = list(funcKeys)

#get number of nodes
subSens = sens[funcKeysList[0]]['Xpts']
nnodes = int(len(subSens)/3)
#print(self.nnodes)

for key in funcKeys:
            
    #get sensitivity w.r.t. mesh aka Xpts sens
    dfdX = sens[key]['Xpts']
    dfdX = dfdX.reshape((nnodes,3))
    
    #thickness desvar gradient, no reordering for this
    #commented out here since don't need this yet
    #dfdH = sens[key]['struct']
    
    #get the nodemap from pytacs
    bdfNodes = np.linspace(1.0,nnodes,nnodes)
    tacsNodeMap = FEASolver.meshLoader.getLocalNodeIDsFromGlobal(bdfNodes,nastranOrdering=True)
    
    #reorder nodes as dfdX_bdf
    dfdX_bdf = np.zeros((nnodes,3))
    for bdfind in range(nnodes):
            #node 1 in bdf, corresponds to 0 bdfind, and plug in bdfind-1 to get tacsInd cause cyclically off by 1
            tacsInd = tacsNodeMap[bdfind] #tacsNodeMap is cyclically off by 1
            dfdX_bdf[bdfind,:] = dfdX[tacsInd,:]
    
    #update function sensitivity
    sens[key]['Xpts'] = dfdX_bdf

##---------Post analysis/print sensitivity to CAPS-----#
#print our dfdX sensitivity to CAPS
#then this will be combined to produce df/dD sensitivity
#available in dynout variables

#where to print .sens file
sensFile = os.path.join(tacsAim.analysisDir, tacsAim.input.Proj_Name+".sens")

#open the file
with open(sensFile, "w") as f:
    
    #write (nfunctions) in first line
    f.write("{}\n".format(nfunc))
    
    #for each function mass, stress, etc.
    for key in funcKeys:
        
        #get the pytacs/tacs sensitivity w.r.t. mesh for that function
        cSens = sens[key]['Xpts']
        
        #write the key,value,nnodes of the function
        f.write(key + "\n")
        f.write("{}\n".format(func[key]))
        f.write("{}\n".format(nnodes))

        #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
        for nodeind in range(nnodes): # d(Func1)/d(xyz)
            
            printNodes = True
            #if we print nodes first, then print them
            if (printNodes):
                bdfind = nodeind + 1
                f.write("{} {} {} {}\n".format(bdfind, cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
            #otherwise just print df/dx, df/dy, df/dz for each element
            else:
                f.write("{} {} {}\n".format(cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))

#run aim postanalysis
tacsAim.postAnalysis()
print("ran postAnalysis()")

##-------------Store Gradient from CAPS dynout---------##
#store the function and full df/dD sensitivities from CAPS aim dynout attributes
#initialize gradient variable
func2 = {} #func2 matches original func
grad = {}
print("starting store Results")
#loop over each pytacs function
for key in funcKeys:
    print("starting function: ",key)
    func2[key] = tacsAim.dynout[key].value
    grad[key] = np.zeros((nvar))
    #print(len(self.grad[key]))

    #loop over each design variable to get the full df/dD gradient
    print("struct sens: ",sens[key]['struct'])
    ind = 0
    thickind = 0
    for desvar in desvars:
        print("storing desvar: ",key, ", ", desvar)
        if ("thick" in desvar):
            #use struct here, #print(self.sens[key]['struct'])
            #struct includes geomDVs in same order
            grad[key][ind] = sens[key]['struct'][ind]
            thickind += 1
        else:
            print(tacsAim.dynout[key].deriv(desvar))
            grad[key][ind] = tacsAim.dynout[key].deriv(desvar)
        ind += 1
#print("finished storing results\n")
print("gradient: ",grad)

##--------Print gradient of each function for demo------##
print("----Results----")
for key in funcKeys:
    print("function {}, value {:.4f}".format(key, func2[key]))
    ind = 0
    for desvar in desvars:
        print("\t deriv for DV {} = {:.4f}".format(desvar, grad[key][ind]))
        ind += 1