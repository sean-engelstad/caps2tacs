import os
import pyCAPS, argparse
from tacs import TACS, elements, constitutive


# START CAPS PROBLEM ========================================
csm_file = "panel.csm"
myCAPS = pyCAPS.Problem(problemName="WorkDir",
                               capsFile = "panel.csm",
                               outLevel = 1)
geometry = myCAPS.geometry
   
#NASTRAN AIM ================================================
#load the Nastran Aim
#egadsAim = myCAPS.analysis.create(aim="egadsTessAIM")
nastranAIM = myCAPS.analysis.create(aim = "nastranAIM",
                                    name = "nastran")
#AIM inputs
projectName = "panel"
nastranAIM.input.Proj_Name = projectName
nastranAIM.input.File_Format = "Free"
nastranAIM.input.Mesh_File_Format = "Large"
nastranAIM.input.Edge_Point_Max = 4
nastranAIM.input.Edge_Point_Min = 4
nastranAIM.input.Analysis_Type = "Static"

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.0E7 ,
                "poissonRatio" : .33,
                "density"      : 0.1}
 
nastranAIM.input.Material = {"Madeupium": madeupium}

shell1  =   {"propertyType"      : "Shell",
          "material"          : "Madeupium",
          "membraneThickness"      : 2.0}
 
shell2  =   {"propertyType"     : "Shell",
          "material"          : "Madeupium",
          "membraneThickness"      : 2.0}
 
nastranAIM.input.Property = {"plate": shell1}
constraint = {"groupName"         : "edge",
              "dofConstraint"     : 123456}
 
nastranAIM.input.Constraint = {"BoundaryCondition": constraint}

load = {"groupName"         : "bottom",
        "loadType"             : "GridForce",
        "forceScaleFactor"     : 20000.0,
        "directionVector"     : [0.0, 0.0, 1.0]}
 
nastranAIM.input.Load = {"appliedForce": load}

value = {"analysisType"         : "Static",
         "analysisConstraint"     : "BoundaryCondition",
         "analysisLoad"         : "appliedForce"}
 
myCAPS.analysis["nastran"].input.Analysis = {"SingleLoadCase": value}

def makeThicknessDV(capsGroup, thickness):
      desvar    = {"groupName" : capsGroup,
              "initialValue" : thickness,
              "lowerBound" : thickness*0.5,
              "upperBound" : thickness*1.5,
              "maxDelta"   : thickness*0.1,
              "fieldName" : "T"}
      return desvar

#thickness design variable section
thickness = shell1["membraneThickness"]
desvar1 = makeThicknessDV("plate", thickness)

#thickness and geometry design variables don't work together yet
thicknessAndGeom = False
if (thicknessAndGeom):
    dvDict = {"thick" : desvar1,
                "plateLength" : {},
                "plateWidth" : {}}
else:
    dvDict = {"thick" : desvar1}

nastranAIM.input.Design_Variable = dvDict

### Don't seem to be able to change nastranAIM.input values after they are set
### How to change the plate thickness in future optimization iterations??
# #trying to change the initial value
# DVname = "thick"
# newThickness = 2 * thickness
# capsGroup = nastranAIM.input.Design_Variable[DVname]["groupName"]
# nastranAIM.input.Property[capsGroup]["membraneThickness"] = newThickness
# nastranAIM.input.Design_Variable[DVname] = {}


nastranAIM.preAnalysis()

print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory
os.chdir(nastranAIM.analysisDir) # Move into test directory
#if (args.noAnalysis == False): os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call
os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")

#generate the bdf file
nastranAIM.postAnalysis()

#nastranAIM.geometry.view()
