import os
import pyCAPS, argparse
from tacs import TACS, elements, constitutive
import numpy as np


# START CAPS PROBLEM ========================================
csm_file = "stiffPanel3.csm"
myCAPS = pyCAPS.Problem(problemName="WorkDir",
                               capsFile = csm_file,
                               outLevel = 1)
geometry = myCAPS.geometry
   
#NASTRAN AIM ================================================
#load the Nastran Aim
nastranAIM = myCAPS.analysis.create(aim = "nastranAIM",
                                    name = "nastran")
#AIM inputs
projectName = "stiffPanel"
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

shell =   {"propertyType"      : "Shell",
          "material"          : "Madeupium",
          "membraneThickness"      : 1.0}

propertyDict = {}
DVdict = {}
thickInd = 1

def makeThicknessDV(capsGroup, thickness):
      desvar    = {"groupName" : capsGroup,
              "initialValue" : thickness,
              "lowerBound" : thickness*0.5,
              "upperBound" : thickness*1.5,
              "maxDelta"   : thickness*0.1,
              "fieldName" : "T"}
      return desvar

nPlates = 80
nStiffs = 4
plateT = np.zeros((nPlates))
stiffT = np.zeros((nStiffs))

#iterate over plate piece capsGroups
for i in range(nPlates):
    ind = i + 1
    #set the thickness, shell property
    thick = ind/11
    plateT[i] = thick
    shell = {"propertyType"      : "Shell",
          "material"          : "Madeupium",
          "membraneThickness"      : thick}

    #store the capsGroup
    group = "plate" + str(ind)

    #add the shell propety for that capsGroup
    propertyDict[group] = shell

    #add the thicknessDV for that capsGroup
    thickness = shell["membraneThickness"]
    DVname = group
    DVdict[DVname] = makeThicknessDV(group, thick)

#iterate over each stiffener capsGroup
for j in range(nStiffs):
    ind = j + 1
    #set the thickness, shell property
    thick = ind/11
    stiffT[j] = thick
    shell = {"propertyType"      : "Shell",
          "material"          : "Madeupium",
          "membraneThickness"      : thick}

    #store the capsGroup
    group = "stiff" + str(ind)

    #add the shell propety for that capsGroup
    propertyDict[group] = shell

    #add the thicknessDV for that capsGroup
    DVname = group
    thickness = shell["membraneThickness"]
    DVdict[DVname] = makeThicknessDV(group, thick)
 
#store the propertyDict we just built
nastranAIM.input.Property = propertyDict

#add the geometric design variables
#currently doens't work when you have both thicknessAndGeom design variables
#thus I turned this off for now
thicknessAndGeom = False
if (thicknessAndGeom):
    geomDVs = ["plateLength", "plateWidth", "stiffHeight"]
    for geomDV in geomDVs:
        DVdict[geomDV] = {}

#store the design variable dictionary
nastranAIM.input.Design_Variable = DVdict

#set the boundary conditions
constraint = {"groupName"         : "edge",
              "dofConstraint"     : 123456}
nastranAIM.input.Constraint = {"BoundaryCondition": constraint}

#set the loads
load = {"groupName"         : "plate",
        "loadType"             : "GridForce",
        "forceScaleFactor"     : 20000.0,
        "directionVector"     : [0.0, 0.0, 1.0]}
nastranAIM.input.Load = {"appliedForce": load}

#set the analysis type
value = {"analysisType"         : "Static",
         "analysisConstraint"     : "BoundaryCondition",
         "analysisLoad"         : "appliedForce"}
myCAPS.analysis["nastran"].input.Analysis = {"SingleLoadCase": value}

#run pre and post analysis to generate bdf and dat files then run nastran
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
