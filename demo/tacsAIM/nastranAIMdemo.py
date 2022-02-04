# -*- coding: utf-8 -*-

import os
import pyCAPS, argparse
import pytacs, pyNastran
from tacs import TACS, elements, constitutive


### SET DIRECTORIES ################
dir_root = os.getcwd()
dir_csm = os.path.join(dir_root,"CSM")

# START CAPS PROBLEM ========================================
csm_file = "f118.csm"
myCAPS = pyCAPS.Problem(problemName="WorkDir",
                               capsFile = csm_file,
                               outLevel = 1)
geometry = myCAPS.geometry
 
### MODIFY DESIGN PARAMETERS ############
#tips, geometry.despmtr[key][key].value;
#geometry.view()

geometry.despmtr["wing"].area = 2*geometry.despmtr["wing"].area
geometry.despmtr["wing"].aspect = 0.5 * geometry.despmtr["wing"].aspect
   
#NASTRAN AIM ================================================
#load the Nastran Aim
nastranAIM = myCAPS.analysis.create(aim = "nastranAIM",
                                    name = "nastran")
#AIM inputs
projectName = "f118"
nastranAIM.input.Proj_Name = projectName
nastranAIM.input.File_Format = "Free"
nastranAIM.input.Mesh_File_Format = "Large"
nastranAIM.input.Edge_Point_Max = 2
nastranAIM.input.Edge_Point_Min = 2
nastranAIM.input.Analysis_Type = "Static"

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.0E7 ,
                "poissonRatio" : .33,
                "density"      : 0.1}
 
nastranAIM.input.Material = {"Madeupium": madeupium}

shell1  =   {"propertyType"      : "Shell",
          "material"          : "Madeupium",
          "membraneThickness"      : 1.0}
 
shell2  =   {"propertyType"     : "Shell",
          "material"          : "Madeupium",
          "membraneThickness"      : 2.0}
 
nastranAIM.input.Property = {"wing:faces": shell1,
                             "fuselage:faces": shell2,
                             "vtail:faces": shell1,
                             "htail:faces": shell1}
constraint = {"groupName"         : "boundary",
              "dofConstraint"     : 123456}
 
nastranAIM.input.Constraint = {"BoundaryCondition": constraint}

load = {"groupName"         : "force",
        "loadType"             : "GridForce",
        "forceScaleFactor"     : 20000.0,
        "directionVector"     : [0.8, -0.6, 0.0]}
 
nastranAIM.input.Load = {"appliedForce": load}

value = {"analysisType"         : "Static",
         "analysisConstraint"     : "BoundaryCondition",
         "analysisLoad"         : "appliedForce"}
 
myCAPS.analysis["nastran"].input.Analysis = {"SingleLoadCase": value}

nastranAIM.preAnalysis()

print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory
os.chdir(nastranAIM.analysisDir) # Move into test directory
#if (args.noAnalysis == False): os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call
os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")

nastranAIM.postAnalysis()

nastranAIM.geometry.view()

