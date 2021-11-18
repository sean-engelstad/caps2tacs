#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 16:52:03 2021

@author: sengelstad6
"""

from __future__ import print_function
import unittest
import time

import os
import glob
import shutil

import sys

import pyCAPS

####################################################################

filename = "panel2"
csmFile = os.path.join("./CSM",filename + ".csm")

myProblem = pyCAPS.Problem('myCAPS', capsFile=csmFile, outLevel=0)

geom = myProblem.geometry

egads = myProblem.analysis.create(aim="egadsTessAIM")

tacsAnalysis = myProblem.analysis.create(aim = "tacsAIM",
                                 name = "tacs")

egads.input.Edge_Point_Min = 15
egads.input.Edge_Point_Max = 20
        
egads.input.Mesh_Elements = "Quad"
        
egads.input.Tess_Params = [.25,.01,15]
        
NumberOfNode = egads.output.NumberOfNode
        
# Link the mesh
tacsAnalysis.input["Mesh"].link(egads.output["Surface_Mesh"])
        
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
        
tacsAnalysis.input.Property = {"leftPlate": shell,
                               "rightPlate": shell}
# Set constraints
constraint = {"groupName" : "edge",
              "dofConstraint" : 123456}
        
tacsAnalysis.input.Constraint = {"edgeConstraint": constraint}

# Set load
load = {"groupName" : "plate",
        "loadType" : "Pressure",
        "pressureForce" : 2.e6}
        
# Set loads
tacsAnalysis.input.Load = {"appliedPressure": load }

tacsAnalysis.input.Design_Variable = {"Length" : {},
                              "Width"  : {}}

tacsAnalysis.preAnalysis()

homeDir = os.getcwd()
csmDir = os.path.join(homeDir,"BDF")
bdfFile = os.path.join(tacsAnalysis.analysisDir, tacsAnalysis.input.Proj_Name + '.bdf')
datFile = os.path.join(tacsAnalysis.analysisDir, tacsAnalysis.input.Proj_Name + '.dat')

bdfFile2 = os.path.join(csmDir, filename + '.bdf')
datFile2 = os.path.join(csmDir, filename + '.dat')

shutil.move(bdfFile,bdfFile2)
shutil.move(datFile,datFile2)
