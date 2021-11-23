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

filename = "stiffPanel4"
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

propertyDict = {}
for i in range(80):
    propertyDict["plate"+str(i+1)] = shell
for j in range(32):
    propertyDict["stiffener"+str(j+1)] = shell
tacsAnalysis.input.Property = propertyDict

# Set constraints
roller = {"groupName" : "edge1",
              "dofConstraint" : 12}
pin = {"groupName" : "edge2",
              "dofConstraint" : 123}
        
tacsAnalysis.input.Constraint = {"roller": roller,
                                 "pin": pin}

# Set load, specifying total load, not line load here
loadValue = 20
numNodesOnEdge = 161
loadPerNode = loadValue / numNodesOnEdge
load = {"groupName"         : "edge1",
        "loadType"             : "GridForce",
        "forceScaleFactor"     : loadPerNode,
        "directionVector"     : [0.0, 0.0, 1.0]}
        
# Set loads
tacsAnalysis.input.Load = {"appliedPressure": load }

tacsAnalysis.input.Design_Variable = {"plateLength" : {},
                              "plateWidth"  : {},
                              "stiffHeight" : {}}

tacsAnalysis.preAnalysis()

#tacsAnalysis.geometry.view()

homeDir = os.getcwd()
csmDir = os.path.join(homeDir,"BDF")
bdfFile = os.path.join(tacsAnalysis.analysisDir, tacsAnalysis.input.Proj_Name + '.bdf')
datFile = os.path.join(tacsAnalysis.analysisDir, tacsAnalysis.input.Proj_Name + '.dat')

bdfFile2 = os.path.join(csmDir, filename + '.bdf')
datFile2 = os.path.join(csmDir, filename + '.dat')

shutil.move(bdfFile,bdfFile2)
shutil.move(datFile,datFile2)
