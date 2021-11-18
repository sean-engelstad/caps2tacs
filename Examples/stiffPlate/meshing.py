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

self.csmFile = csmFile

self.myProblem = pyCAPS.Problem('myCAPS', capsFile=csmFile, outLevel=0)

self.geom = self.myProblem.geometry
self.deskeys = self.geom.despmtr.keys()

self.egads = self.myProblem.analysis.create(aim="egadsTessAIM")

self.tacs = self.myProblem.analysis.create(aim = "tacsAIM",
                                 name = "tacs")

#filename = os.path.join("stiffPanel.csm")
filename = os.path.join("panel2.csm")
#mkdir(self.problemName + str(self.iProb))
myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=filename, outLevel=0); self.__class__.iProb += 1

egads = myProblem.analysis.create(aim="egadsTessAIM")

tacsAnalysis = myProblem.analysis.create(aim = "tacsAIM",
                                 name = "tacs")

egads.input.Edge_Point_Min = 10
egads.input.Edge_Point_Max = 15

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

# Run Small format
tacsAnalysis.preAnalysis()

tacsAnalysis.geometry.view()

#time.sleep(20)

bdf_file = 'nastran_CAPS.bdf'
dat_file = 'nastran_CAPS.dat'

orig_dir = os.getcwd()

caps_dir = os.path.join(orig_dir,'myCAPS1/Scratch/tacs')

bdf_dir = os.path.join(caps_dir,bdf_file)
dat_dir = os.path.join(caps_dir,dat_file)

bdf_newdir = os.path.join(orig_dir,bdf_file)
data_newdir = os.path.join(orig_dir,dat_file)

shutil.move(bdf_dir,bdf_newdir)
shutil.move(dat_dir,data_newdir)
