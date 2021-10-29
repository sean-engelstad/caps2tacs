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

class TestTACS(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.problemName = "workDir_tacs"
        cls.iProb = 1
        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

    def test_Plate(self): #probably want to make a function df/dX(desvarX) which runs TACS each time, then put that into parOpt

        filename = os.path.join("stiffPanel.csm")
        #mkdir(self.problemName + str(self.iProb))
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=filename, outLevel=0); self.__class__.iProb += 1

        mesh = myProblem.analysis.create(aim="egadsTessAIM")
        
        tacs = myProblem.analysis.create(aim = "tacsAIM",
                                         name = "tacs")

        mesh.input.Edge_Point_Min = 3
        mesh.input.Edge_Point_Max = 4

        mesh.input.Mesh_Elements = "Quad"

        mesh.input.Tess_Params = [.25,.01,15]
        
        mesh.runAnalysis()
        
        mesh.geometry.view()
        
        # Link the mesh
        tacs.input["Mesh"].link(mesh.output["Surface_Mesh"])

        # Set analysis type
        tacs.input.Analysis_Type = "Static"

        # Set materials
        madeupium    = {"materialType" : "isotropic",
                        "youngModulus" : 72.0E9 ,
                        "poissonRatio": 0.33,
                        "density" : 2.8E3}

        tacs.input.Material = {"Madeupium": madeupium}

        # Set properties
        shell  = {"propertyType" : "Shell",
                  "membraneThickness" : 0.006,
                  "material"        : "madeupium",
                  "bendingInertiaRatio" : 1.0, # Default
                  "shearMembraneRatio"  : 5.0/6.0} # Default

        tacs.input.Property = {"plate": shell}

        # Set constraints
        constraint = {"groupName" : "plateEdge",
                      "dofConstraint" : 123456}

        tacs.input.Constraint = {"edgeConstraint": constraint}

        # Set load
        load = {"groupName" : "plate",
                "loadType" : "Pressure",
                "pressureForce" : 2.e6}

        # Set loads
        tacs.input.Load = {"appliedPressure": load }

        tacs.input.Design_Variable = {"plateLength" : {},
                                      "plateWidth"  : {},
                                      "stiffHeight" : {}}

        # Run Small format
        tacs.preAnalysis()
        
        tacs.geometry.view()
        time.sleep(20)
        
        # Create a dummy sensitivity file
        filename = os.path.join(tacs.analysisDir, tacs.input.Proj_Name+".sens")
        with open(filename, "w") as f:
            f.write("2 {}\n".format(NumberOfNode)) # Two functionals Number Of Nodes
            f.write("Func1\n") # 1st functional
            f.write("42\n")    # Value of Func1
            for i in range(NumberOfNode): # d(Func1)/d(xyz)
                f.write("{} {} {}\n".format(i, 2*i, 3*i))

            f.write("Func2\n") # 2nd functiona;
            f.write("21\n")    # Value of Func2
            for i in range(NumberOfNode): # d(Func2)/d(xyz)
                f.write("{} {} {}\n".format(3*i, 2*i, i))

        tacs.postAnalysis()

        Func1             = tacs.dynout["Func1"].value
        Func1_plateLength = tacs.dynout["Func1"].deriv("plateLength")
        Func1_plateWidth  = tacs.dynout["Func1"].deriv("plateWidth")
        
        self.assertEqual(Func1, 42)

        Func2             = tacs.dynout["Func2"].value
        Func2_plateLength = tacs.dynout["Func2"].deriv("plateLength")
        Func2_plateWidth  = tacs.dynout["Func2"].deriv("plateWidth")

        self.assertEqual(Func2, 21)

        print(filename)

if __name__ == '__main__':
    unittest.main()
