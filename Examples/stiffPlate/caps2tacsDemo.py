#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 16:52:03 2021

@author: sengelstad6
"""

from __future__ import print_function

import unittest

import os, glob, shutil, sys

import pyCAPS

import numpy as np

from tacs.pytacs import pyTACS

from tacs import functions

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

        mesh.input.Edge_Point_Min = 10
        mesh.input.Edge_Point_Max = 20

        mesh.input.Mesh_Elements = "Quad"

        mesh.input.Tess_Params = [.25,.01,15]
        
        # Link the mesh
        tacs.input["Mesh"].link(mesh.output["Surface_Mesh"])

        # Set analysis type
        tacs.input.Analysis_Type = "Static"

        # Set materials/home/sengelstad6/git/caps2tacs/workDir_tacs1/Scratch/tacs/nastran_CAPS.bdf
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
        #need to add the strength here

        tacs.input.Property = {"plate": shell,
                               "stiffener": shell}

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
        
        #tacs.geometry.view()
        
        #time.sleep(20)
        
        bdf_file = 'nastran_CAPS.bdf'
        dat_file = 'nastran_CAPS.dat'
        
        orig_dir = os.getcwd()
        
        caps_dir = os.path.join(orig_dir,'workDir_tacs1/Scratch/tacs')
        
        bdf_dir = os.path.join(caps_dir,bdf_file)
        dat_dir = os.path.join(caps_dir,dat_file)
        
        bdf_newdir = os.path.join(orig_dir,bdf_file)
        data_newdir = os.path.join(orig_dir,dat_file)
        
        shutil.move(bdf_dir,bdf_newdir)
        shutil.move(dat_dir,data_newdir)
        
        #pytacs
        structOptions = {'writeSolution': True, }


        datFile = os.path.join(os.path.dirname(__file__), 'nastran_CAPS.dat')
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
            
        print(funcs)
        print(funcsSens)
        
        
        #might not have the node order right
        dfdX = funcsSens['load_set_001_wing_mass']['Xpts']
        #sens w.r.t. x1,y1,z1,x2,y2,z2,...
        #but 1,2,3,... might not match bdf nodes
        
        #reorder nodes, start at
        TACSnodeMap = FEASolver._getGlobalToLocalNodeIDDict()
        dfdX = dfdX.reshape(242,3)
        dfdX_bdf = np.zeros((242,3))
        for i in range(242):
            #i is bdf node, tacsNode is globalNode
            tacsNode = TACSnodeMap[i]
            dfdX_bdf[i,:] = dfdX[tacsNode,:]
            
        #now use dfdX_bdf for the sensitivity back to ESP/CAPS
if __name__ == '__main__':
    unittest.main()
