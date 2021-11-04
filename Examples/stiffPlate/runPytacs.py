#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 12:11:40 2021

@author: sengelstad6
"""

import os, shutil
from tacs.pytacs import pyTACS

structOptions = {'writeSolution': True, }

datFile = os.path.join(os.path.dirname(__file__), 'nastran_CAPS.dat')
# Load BDF file
FEASolver = pyTACS(datFile, options=structOptions)
# Set up TACS Assembler
FEASolver.createTACSAssembler()
# Read in forces from BDF and create tacs struct problems
SPs = FEASolver.createTACSProbsFromBDF()
print("ran pytacs")
# Solve each structural problem and write solutions
for caseID in SPs:
    FEASolver(SPs[caseID])
    FEASolver.writeSolution(outputDir=os.path.dirname(__file__))
    
