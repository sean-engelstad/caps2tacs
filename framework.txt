#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 15:03:32 2021

@author: sengelstad6
"""

#USES THE FOLLOWING:
#ESP/CAPS .csm file for desvar, geom, and generate .bdf
#pytacs to import .bdf file into tacs and use tacs
    #tacs, pyNASTRAN, baseclasses

# general optimization function:

# 1 - receive design variables, desvar

# 2 - build geometry from .csm file, desvar

# 3 - compute d(mesh)/d(desvar), md_grad from CAPS
    #not sure how to do this asked galbramc@mit.edu

# 4 - build .bdf file and mesh with NASTRAN AIM
    
# 5 - run pytacs to import .bdf file
    
# 6 - query the d(func)/d(mesh), the fm_grad from pytacs
    
# 7 - query the current function values from pytacs
     
# 8 - output function and gradient for ParOpt
    
#Port these functions objfunc, gradfunc into ParOpt
# RESULT: optimized desvar for the optimization problem